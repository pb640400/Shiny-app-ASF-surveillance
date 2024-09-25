#include <R.h>
#include "Rmath.h"


/* Macros for implementing simple operations */
#define MAX(a,b) (((a)>(b))?(a):(b))  /* Maximum */
#define MIN(a,b) (((a)<(b))?(a):(b))
#define ZEROPAD(array_name,ind) ((ind)<0?(0):(array_name[ind]))  /* Produces 0 if requested array index is less than 0 */


#define probabilistic_trans 1  /* Use distributions (1) or just use mean values (0) */
#define initial_infectious 0   /* Start individual as infectious (1) or as latent (0) */
#define use_survival_period 0  /* Assume a distribution for the amount of time that birds live (use a survival/infected period as opposed to an infectious period) */


/* Define macros to make more readible code */
/* Used in 2D array, "state_array" */
/* States are:
0: Susceptible
1: Latent leading to recovery
2: Latent leading to mortality
3: Infectious leading to recovery
4: Infectious leading to mortality
5: Recovered
6: Mortality
7: Clinical signs
8: Blood positive
9: Severe clinical signs
10: Mortality still in the pen
*/
#define susceptible_st 0   /* Note that these states contain the current number of individuals */
#define latent_r_st 1
#define latent_m_st 2
#define infectious_r_st 3
#define infectious_m_st 4
#define recovered_st 5     /* Note that this state is cumulative (individuals do not leave) */
#define mortality_st 6     /* Note that this state is cumulative (individuals do not leave) */
#define sero_st 7
#define blood_st 8
#define severe_clin_st 9
#define mortdelay_st 10

/* Used for creating the 3D array, "state_array," to store aboved time states over time */
#define num_states 11

/* Macros for easily accessing state information from 3D array (number of individuals in each compartment in each pen), indexed by dt */
#define N_S(tpen,ind) (state_array[tpen][susceptible_st][ind])   /* Number of susceptible */
#define N_Er(tpen,ind) (state_array[tpen][latent_r_st][ind])     /* Number of latent individuals, eventually will be "recovered" */
#define N_Em(tpen,ind) (state_array[tpen][latent_m_st][ind])     /* Number of latent individuals, eventually will be mortalities */
#define N_E(tpen,ind) (N_Er(tpen,ind)+N_Em(tpen,ind))      /* Number of latent individuals */
#define N_Ir(tpen,ind) (state_array[tpen][infectious_r_st][ind])        /* Number of infected individuals (or infectious if the proper flag is set), eventually will be "recovered" */
#define N_Im(tpen,ind) (state_array[tpen][infectious_m_st][ind])  /* Number of infected individuals (or infectious if the proper flag is set), eventually will be "dead */
#define N_I(tpen,ind) (N_Ir(tpen,ind)+N_Im(tpen,ind))   /* Total number of infected individuals (or infectious if the proper flag is set) */
#define N_R(tpen,ind) (state_array[tpen][recovered_st][ind])    /* Number of "recovered" individuals. */
#define N_M(tpen,ind) (state_array[tpen][mortality_st][ind])    /* Number of mortalities. Since individuals do not leave, this is the cumulative number of mortalities. */
#define N_N(tpen,ind) (N_S(tpen,ind)+N_E(tpen,ind)+N_I(tpen,ind)+N_R(tpen,ind)) /* Total number of living individuals */
#define N_Ser(tpen,ind) (state_array[tpen][sero_st][ind]) /*Number with clinical signs*/
#define N_B(tpen,ind) (state_array[tpen][blood_st][ind]) /*Number blood positive */
#define N_SClin(tpen,ind) (state_array[tpen][severe_clin_st][ind]) /*Number with severe clinical signs*/
#define N_Md(tpen,ind) (state_array[tpen][mortdelay_st][ind]) /*Number of dead pigs not removed from pen*/

/* Macros producing pointer to beginning of a 2D array for each state (the names are the same as defined above) for a given pen */
/* Allows passing the arrays to functions */
#define array_N_S(tpen) state_array[tpen][susceptible_st]
#define array_N_Er(tpen) state_array[tpen][latent_r_st]
#define array_N_Em(tpen) state_array[tpen][latent_m_st]
#define array_N_Ir(tpen) state_array[tpen][infectious_r_st]
#define array_N_Im(tpen) state_array[tpen][infectious_m_st]
#define array_N_R(tpen) state_array[tpen][recovered_st]
#define array_N_M(tpen) state_array[tpen][mortality_st]
#define array_N_Ser(tpen) state_array[tpen][sero_st]
#define array_N_B(tpen) state_array[tpen][blood_st]
#define array_N_SClin(tpen) state_array[tpen][severe_clin_st]
#define array_N_Md(tpen) state_array[tpen][mortdelay_st]


/* Create structure to store information on gamma distributions */
typedef struct {
	/* Can directly simulate random variable */
	double shape;
	double scale;
	/* Allow for truncating, set to a negative number if not truncating at bound */
	double min;
	double max;
	/* Defined to allow for a faster look-up table approach */
	double * values;
	int * tindex;
} gamma_parm;

/* Create structure to store information on distributions for state transition times */
typedef struct {
	gamma_parm latent;
	gamma_parm infectious;
} cmpt_parms;

//Boolean type not native to C languange, must be user defined.
typedef int bool;
#define true 1
#define false 0

/* Global variable to allow swapping between look-up table or directly simulating random variables */
int USE_LUT;   /* TRUE to use look-up table (inverese CDF with discretization), FALSE to directly generate gamma distributions */




/* Function to perform disease state transitions for virus-exposed individuals */
/* Inputs:
-N: number of individuals to distribute
-parms: parameters for latent and infectious period distributions for individuals that die and recover
-s_parms: parameters for the distribution modeling the time to appearance of clinical signs
-c_parms: parameters for the distribution modeling the time clinical signs last
-b_parms: parameters for the distribution modeling the time to viremia
-t_to_sc_parms:
-legnth_sc_parms:
-md_parms:
-is_dying:
-p_sero:
-dt: Simulation discretization
-curr_ind: Denotes the time index that the simulation is currently on
-max_ind: end index for arrays storing state information so that C does not attempt to change variables off the end of the array
-lat_array[]: Array for latent state
-inf_array[]: Array for infectious state
-sero_array[]: Array for clinical signs state
-abs_array[]: Array for absorbing state (mortality or remove/recover)
-bpos_array[]: Array for viremic (blood positive) state
-sclin_array[]: array for severe clinical signs state
-mortdelay_array[]: array for mortality in a pen
-distrib_infectious: boolean for if newly infected birds start in latent or infectious period
Modifies: lat_array[], inf_array[], sero_array[], abs_array[], bpos_array[], sclin_array[], mortdelay_array[]
*/
void gamma_dist_event(int N, cmpt_parms parms, gamma_parm s_parms,gamma_parm c_parms, gamma_parm b_parms, gamma_parm t_to_sc_parms, gamma_parm length_sc_parms,
		gamma_parm md_parms, bool is_dying, double p_sero, double dt, int curr_ind, int max_ind, int lat_array[], int inf_array[], int sero_array[], int abs_array[],
		int bpos_array[], int sclin_array[], int mortdelay_array[], int distrib_infectious) {

	int ii, jj;  /* Loop variables */
	double t_latent, t_infectious, t_sero, t_clin, t_bloodpos, t_to_sclin, t_length_sclin, t_mortdelay;  /* Used to store the values from the pseudo-random number generator */
	int n_dt_latent, n_dt_infectious, n_dt_sero, n_dtint_clin, n_dt_bloodpos, n_dt_ttosclin, n_dt_lengthsclin, n_dt_mortdelay; /* Used to store the number of time steps (dt) for disease states*/

	/* Assumes that we can write to curr_ind+1 -- we checked before running function */
	/* Writing to array */
	lat_array[curr_ind] += N; /* All infected, should enter the latent period */


	for (ii = 0; ii < N; ii++) {  /* Perform disease state transitions for each new infection */

		/* This chunk randomly draws lengths for the infectious and latent periods and converts them into number of time steps */
		if (use_survival_period) {
			if (probabilistic_trans) {  /* Simulate with gamma distributions (or truncated gamma distributions) */
				do {
					t_latent = rgamma(parms.latent.shape, parms.latent.scale);  /* Use PRNG to generate amount of time in the latent state */
				} while (t_latent >= parms.latent.max && parms.latent.max > 0);  /*Truncating step. If latent period outside bounds, new number generated */
				do {
					t_infectious = rgamma(parms.infectious.shape, parms.infectious.scale);  /* Use PRNG to generate amount of time in the "infected" state */
				} while (t_infectious <= parms.infectious.min && parms.infectious.min > 0);  /* Truncating step for infectious period */
			}
			else {  /* Just use the mean values. Runs faster and closer to GLM/GEE approaches */
				t_latent = parms.latent.shape*parms.latent.scale;
				t_infectious = parms.infectious.shape*parms.infectious.scale;
			}
			/* Directly generate transitions from PRNGs */
			n_dt_latent = ((int)(MIN(t_latent, t_infectious) / dt + .5));  /* Round time to the nearest discretization */
			n_dt_infectious = ((int)(t_infectious / dt + .5));

			if (distrib_infectious) {  /* Special case for starting with an infectious bird */
				n_dt_latent = 0;
				n_dt_infectious = ((int)(MAX((t_infectious - t_latent), 0) / dt + .5));
			}
		}
		else { /* add latent and infectious periods */
			// This is the approach used in the Shiny app

			do {
				t_latent = rgamma(parms.latent.shape, parms.latent.scale);  /* Use PRNG to generate amount of time in the latent state */
			} while (t_latent > parms.latent.max && parms.latent.max > 0);  /*Truncating step. If latent period outside bounds, new number generated */
			do {
				t_infectious = rgamma(parms.infectious.shape, parms.infectious.scale);  /* Use PRNG to generate amount of time in the "infected" state */
			} while (t_infectious > parms.infectious.max && parms.infectious.max > 0);  /* Truncating step for infectious period */
			t_sero = rgamma(s_parms.shape, s_parms.scale); // generate time to clinical signs
            t_clin = rgamma(c_parms.shape, c_parms.scale); // generate length of clinical signs
            //generate time to blood pos relative to onset of infectious period
            do {
            	t_bloodpos = rnorm(b_parms.shape, b_parms.scale); //rnorm(mu, sigma);
            	//Truncating step. The condition (t_bloodpos + t_latent) < 0 is added to prevent a pig from being blood positive prior to exposure.
            } while ((t_bloodpos < b_parms.min) || (t_bloodpos > b_parms.max) || ((t_bloodpos + t_latent) < 0));
            //generate time to severe clinical time and length of severe clinical signs
            t_to_sclin = rgamma(t_to_sc_parms.shape, t_to_sc_parms.scale);
            t_length_sclin = rgamma(length_sc_parms.shape, length_sc_parms.scale);

			/* Convert to simulation times */
			n_dt_latent = ((int)(t_latent / dt + .5));  /* Round time to the nearest discretization */
			n_dt_infectious = ((int)((t_latent + t_infectious) / dt + .5));
			n_dt_sero = ((int)(t_sero / dt + .5)); //time until onset of mild clinical signs
			n_dt_bloodpos = ((int)(fabs(t_bloodpos) / dt + .5)); //fabs is absolute value for input of type double
			n_dt_ttosclin = ((int)(t_to_sclin / dt + 0.5)); //time until unset of severe clinical signs
			n_dtint_clin= ((int)((t_clin)/ dt + .5)); //duration of mild clinical signs
			n_dt_lengthsclin = ((int)(t_length_sclin / dt + .5)); //duration of severe clinical
			
			if (distrib_infectious) {  /* Special case for starting with an infectious bird */
				n_dt_latent = 0;
				n_dt_infectious = ((int)(t_infectious / dt + .5));
			}
		}

		// Update disease state arrays
		/* Check to make sure we are not already running off the simulation array */
		/* If so, do not need to calculate infectious to removed/mortality transition */
		if (n_dt_latent <= (max_ind - curr_ind)) {
			lat_array[curr_ind + n_dt_latent]--;  /* Denote the time leaving the latent compartment */
			inf_array[curr_ind + n_dt_latent]++;  /* Enter infectious compartment at the same time */

			if (n_dt_infectious <= (max_ind - curr_ind)) {
				inf_array[curr_ind + n_dt_infectious]--;  /* Denote the time leaving the infectious compartment */
				abs_array[curr_ind + n_dt_infectious]++;  /* Denote the time entering the absorbing compartment (mortality or recovered */


				// Dead pigs that remain in pen for some amount of time
				if(is_dying){
					t_mortdelay = runif(md_parms.shape, md_parms.scale);
					if(t_mortdelay > 0){

						mortdelay_array[curr_ind + n_dt_infectious]++;
						n_dt_mortdelay = ((int)((t_mortdelay + t_latent + t_infectious) / dt + .5));

						if(n_dt_mortdelay <= (max_ind - curr_ind)){
						mortdelay_array[curr_ind + n_dt_mortdelay]--;
						}
					}
				}
			}
		}


		// blood clinical signs code for time pig becomes blood positive in terms of the time prior to the onset of the
		// infectious period. Pig is removed from state at the end of the infectious period
		if (t_bloodpos < 0){ //pig becomes blood positive before onset of infectiousness
		// with this loop there is an else condition that we are missing potentially
		// what happens if there is a case where blood possitive occurs before sim time?
		
			if((curr_ind + n_dt_latent - n_dt_bloodpos) < 0){
				Rprintf("this loop shound not be running");			
			}
			if(((curr_ind + n_dt_latent - n_dt_bloodpos) >= 0) && ((curr_ind + n_dt_latent - n_dt_bloodpos) <= max_ind)){
				bpos_array[(curr_ind + n_dt_latent - n_dt_bloodpos)]++;
				// individual leaves blood positive state at the end of infectiousness
				if((curr_ind + n_dt_infectious) <= max_ind){
					bpos_array[(curr_ind + n_dt_infectious)]--;
				}
			}
		} else { //pig becomes blood positive at onset of infectiousness
			if (n_dt_latent <= (max_ind - curr_ind)) {
				bpos_array[(curr_ind + n_dt_latent )]++;
				// individual leaves blood positive state at the end of infectiousness
				if((curr_ind + n_dt_infectious) <= max_ind){
					bpos_array[(curr_ind + n_dt_infectious)]--;
				}
			}
		}

		// pigs that die have clinical signs until death; pigs that recover have clinical
		// signs for n_dtint_clin time periods. The earliest a pig can develop clinical signs is at the start of the infectious period
		if (runif(0, 1) < p_sero){ //check if individual develops clinical signs
			if (!is_dying){ //individuals to recover are handled differently from individuals to die
				if (n_dt_sero >= n_dt_latent){ // if clinical signs develop after the start of the infectious period
					if (n_dt_sero <= (max_ind - curr_ind)) {
						sero_array[curr_ind + n_dt_sero]++;
					}
					// removing after clinical period
					if ((n_dtint_clin+n_dt_sero) <= (max_ind - curr_ind)) {
						sero_array[curr_ind + (n_dtint_clin+n_dt_sero)]--;
					}
					
				}
				else {  // earliest clinical signs can develop is at the start of the infectious period
					if (n_dt_latent <= (max_ind - curr_ind)){
						sero_array[curr_ind + n_dt_latent]++;
						if ((n_dtint_clin+n_dt_latent) <= (max_ind - curr_ind)) {
						sero_array[curr_ind + (n_dtint_clin+n_dt_latent)]--;
						}
					}
				}
				
			}
			else { //pigs that die have clinical signs until the end of the infectious period (death)
				if (n_dt_sero >= n_dt_latent && n_dt_sero < n_dt_infectious){
					if (n_dt_sero <= (max_ind - curr_ind)) {
						sero_array[curr_ind + n_dt_sero]++;
					}
					if (n_dt_infectious <= (max_ind - curr_ind)){
						sero_array[curr_ind + n_dt_infectious]--;
					}
				}
				if (n_dt_sero < n_dt_latent) { // earliest clinical signs can develop is the start of the infectious period
					if (n_dt_latent <= (max_ind - curr_ind)){
						sero_array[curr_ind + n_dt_latent]++;
					}
					if (n_dt_infectious <= (max_ind - curr_ind)){
						sero_array[curr_ind + n_dt_infectious]--;
					}
				}
			}
		} //end if clinical signs develop

		//handles entry/exit into severe clinical signs state
		// If pig does not die, it enters into severe clincal signs states n_dt_ttosclin time steps after exposure
		// (unless n_dt_ttosclin occurs before transition out of the latent state, in which case the pig enters
		// at the time of transition out of the latent state) and remains there n_dt_lengthsclin time steps.
		// If a pig does die, it enters severe clinical signs state as above and exits at the end of the infectious
		// period.
		if (!is_dying){
			if (n_dt_ttosclin >= n_dt_latent){
				if (n_dt_ttosclin <= (max_ind - curr_ind)) {
					sclin_array[curr_ind + n_dt_ttosclin]++;
				}
				// removing after clinical period
				if ((n_dt_lengthsclin + n_dt_ttosclin) <= (max_ind - curr_ind)){
					sclin_array[curr_ind + (n_dt_lengthsclin + n_dt_ttosclin)]--;
				}
			}
			else {
				if (n_dt_latent <= (max_ind - curr_ind)){
					sclin_array[curr_ind + n_dt_latent]++;
					if ((n_dt_lengthsclin+n_dt_latent) <= (max_ind - curr_ind)) {
					sclin_array[curr_ind + (n_dt_lengthsclin+n_dt_latent) ]--;
					}
				}
			}

		}
		else { //pig dies
			if (n_dt_ttosclin >= n_dt_latent && n_dt_ttosclin < n_dt_infectious){
				if (n_dt_ttosclin  <= (max_ind - curr_ind)) {
					sclin_array[curr_ind + n_dt_ttosclin ]++;
				}
				if (n_dt_infectious <= (max_ind - curr_ind)){
					sclin_array[curr_ind + n_dt_infectious]--;
				}
			}
			if (n_dt_ttosclin < n_dt_latent) {
				if (n_dt_latent <= (max_ind - curr_ind)){
					sclin_array[curr_ind + n_dt_latent]++;
				}
				if (n_dt_infectious <= (max_ind - curr_ind)){
					sclin_array[curr_ind + n_dt_infectious]--;
				}
			}
		}


	}//end for loop
}



/* Function to run a single SEIR simulation */
/* Inputs:
*state_array[]: array of pointers to rv[], creating a 2D array for storing state trajectories
*cpen_distarray[]: array of pointers to pen_distarr for accessing which pens are adjacent
*croom_distarray[]: for accessing which rooms are adjacent/transmission partners
beta = within-pen transmission rate
betabetween = between pen transmission rate
betapeople = within room distance independent transmission
mortp = shapes and scales of gamma distributed latent and infectious periods for pigs that die following infection
recovp = shapes and scales of gamma distributed latent and infectious periods for pigs that recover following infection
serop = distribution parameters for time to clinical signs
clinp = distribution parameters for length of clinical signs
bposp = distribtion parameters for when pigs are viremic relative to onset of infectious period
tsclinp = distribution parameters for time to severe clinical signs
lsclinp = distribution parameters for length of severe clinical signs
mortdelayp = distribution parameters for how long mortality may remain in a pen
morttrans_mult = scalar for increasing or decreasing transmissiblity of mortality
dt = day time step of simulation
num_dt = total number of time steps
num_pens = number of pens
p_mort = probability of mortality following infection
p_sero = probability of developing clinical signs following infection
N0 = total number of pigs
init_cond[] = array with the number of susceptible pigs in the initially exposed pen, number of initially exposed pigs, and the pen number of the initially exposed pen
pensizearr[] = number of pigs per pen
trans_type = transmission scenario
betaroom = between room tranmsission rate
roomofpen[] = room number for each pen
is_pen_edge[] = indicator for if a pen is an edge pen
nrooms = number of rooms in the barn
Modifies: rv[]. (*state_array[] is an array of points that point to the originally passed "rv[]" 1D arrays, allowing them to be treated as 2D arrays)
*/
void MC_SEIR_run(int *state_array[][num_states],double *cpen_distarray[],double *croom_distarray[], double beta, double betabetween,double betapeople,
		cmpt_parms mortp, cmpt_parms recovp, gamma_parm serop,gamma_parm clinp, gamma_parm bposp, gamma_parm tsclinp, gamma_parm lsclinp, gamma_parm mortdelayp,
		double morttrans_mult, double dt, int num_dt, int num_pens, double p_mort, double p_sero, int N0, int init_cond[], int pensizearr[],int trans_type,
		double betaroom,int roomofpen[],int is_pen_edge[],int nrooms) {

	// Initalized variable include number latently infected, number latently infected to die, number latently infected to recover, time iterator, number initially infected to die, iterators
	int S_E, S_Em, S_Er, tt, dN_Sm,ii,pp,kk,rr;
	// Total number of infectious pigs in the barn
	int Itot;
	// Total number of alive pigs in the barn
	int Ntot;
	// Probability of a susceptible pig becoming infected
	double binprob;
	// Total number of infectious pigs in each room
    int Itot_room[nrooms];
    // Total number of alive pigs in each room
	int Ntot_room[nrooms];
	// Stores transmission probabilities
	double temproabA,temproabB,temproabC,temproabD;


	// Number of susceptibles animals in initially infected pen
	N_S((init_cond[2]-1),0) = init_cond[0];
	// loops through pens and initalizes the number of susceptible pigs in all other pens than the initially infected pen
	for (ii = 0; ii< num_pens; ii++){
		if (ii!=(init_cond[2]-1)){
		N_S(ii,0) = pensizearr[ii];
		}
	}
	
	// /* Number initially infected, split into recovery/mortality
	 dN_Sm = rbinom(init_cond[1], p_mort);


	// Determines disease state transitions for initially infected pigs that die
	gamma_dist_event((dN_Sm), mortp, serop,clinp, bposp, tsclinp, lsclinp, mortdelayp, true, p_sero, dt, 0, num_dt, array_N_Em(init_cond[2]-1), array_N_Im(init_cond[2]-1), array_N_Ser(init_cond[2]-1), array_N_M(init_cond[2]-1), array_N_B(init_cond[2]-1), array_N_SClin(init_cond[2]-1), array_N_Md(init_cond[2]-1), initial_infectious);
	// Determines disease state transitions for initially infected pigs that recover
	gamma_dist_event((init_cond[1] - dN_Sm), recovp, serop,clinp, bposp, tsclinp, lsclinp, mortdelayp, false, p_sero, dt, 0, num_dt, array_N_Er(init_cond[2]-1), array_N_Ir(init_cond[2]-1), array_N_Ser(init_cond[2]-1), array_N_R(init_cond[2]-1), array_N_B(init_cond[2]-1), array_N_SClin(init_cond[2]-1), array_N_Md(init_cond[2]-1), initial_infectious);
    
	
	
	// /* Main loop to move through time */
	for (tt = 1; tt <= num_dt; tt++) {

		Itot=0;
		Ntot=0;
		for (rr = 0; rr< nrooms; rr++){
			Itot_room[rr]=0;
			Ntot_room[rr]=0;
		}


		// Update the number of infectious and alive pigs in the barn and in each room at time tt
		for (kk = 0; kk< num_pens; kk++){
			Itot+=N_I(kk,(tt-1));
			Ntot+=N_N(kk,(tt-1));
			Itot_room[roomofpen[kk]-1]+=N_I(kk,(tt-1));
			Ntot_room[roomofpen[kk]-1]+=N_N(kk,(tt-1));
		}


		// loop through the pens
		for (pp = 0; pp< num_pens; pp++){

			if (N_N(pp,tt - 1) > 0){ 	// prevents division by 0
				double temp_betaroom=0;
				int hh=0;
				double temp_denom_sroom=0;
				double tempbetabetween=0;

				// determines tranmission probabilities based on selected transmission scenario
				switch(trans_type) {

				case 1  : //similar to guinat na+nb in denominator of force of infection
		
					temproabA=0;
					temproabB=0;
	
					temproabA=-1 * beta*((double)(N_I(pp,tt - 1)) / ((double)N_N(pp,tt - 1)))*dt;
					// Number of infectious pigs in pens other than pp divided by the total number of pigs in all pens
					temproabB=-1 * betabetween*((double)(Itot-N_I(pp,tt - 1))/(double)Ntot)*dt;
					binprob = (1 - exp(temproabA+temproabB));
					break;


				case 2  :// this is with only nb in denominator similar to royal vetmed
		
					if ((Ntot-N_N(pp,(tt-1))) > 0){ //are there any animals in other pens to avoid division by zero
						temproabA=0;
						temproabB=0;
			
						temproabA=-1 * beta*((double)(N_I(pp,tt - 1)) / ((double)N_N(pp,tt - 1)))*dt;
						// Number of infectious pigs in other pens other than pp divided by the total number of pigs in other pens
						temproabB=-1 * betabetween*((double)(Itot-N_I(pp,tt - 1))/(double)(Ntot-N_N(pp,(tt-1))))*dt;
						binprob = (1 - exp(temproabA+temproabB));
  
					} else {
						binprob = 0;
					}

					break;
		
		
				case 3 : //hohle density just based on number of infectious
		
					temproabA=0;
					temproabB=0;
			
					temproabA=-1 * beta*((double)(N_I(pp,tt - 1)) )*dt;
					// force of infection coming from all other pens. Frequency dependent transmission
					temproabB=-1 * betabetween*((double)(Itot-N_I(pp,tt - 1)))*dt;
					binprob = (1 - exp(temproabA+temproabB));
      			
					break;
		
		
				case 4 : //rooms and stuff and as implemented by ausvet but without pen level kernel. Need to verify logic regarding denominator
		
					temproabA=0;
					temproabB=0;
					temproabA=-1 * beta*((double)(N_I(pp,tt - 1)) / ((double)N_N(pp,tt - 1)))*dt;
					temp_betaroom=0;
		
					for (hh = 0; hh< num_pens; hh++){//loop to sum up source of infection terms for different source pens
						temp_betaroom=0;
						// using between pen beta if same room or between room beta otherwise
						if(roomofpen[hh]!=roomofpen[pp]){
							temp_betaroom=betaroom;
							temproabB+=-1 * temp_betaroom*((double)(N_I(hh,tt - 1))/(double)Ntot)*dt;//using source pen betas in the different rooms force of infection
						}else{
							temp_betaroom=betabetween;
							temproabB+=-1 * temp_betaroom*((double)(N_I(hh,tt - 1))/(double)Ntot)*dt;//using source pen betas in the same room force of infection
						}
					} // end hh loop for source pen
		
					binprob = (1 - exp(temproabA+temproabB));
		
					break;


				case 5 :
					// perez and klinkenberg email
					//Within-room FOI: βwSroom t (Is room t + Ic room t)/Nroom t

					//Between-room FOI: βbSherd t (Is herd t + Ic herd t )/Nherd t
					//interpreting pen as a room
					temproabA=0;
					temproabB=0;
	
		
					temproabA=-1 * beta*((double)(N_I(pp,tt - 1)) / ((double)N_N(pp,tt - 1)))*dt;
					// Number of infectious birds in pens other than pp divided by the total number of birds in all pens
					// having itot at the top is the main kernel that is used in perezes work
					// this is also what Dr Klinkenberg said with between pen contacts representing a different contact mechanism such as aerosol
					temproabB=-1 * betabetween*((double)(Itot)/(double)Ntot)*dt;
					binprob = (1 - exp(temproabA+temproabB));
		
					break;


				case 6 :
					// Nielsen with scaling factor of 1 for adjecent pens and zero otherwise with one room seems best.
					// Nielson between pen term is sum (Iadjacent)/Ntotal when rho is 1 for adjacent pens and 0 otherwise.
					// Note this can also accomodate other nielson kernals if we give the correct pen distance matrix in R
					temproabA=0;
					temproabB=0;
					// within-pen transmission
					temproabA=-1 * beta*((double)(N_I(pp,tt - 1)) / ((double)N_N(pp,tt - 1)))*dt;

					//the contacts are only divided among adjacent pens and not the whole premises. So the denominator should be the number of
					//animals in adjacent pens. This is a variable to capture the number of animals in the adjacent pens in the same room
					temp_denom_sroom=0;
					for (hh = 0; hh< num_pens; hh++){//loop to sum up source of infection terms for different source pens

						temp_betaroom=0;

						// using between pen beta if same room or between room beta otherwise
						if(roomofpen[hh]!=roomofpen[pp]){
							temp_betaroom=betaroom;
							temproabB+=-1 * temp_betaroom*((double)(N_I(hh,tt - 1))/(double)Ntot)*dt;
						}else{
							// adding up denominator for adjacent pens in same room. i.e contacts are divided among the available adjacent pens.
							temp_denom_sroom+=cpen_distarray[pp][hh]*(double)N_N(hh,tt - 1);
							temp_betaroom=betabetween;
							temproabB +=-1 * temp_betaroom*cpen_distarray[pp][hh]*((double)(N_I(hh,tt - 1)))*dt;
						}
		
					} // end hh loop for source pen
					// adding the recipient pen term to be consistent with nielson
					temp_denom_sroom+=(double)(N_N(pp,tt - 1));
					temproabB=temproabB/temp_denom_sroom;
					binprob = (1 - exp(temproabA+temproabB));

					break;


				case 7:
					//This is the default transmission scenario used for ASFV used in the Shiny app. For details on formulation see:
					// Ssematimba A, Malladi S, Bonney PJ, Charles KM, Boyer TC, Goldsmith T, Cardona CJ, Corzo CA, Culhane MR. African swine fever detection and
					// transmission estimates using homogeneous versus heterogeneous model formulation in stochastic simulations within pig premises. Open Veterinary
					// Journal. 2022;12(6):787-96.
		
					// within-pen force of infection
					temproabA=0;
					// between pen direct contact force of infection
					temproabB=0;
					// between pen indirect contact force of infection
					temproabC=0;
					// between room force of infection
					temproabD=0;

					// force of infection for within-pen transmission
					temproabA=-1 * beta*( ((double)(N_I(pp,tt - 1)) + morttrans_mult*(double)(N_Md(pp,tt-1))) / ((double)(N_N(pp,tt - 1)) + (double)(N_Md(pp,tt-1))) )*dt;

					//tracks the number of alive pigs in pens adjacent to pp to divide by the proper term for direct contact between pen force of infection
					temp_denom_sroom=0;
		
		
					for (hh = 0; hh< num_pens; hh++){//loop to sum up source of infection terms for different source pens
						temp_betaroom=0;
						tempbetabetween=0;

						// using between pen beta if same room or between room beta otherwise
						if(roomofpen[hh]!=roomofpen[pp]){ //pen hh and pp are in different rooms
			
							temp_betaroom=betaroom;
							if((Ntot-Ntot_room[roomofpen[pp]-1]) > 0){ //makes sure no division by 0
								// between room force of infection
								temproabD+=-1 * temp_betaroom*((double)(N_I(hh,tt - 1))/(double)(Ntot-Ntot_room[roomofpen[pp]-1]))*dt;
							}

						}else{ //pen hh and pp are in the same room

							// adding up denominator for adjacent pens in same room, i.e, contacts are divided among the available adjacent pens.
							//Cpendist array will be 1 if pens are adjacent. Therefore just adding the terms for adjacent pens. I.e., the nose to nose contacts
							//are divided between the adjacent pens (force of infection for pen B is Betanose*(IA+IC)/(NA+NB))
							temp_denom_sroom+=cpen_distarray[pp][hh]*(double)N_N(hh,tt - 1);

							// between pen direct contact transmission rate
							tempbetabetween=betabetween;
			
							//at this stage temp proabB is just the numerator
							// if the recipent pen is edge beta between is reduced by half. One wall=half transmission!
							if(is_pen_edge[pp]==0){
								temproabB +=-1 * tempbetabetween*cpen_distarray[pp][hh]*((double)(N_I(hh,tt - 1))) *dt;
							} else {
								temproabB +=-1 * tempbetabetween*cpen_distarray[pp][hh]*((double)(N_I(hh,tt - 1)))*dt/2;
							}
			
							// between pen fomite-mediated transmission
							if((double)(Ntot_room[roomofpen[pp]-1])>0){

								// contribution of hh to the force of infection depends on number of infectious pigs in pen hh / total number of pigs in the room
								temproabC +=-1 * betapeople*((double)(N_I(hh,tt - 1))/(double)(Ntot_room[roomofpen[pp]-1]))*dt;

							}
						}
		
					} // end hh loop for source pen

					//Dividing by denominator for adjacent pen nose-to-nose force of infection
					if(temp_denom_sroom>0){
						temproabB=temproabB/temp_denom_sroom;
					}

					// total probability of infection
					// temproabD is between room, temproabB is adjacent pen, temproabA is within herd and temproabC is distance independent pathways related to beta people.
					binprob = (1 - exp(temproabA+temproabB+temproabC+temproabD));
		
					break;


				default :
					binprob = 0;
				}//switch end
		
			
				// Simulate the number of susceptible pigs in pen pp infected in time period tt-1
				S_E = rbinom(N_S(pp,tt - 1), binprob);

			}// end if for whether there are any animals in the pen
			else {
				S_E = 0;
			}
		

		// /* split the newly infected into recover/mortality */
		dN_Sm = rbinom(S_E, p_mort);


		 if (tt < num_dt) {
			// /* disease state transitions for eventual mortalities */
			 gamma_dist_event((dN_Sm), mortp, serop,clinp, bposp, tsclinp, lsclinp, mortdelayp, true, p_sero, dt, tt, num_dt, array_N_Em(pp), array_N_Im(pp), array_N_Ser(pp), array_N_M(pp), array_N_B(pp), array_N_SClin(pp), array_N_Md(pp), 0);

			// /* disease state transitions for eventual recoveries */
		     gamma_dist_event((S_E - dN_Sm), recovp, serop,clinp, bposp, tsclinp, lsclinp, mortdelayp, false, p_sero, dt, tt, num_dt, array_N_Er(pp), array_N_Ir(pp), array_N_Ser(pp), array_N_R(pp), array_N_B(pp), array_N_SClin(pp), array_N_Md(pp), 0);
		 }

		 // /* Update states */
		//Note: These are macros defined at the beginning of the program. Number of birds in each state per time period stored in 3-dimensional "state_array"
		 N_S(pp,tt) = N_S(pp,tt - 1) - S_E;
		 N_Er(pp,tt) += N_Er(pp,tt - 1);
		 N_Em(pp,tt) += N_Em(pp,tt - 1);
		 N_Ir(pp,tt) += N_Ir(pp,tt - 1);
		 N_Im(pp,tt) += N_Im(pp,tt - 1);
		 N_R(pp,tt) += N_R(pp,tt - 1);
		 N_M(pp,tt) += N_M(pp,tt - 1);
		 N_Ser(pp,tt) += N_Ser(pp,tt - 1);
		 N_B(pp,tt) += N_B(pp,tt-1);
		 N_SClin(pp,tt) += N_SClin(pp,tt-1);
		 N_Md(pp,tt) += N_Md(pp,tt-1);
		} //end for loop through pens
	 } //end for loop through time

}






/* External function for R */
/* Run transmission simulation */
/* Inputs:
rv[]: "Return vector", for state trajectories. Transformed into a 3D array in C and a 1D vector in R
*n_runs: number of Monte Carlo simulations to perform
*n_pens: total number of pens
mort_lat[]: shape and scale of gamma distributed latent period for pigs to die
recov_lat[]: shape and scale of gamma distributed latent period for pigs that recover
mort_in[]: shape and scale of gamma distributed infectious period for pigs that die
recov_inf[]: shape and scale of gamma distributed infectious period for pigs that recover
mort_comp_bounds[]: min/max bounds for latent and infectious periods of pigs that die
recov_comp_bounds[]: min/max bounds for latent and infectious periods of pigs that recover
bpos_comp_bounds[]: min/max bounds of when pigs are viremic relative to infectiousness
sero_par[]: shape and scale of gamma distributed time to appearance of clinical signs following infection
clin_par[]: shape and scale of gamma distributed length of clinical signs
bpos_par[]: mean and standard deviation of normal distributed time to viremia relative to start of the infectious period
tsclin_par[]: shape and scale of gamma distributed time to appearance of severe clinical signs following infection
lsclin_par[]: shape and scale of gamma distributed length of severe clinical signs
mortdelay_par[]: min and max of uniform distribution for how long pigs remaining in a pen following death
*morttrans_mult_in: scaling factor increasing or decreasing contribution of dead pigs left in the pen to transmission of virus
*dt_in: day time step
*num_dt_in: total number of simulation time steps
*p_mort_in: probability a pig dies following infection
*p_sero_in: probability a pig develops clinical signs following infection
*init_inf_in: number of initially infected pigs
*init_pen: pen location of initially infected pigs
*beta_in: within-pen transmission rate
*betabet_in: between pen transmission rate
*betapeople_in: distance independent transmission rate
*N0_in: total number of pigs
*nrooms: total number of rooms
roomofpen[]: room number of each pen
is_pen_edge[]: indicator for if a pen is on the edge of a row
pendistarr[]: indicator for if pens are adjacent
room_distarr[]: indicator for if rooms are adjacent/transmission partners
*trans_type: transmission scenario
*betaroom_in: between-room transmission rate
Modifies:
rv[]: used to return state trajectories
*/
void MC_SEIR_parms(int rv[], int *n_runs,int *npens,int pensizearr[], double mort_lat[], double recov_lat[], double mort_inf[], double recov_inf[], double mort_comp_bounds[],
		double recov_comp_bounds[], double bpos_comp_bounds[], double sero_par[], double clin_par[], double bpos_par[], double tsclin_par[], double lsclin_par[], double mortdelay_par[], double *morttrans_mult_in, double *dt_in, int *num_dt_in, double *p_mort_in, double *p_sero_in,
		int *init_inf_in, int *init_pen, double *beta_in, double *betabet_in,  double *betapeople_in,int *N0_in,int *nrooms,int roomofpen[],int is_pen_edge[], double pen_distarr[],double room_distarr[],int *trans_type,double *betaroom_in) {

	USE_LUT = 0; /* Not using the look-up table, directly generating random distributions */

	/* Assign values to structures storing information on distributions for the latent period, the infected/survival period (or infectious if the flag is set), etc. */
	/* Additionally, while it has not been tested, there are components of the structure for truncating the distributions. If using this feature, would need to make sure that (or correct for) it does not significantly change the mean of the distribution */
	/* cmpt_parms is a user defined strutcure. See program header */
	cmpt_parms recovp, mortp;  /* Initialize structures for the parameters for the distributions for individuals that will eventually "recover" (or as modeled, become closer to persistently infectious) and ones that will become mortalities  */
	mortp.latent.shape = mort_lat[0];
	mortp.latent.scale = mort_lat[1];
	mortp.infectious.shape = mort_inf[0];
	mortp.infectious.scale = mort_inf[1];
	recovp.latent.shape = recov_lat[0];
	recovp.latent.scale = recov_lat[1];
	recovp.infectious.shape = recov_inf[0];
	recovp.infectious.scale = recov_inf[1];
	mortp.latent.max = mort_comp_bounds[1]; 
	mortp.infectious.max = mort_comp_bounds[3];
	recovp.latent.max = recov_comp_bounds[1];
	recovp.infectious.max = recov_comp_bounds[3];
  
	gamma_parm serop; //Initialize structure for time to clinical signs
	serop.shape = sero_par[0];
	serop.scale = sero_par[1];
	gamma_parm clinp; //Initialize structure for ending clinical distribution
	clinp.shape = clin_par[0];
	clinp.scale = clin_par[1];
	gamma_parm bposp; // initialize structure for when pig is viremic
	bposp.shape=bpos_par[0];
	bposp.scale=bpos_par[1];
	bposp.min=bpos_comp_bounds[0];
	bposp.max=bpos_comp_bounds[1];
	gamma_parm tsclinp, lsclinp; // initialize structures for time to and length of severe clinical signs
	tsclinp.shape = tsclin_par[0];
	tsclinp.scale = tsclin_par[1];
	lsclinp.shape = lsclin_par[0];
	lsclinp.scale = lsclin_par[1];
	gamma_parm mortdelayp; // initialize structure for how long mortality may remain in a pen
	mortdelayp.shape = mortdelay_par[0];
	mortdelayp.scale = mortdelay_par[1];

	/* Copy other inputted constants (*beta_in, etc. are inputs to the function) */
	double dt = *dt_in;        /* Time discretization, typically would be a fraction of a day */
	double p_mort = *p_mort_in;   /* The mortality fraction; the probability of mortality vs. recovery of infected pigs */
	double p_sero = *p_sero_in; //The probability of pigs developing clinical signs post infection
	double morttrans_mult = *morttrans_mult_in; //Parameter that can increase or decrease the transmission potential of dead pigs not yet removed from pens
	// transmission parameters
	double beta = *beta_in;
	double betabetween = *betabet_in;
	double betapeople = *betapeople_in;
	double betaroom = *betaroom_in;
	int N0 = *N0_in; //total number of pigs

	int init_cond[3]; //Initializes array for 1) the number of susceptible pigs in the initially exposed pen; 2) the number of initially exposed pigs; 3) the initially exposed pen number
	

	//getting room and pen distances in more usable data structure
	double *cpen_distarray[*npens];
	int  iii, jjj;
	int kk;
	//asigning the pendistance 2d array in C
	for (kk = 0; kk < *npens; kk++) {
		cpen_distarray[kk] = pen_distarr + kk *(*npens);
	}
	double *croom_distarray[*nrooms];
	for (kk = 0; kk < *nrooms; kk++) {
		croom_distarray[kk] = room_distarr + kk*(*nrooms);
	}
	
	/* Create a 2D array for storing information about states */
	/* Note that macros were defined at the beginning of the file for easy access to this 2D array */
	/* Also note that this array was defined in R and is accessible to R once the file is finished running if one wants a single simulation (the array is placed back into a matrix in R--it just is easier to pass a matrix as 1D array) */
	int *state_array[*npens][num_states];
	int ii,jj;
	for (ii = 0; ii < *npens; ii++) {
		for (jj = 0; jj < num_states; jj++) {
			// A call to rv is the memory location of 0th element in the array.
			// Adding some term "k" to rv selects the address of the "k"th element of rv.
			// (Aside: if we wanted the value stored at that location we would use *(rv+k). rv acts as a pointer to the 0th element)
			// So the below has pointers in the 2-D pointer array pointing to different memory locations of rv.
			state_array[ii][jj] = rv + ii*(*num_dt_in+1)*(num_states)+jj*(*num_dt_in + 1);
		}
	}

	/* Get the state of the pseudo-random number generator (PRNG) from R */
	GetRNGstate();

	/* For each simulation iteration */
	for (ii = 0; ii<*n_runs; ii++) {

		init_cond[0] = pensizearr[*init_pen-1] - *init_inf_in;
		init_cond[1] = *init_inf_in;
        init_cond[2] = *init_pen;

		memset(rv, 0, sizeof(int)*(*num_dt_in + 1)*(num_states)*(*npens));  /* Zero memory in order to perform a new simulation */
		// simulates disease transmission for a single iteration
		MC_SEIR_run(state_array,cpen_distarray,croom_distarray, beta,betabetween, betapeople, mortp, recovp, serop,clinp, bposp, tsclinp, lsclinp, mortdelayp, morttrans_mult, dt, *num_dt_in, *npens, *p_mort_in, *p_sero_in, N0, init_cond,pensizearr,*trans_type,betaroom,roomofpen,is_pen_edge,*nrooms);  /* Run a simulation */

	}
	
	PutRNGstate();  /* Pass the state of the PRNG back to R */
}


