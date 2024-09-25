# This file contains functions for simulating disease transmission and formatting transmission model output


# Function returning results for states discretized by sim_parms$readdt for all pens given transmission results from a single iteration
state_list_names<-c("sus","latrec","latdead","infrec","infdead","recovered","dead","clinical","bloodpos","severeclin","deadinpen")



# Inputs: betawithin = within-pen transmission rate; betabetween = between-pen transmission rate;
#   betaroom = between room transmission rate; barnsize = number of pigs in the barn;
#   sim_parms = simulation parameters; betapeople = distance independent within barn transmission rate;
#   in_barn_pen_arraylist = list with arrays with 1) the room of each pen, 2) whether a pen is an edge pen or not,
#   3) which pens are adjacent, 4) which rooms are adjacent
# Outputs: formatted transmission model simulation run with a list per pen of the number of pigs over time in each disease state
get_single_sim <- function(betawithin, betabetween,betaroom, barnsize,sim_parms=sim_parms,betapeople=0,in_barn_pen_arraylist=in_barn_pen_arraylist){
  # array containing the number of pigs in each pen
  pensizearray<-rep(round(sim_parms$flock_size_pars/sim_parms$npens),sim_parms$npens)
  #infection starts in random pen
  sim_parms$init_pen<-sample(c(1:sim_parms$npens), size = 1)
  
  # determines shape and scale of gamma distribution for time to clinical signs for georgia strain
  if(sim_parms$is_georgia_strain){
    
    gamma_mean <- -1
    while(gamma_mean < sim_parms$georgia_clinmean_bounds[1] || gamma_mean > sim_parms$georgia_clinmean_bounds[2]){
      gamma_mean <- rnorm(1, mean = sim_parms$georgia_clinmean[1], sd = sim_parms$georgia_clinmean[2])
    }
    gamma_var <- 0.7358
    
    georgia_shape = (gamma_mean^2) / gamma_var
    georgia_scale = (gamma_var) / gamma_mean

    
    sim_parms$sero_pars = c(georgia_shape, georgia_scale)
  }

  # Calls C-function for simulating transmission model run
  return_val <- .C("MC_SEIR_parms", rv=as.integer(vector("integer", (sim_parms$ndays/sim_parms$time_step + 1)*11*sim_parms$npens)), n_runs=as.integer(sim_parms$n_runs), 
                   as.integer(sim_parms$npens),as.integer( pensizearray),mort_lat=as.double(sim_parms$mort_lat_pars), recov_lat=as.double(sim_parms$recov_lat),
                   mort_inf=as.double(sim_parms$mort_inf_pars), recov_inf=as.double(sim_parms$recov_inf), mort_comp_bounds=as.double(sim_parms$mort_par_bounds), 
                   recov_comp_bounds=as.double(sim_parms$recov_par_bounds), bpos_comp_bounds=as.double(sim_parms$bpos_par_bounds), 
                   sero_par=as.double(sim_parms$sero_pars), clin_par=as.double(sim_parms$clinpar),
                   bpos_par=as.double(sim_parms$bpospar), tsclin_par=as.double(sim_parms$tsclinpar), lsclin_par=as.double(sim_parms$lsclin_par),
                   mortdelay_par=as.double(sim_parms$mortdelaypar), morttrans_mult_in=as.double(sim_parms$morttrans_mult),
                   dt_in=as.double(sim_parms$time_step), num_dt_in=as.integer(sim_parms$ndays/sim_parms$time_step), 
                   p_mort_in=as.double(sim_parms$p_mort), p_sero_in=as.double(sim_parms$p_sero), init_inf_in=as.integer(sim_parms$init_inf),
                   initinf_pen=as.integer(sim_parms$init_pen), beta_in=as.double(betawithin),betabet_in=as.double(betabetween),betapeople_in=as.double(betapeople),
                   N0_in=as.integer(barnsize),nrooms=as.integer(sim_parms$nrooms),roomofpen=as.integer(in_barn_pen_arraylist$room_of_pen_arr),is_pen_edge=as.integer(in_barn_pen_arraylist$is_edge_pen_arr),
                   pen_distarr=as.double(in_barn_pen_arraylist$pen_distmult_linvec),room_distarr=as.double(in_barn_pen_arraylist$room_distmult_linvec),trans_type=as.integer(sim_parms$trans_type),betaroom_in=as.double(betaroom))$rv 
  
  # Processes C-function transmission run output
  return_val_processed <- process_single_sim(return_val,sim_parms =sim_parms  )
  
  return(return_val_processed)
}


# Inputs: single_sim_in = single iteration of transmission model output from C function; sim_parms = simulation parameters
# Outputs: List with the number of pigs in each disease state over time post virus exposure for each pen in the barn
process_single_sim <- function(single_sim_in,sim_parms=sim_parms){
  # define names for the different disease states
  statelist<-as.list(1:11)
  names(statelist)<-c("sus","latrec","latdead","infrec","infdead","recovered","dead","clinical","bloodpos","severeclin","deadinpen")
  
  
  # Inputs: inpen = numeric id of pen in the barn; instate = the disease state name; readdt = time step at which to save data
  # Output: matrix with the number of pigs over time in the given pen and disease state
  get_resultsfor_pen_state_readdt<-function(inpen=3,instate=statelist$sus,readdt=sim_parms$readdt){
    # starting index for selecting the data of interest from single_sim_in
    tstart<-(inpen-1)*(sim_parms$ndays/sim_parms$time_step+1)*11+ (sim_parms$ndays/sim_parms$time_step+1)*(instate-1)+1 
    # ending index for selecting the data of interest from single_sim_in
    tend<- (inpen-1)*(sim_parms$ndays/sim_parms$time_step+1)*11+ (sim_parms$ndays/sim_parms$time_step+1)*(instate) 
    
    # array of indexes for selecting data by desired time step (sim_parms$readdt)
    tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)+1
    # selects desired output from single_sim_in
    outmatrix<-single_sim_in[tstart:tend][tseq]
    
    return(outmatrix)
  }
  
  # Input: inpen = numeric id of pen in the barn
  # Output: List with the number of pigs in each disease state over time post virus exposure for the given pen 'inpen'. 
  get_resultsfor_onepen<-function(inpen=3){
    outlist<-lapply(1:length(statelist),function(g){get_resultsfor_pen_state_readdt(inpen,statelist[[g]])}) 
    names(outlist)<-names(statelist)
    return(outlist)
  }
  
  # List with the number of pigs in each disease state over time post virus exposure for each pen in the barn
  simresults_list<-lapply(1:sim_parms$npens,function(h){get_resultsfor_onepen(h)})
  
  return(simresults_list)
}



# Input: in_sim_resoneiter = list with the number of pigs in each disease state over time post virus exposure for 
#   each pen in the barn as formatted by the function process_single_sim
# Output: returns a list with the total number of pigs in each disease over time in the barn 
summarize_prem_oneiter<-function(in_sim_resoneiter){
  # Define disease state names
  statelist<-as.list(1:11)
  names(statelist)<-c("sus","latrec","latdead","infrec","infdead","recovered","dead","clinical","bloodpos","severeclin","deadinpen")
  
  # Input: stateno = iterator of a disease state
  # Output: 
  getonestate<-function(stateno=1){
    # extracts matrix with the number of pigs in a given disease state over time in each pen
    tout<-do.call("rbind",lapply(1:length(in_sim_resoneiter),function(g){in_sim_resoneiter[[g]][[stateno]]}))
    # gets total number of pigs in a given disease state over time in the barn
    tout<-apply(tout,2,sum)
    return(tout)
  }
  
  # Gets total number of pigs in each disease state over time in the barn
  out_list<-lapply(1:length(statelist),getonestate)
  names(out_list)<-names(statelist)
  return(out_list)
}


