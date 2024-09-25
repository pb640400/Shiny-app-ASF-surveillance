# This file contains functions for simulating routine mortality and morbidity as well as simulating sampling and testing
# of blood and/or tissue samples from sick/apparently healthy or dead pigs, respectively

library(extraDistr) #for rmultihyper distribution



# This function simulates routine mortality for each pen, day post virus exposure run in the transmission model,
# and simulation iteration. Has form normmort[[x]][y,z] where x is the iteration, y is the pen, and z is the day.
# !!! Note: the function provided here is an approximate version of the normal mortality function used in the app
# which was not shared due to its reliance on actual mortality data from swine producers and our inability to share this
# data to maintain confidentiality !!!
# Inputs: sim_parms = simulation parameters; pensizearray = number of pigs in each pen
# Outputs: list with the simulated routine mortality for each pen on each simulation day for each simulation iteration
sim_normmort_func <- function(sim_parms, pensizearray){
  
  # Checks whether to use normal or high routine mortality rate
  if(sim_parms$use_high_mortdata){
    
    sim_pennormmort_approx <- function(){
      
      # simulates mortality on each simulation day for each pen in the barn
      temp_penmort <- sapply(1:sim_parms$npen, function(x){
        # simulates a high daily mortality rate for routine mortality at the pen level for each simulation day
        tempmortrate <- rexp(n = (sim_parms$ndays + 1), rate = 1/0.004) + 0.010
        # converts the mortality rate into an actual number of dead pigs
        tempmort <- round(pensizearray[x] * tempmortrate)
        return(tempmort)
      })
      
      return(t(temp_penmort))
    }
    
  }else{
    
    sim_pennormmort_approx <- function(){
      
      # simulates mortality on each simulation day for each pen in the barn
      temp_penmort <- sapply(1:sim_parms$npen, function(x){
        # simulates an average daily mortality rate for routine mortality at the pen level for each simulation day
        tempmortrate <- rexp(n = (sim_parms$ndays + 1), rate = 1/0.003)
        # converts the mortality rate into an actual number of dead pigs
        tempmort <- round(pensizearray[x] * tempmortrate)
        return(tempmort)
      })
      
      return(t(temp_penmort))
    }
    
  }
  
  # simulates pen-level mortality across the simulation iterations
  all_normmort <- lapply(1:sim_parms$num_iterations, function(y){sim_pennormmort_approx()})
  
  return(all_normmort)
}



# Function for simulating the number of pigs with clinical signs due to non-ASF related causes (normal morbidity) for
#   each pen at each simulation time point and for each simulation iteration
# Inputs: transmission_results = number of pigs in each disease state over time post virus exposure in each pen for each
#   simulation iteration; in_sim_parms = simulation parameters; pensizearray = number of pigs per pen
# Outputs: list with the simulated routine morbidity with form all_normsick[[x]][y,z] where x is iteration, y is pen, z is simulation timepoint
sim_normsick_func <- function(transmission_results=transmission_results, in_sim_parms=in_sim_parms,pensizearray){
  
  # Simulates the routine morbidity for a particular sim time and iteration and pen
  # Inputs: time_test_readtindex = simulation time; initer = iteration number; inpen = pen number; transmission_results =
  #   pen-level transmission model output
  # Outputs: the number of pigs with clincal signs due to non-ASF related causes in the given pen, at the given time point,
  #   and for the given iteration number
  sim_singlenormsick <- function(time_test_readtindex = 110, initer = 5, inpen = 5,transmission_results=transmission_results){
    # Number of pigs alive at simulation time point
    temp_alive <- pensizearray[inpen]-transmission_results[[initer]][[inpen]][[7]][time_test_readtindex]
    # Number of infectious pigs at the simulation time point
    temp_inf <- transmission_results[[initer]][[inpen]]$infrec[time_test_readtindex] + transmission_results[[initer]][[inpen]]$infdead[time_test_readtindex]
    # Number of non-infectious pigs still alive
    alive_nonsick <- temp_alive - temp_inf
    
    # Simulates the number of the living, non-infectious pigs that have clinical signs due to routine causes
    if(alive_nonsick > 0){
      # simulates the probability of having clinical signs due to routine causes
      tproab=rspert(n = 1,x.min = in_sim_parms$normsickparms[1], x.max=in_sim_parms$normsickparms[3],x.mode = in_sim_parms$normsickparms[2])
      # simulates the number of living-noninfectious pigs with clinical signs due to routine causes based on tproab
      out<-rbinom(n=1,size = alive_nonsick,prob =  tproab)
      return(out)
    } else {return(0)}
  }
  
  
  # simulates the routine morbidity for each pen, simulation time point, and simulation iteration
  # all_normsick[[x]][y,z] --> x is iteration, y is pen, z is simulation timepoint
  all_normsick <- lapply(1:in_sim_parms$num_iterations, function(x){
    pen_normsick <- sapply(1:in_sim_parms$npens, function(y){
      t_normsick <- sapply(1:(in_sim_parms$ndays / in_sim_parms$readdt + 1), function(z){
        sim_singlenormsick(time_test_readtindex = z, initer = x, inpen = y, transmission_results = transmission_results)
      })
    })
    return(t(pen_normsick))
  })
  
  return(all_normsick)
}



# This function simulates blood (or oral) and/or tissue sampling, sample pooling, and testing by PCR for a given simulation iteration
# and time post virus exposure when sampling occurs
# Inputs: time_test_readtindex = simulation time index; initer = iteration number; sample_type = number indicating the
#   sampling prioritization scheme for dead/tissue sampling (details given below); insamplesize = total number of samples
#   to take; transmission_results = pen-level simulated transmission model output; sim_normmort = simulated routine mortality; 
#   sim_normsick = simulated routine morbidity; transmission_results_premlevel = premises-level simulated transmission
#   model output; boolretnumsamp = T/F for returning number of samples taken from sick, dead and live pigs; 
#   in_selectedpensvec = allows user to specify pens from which samples are taken; pensizearry = array with the number
#   of pigs in each pen; in_sim_parms = simulation parameters
# Outputs: If boolretnumsamp == TRUE, returns 1 if detected or 0 if not, the number of pooled samples tested by PCR,
#   and the number of samples taken from sick, dead, and live/apparently healthy pigs. If boolretnumsamp == FALSE,
#   only returns 1 or 0 if detected or not

# sample_type = 1: dead, sick, live sampling; sample_type = 2: sick, dead, live sampling; sample_type = 3: dead then live sampling;
# sample_type = 4: dead only; sample size =5: dead then sick sampling; sample type =6: sick then non sick live, sample_type = 7:
# random live sampling
sample_blood_iter_simtime<-function(time_test_readtindex = 49, initer = 5, 
                                    sample_type=1,insamplesize,transmission_results=transmission_results,
                                    sim_normmort=sim_normmort,
                                    sim_normsick=sim_normsick,
                                    transmission_results_premlevel=transmission_results_premlevel,
                                    boolretnumsamp=FALSE,
                                    in_selectedpensvec=NA,
                                    pensizearray=pensizearray, in_sim_parms){
  
  # Get number of dead pigs in each pen due to ASF at given simulation time point
  pendismort <- sapply(1:in_sim_parms$npens, function(x){get_pen_dismort(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results,in_sim_parms = in_sim_parms)})
  # Get number of sick pigs in each pen due to ASF at given simulation time point
  pendissick <- sapply(1:in_sim_parms$npens, function(x){get_pen_mildclin(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results)})
  # Get amount of routine morbidity in each pen at given simulation time point
  pennormsick <- sapply(1:in_sim_parms$npens, function(x){get_pen_normsick(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normsick = sim_normsick)})
  # Get amount of routine mortality in each pen at given simulation time point
  pennormmort <- sapply(1:in_sim_parms$npens, function(x){get_pen_normmort(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normmort = sim_normmort,in_sim_parms = in_sim_parms )})
  
  
  if (is.na(in_selectedpensvec[1])){
    # Gets total number of disease dead, disease sick, routine sick and routine dead at the barn level
    premdismort <- sum(pendismort); premdissick <- sum(pendissick); premnormsick <- sum(pennormsick); premnormmort <- sum(pennormmort)
  } else {
    # Gets total number of disease dead, disease sick, routine sick and routine dead from user-inputted pens to target for sampling
    premdismort <- sum(pendismort[in_selectedpensvec]); premdissick <- sum(pendissick[in_selectedpensvec]); premnormsick <- sum(pennormsick[in_selectedpensvec]); premnormmort <- sum(pennormmort[in_selectedpensvec])
  }

  
  if (is.na(in_selectedpensvec[1])){
    # all pens in the barn used here
    if(in_sim_parms$samp_ind_viablood_not_oral){
      # blood samples taken
      # gets prevalence of viremic pigs
      premlindprev <- get_prem_lbloodprev(time_test_readtindex = time_test_readtindex, initer = initer, sample_type = sample_type, num_dissickobs=premdissick, transmission_results_premlevel =transmission_results_premlevel,pensizearray =  pensizearray)
    }else{
      # oral samples taken
      # gets prevalence of oral-positive pigs
      premlindprev <- get_prem_lindoralprev(time_test_readtindex = time_test_readtindex, initer = initer, sample_type = sample_type, num_dissickobs=premdissick, transmission_results_premlevel =transmission_results_premlevel,pensizearray =pensizearray  )  
    }
    
  } else {
    # only selected pens used here
    if(in_sim_parms$samp_ind_viablood_not_oral){
      # blood samples taken
      # gets prevalence of viremic pigs (algorithm changes with the sample type)
      premlindprev <- get_prem_lbloodprev_subsetpens(time_test_readtindex = time_test_readtindex, initer = initer, sample_type = sample_type, num_dissickobs=premdissick, transmission_results = transmission_results,in_pen_subset = in_selectedpensvec,pensizearray=pensizearray)
    }else{
      # oral samples taken
      # gets prevalence of oral-positive pigs
      premlindprev <- get_prem_lindoralprev_subsetpens(time_test_readtindex = time_test_readtindex, initer = initer, sample_type = sample_type, num_dissickobs=premdissick, transmission_results = transmission_results,in_pen_subset = in_selectedpensvec ,pensizearray=pensizearray)  
    }
    
  } 
  
  
  # initialize PCR test results
  temppremresults<-NA
  if(sample_type == 1){
    # Simulates PCR results with dead, sick, then apparently healthy sampling prioritization
    temppremresults <- perform_pcr_surv_dead_firstsick_livepooled(lprevaldbl = premlindprev, normdead = premnormmort, disdead = premdismort,
                                                                  normsick = premnormsick, dissick = premdissick, tempsampsize = insamplesize,boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 2){
    # Simulates PCR results with sick, dead, then apparently healthy pig sampling prioritization
    temppremresults <- perform_pcr_surv_sickfirstdead_livepooled(lprevaldbl = premlindprev, normdead = premnormmort, disdead = premdismort,
                                                                 normsick = premnormsick, dissick = premdissick, tempsampsize = insamplesize, boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 3){
    # Simulates PCR results with dead then random live pig sampling prioritization 
    temppremresults <- perform_pcr_surv_dead_firstsick_livepooled(lprevaldbl = premlindprev, normdead = premnormmort, disdead = premdismort,
                                                                  normsick = 0, dissick = 0, tempsampsize = insamplesize, boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 4){
    # Simulates PCR results with dead pig only sampling
    temppremresults <- perform_pcr_surv_dead_pooled(lprevaldbl = premlindprev, normdead = premnormmort, disdead = premdismort,
                                                tempsampsize = insamplesize, boolretnsamp=TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 5){
    # Simulates PCR results with dead then sick pig sampling prioritization
    temppremresults <- perform_pcr_surv_deadfirst_sick_pooled(lprevaldbl = premlindprev, normdead = premnormmort, disdead = premdismort,
                                                                  normsick = premnormsick, dissick = premdissick, tempsampsize = insamplesize,boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 6){
    # Simulates PCR results for sick then apparently healthy sampling prioritization
    temppremresults <- perform_pcr_surv_sickfirstdead_livepooled(lprevaldbl = premlindprev, normdead = 0, disdead = 0,
                                                                 normsick = premnormsick, dissick = premdissick, tempsampsize = insamplesize, boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 7){
    # Simulates PCR results for random live pig sampling
    temppremresults <- perform_pcr_surv_sickfirstdead_livepooled(lprevaldbl = premlindprev, normdead = 0, disdead = 0,
                                                                 normsick = 0, dissick = 0, tempsampsize = insamplesize, boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  }
  
  
  # returns desired results
  if(boolretnumsamp==FALSE){
    if(premresults > 0){
      return(1)
    } else{return(0)}
  }else{
    return(temppremresults)
  }
}


# This function simulates blood (or oral) and/or tissue sampling, sample pooling, and testing by PCR for a given simulation iteration
# and time post virus exposure when sampling occurs. Allows for mortality from previous days to be sampled and tested at 
# the same time
# Inputs: time_test_readtindex = simulation time index; initer = iteration number; sample_type = number indicating the
#   sampling prioritization scheme for dead/tissue sampling (details given below); insamplesize = total number of samples
#   to take; transmission_results = pen-level simulated transmission model output; sim_normmort = simulated routine mortality; 
#   sim_normsick = simulated routine morbidity; transmission_results_premlevel = premises-level simulated transmission
#   model output; instordays = number of days mortality was stored; inmaxdeadsamples = upper bound for dead pig samples;
#   boolretnumsamp = T/F for returning number of samples taken from sick, dead and live pigs; 
#   in_selectedpensvec = allows user to specify pens from which samples are taken; pensizearry = array with the number
#   of pigs in each pen; in_sim_parms = simulation parameters
# Outputs: If boolretnumsamp == TRUE, returns 1 if detected or 0 if not, the number of pooled samples tested by PCR,
#   and the number of samples taken from sick, dead, and live/apparently healthy pigs. If boolretnumsamp == FALSE,
#   only returns 1 or 0 if detected or not

# sample_type = 1: dead, sick, live sampling; sample_type = 2: sick, dead, live sampling; sample_type = 3: dead then live sampling;
# sample_type = 4: dead only; sample size =5: dead then sick sampling; sample type =6: sick then non sick live, sample_type = 7:
# random live sampling
sample_blood_iter_simtime_refri<-function(time_test_readtindex = 12, initer = 5, 
                                    sample_type=1,insamplesize=in_sim_parms$sampsize,transmission_results=transmission_results,
                                    sim_normmort=sim_normmort,
                                    sim_normsick=sim_normsick,
                                    transmission_results_premlevel=transmission_results_premlevel,instordays=5,inmaxdeadsamples=500,boolretnumsamp=FALSE,
                                    in_selectedpensvec=NA,pensizearray, in_sim_parms){

  # array with timepoints when daily mortality stored
  test_time_arr <- (time_test_readtindex-(0:(instordays-1))/in_sim_parms$readdt)
    

  # Checks if array of pen numbers inputted. If yes, only these pens are used for sampling
  if (is.na(in_selectedpensvec[1])){
    # array with ASF mortality for each day in test_time_arr
    dismort_vec <- sapply(test_time_arr, function(x){
      if(x>0){
        temp_pen_dismort <- sapply(1:in_sim_parms$npens, function(y){
          get_pen_dismort(time_test_readtindex = x, initer = initer, inpen = y, transmission_results = transmission_results,in_sim_parms = in_sim_parms)
        })
        return(sum(temp_pen_dismort))}else{
          return(0)
        }
    })
  
    # average amount of simulated normal mortality for the first 7 days of the simulation. Used to simulate routine
    # mortality for days prior to time zero (time of exposure)
    tfirstweekmeannormmortrate<-mean(apply(sim_normmort[[initer]][, 1:7],2,sum))
    # simulates routine mortality for each day where mortality is stored
    normmort_vec <- sapply(test_time_arr, function(x){
      if(x>0){
        temp_pen_normmort <- sapply(1:in_sim_parms$npens, function(y){
          get_pen_normmort(time_test_readtindex = x, initer = initer, inpen=y, sim_normmort = sim_normmort,in_sim_parms = in_sim_parms)
        })
        return(sum(temp_pen_normmort))}else{
          #taking mean of first one week as a poission rate for days prior to zero
          return(rpois(lambda = tfirstweekmeannormmortrate,n = 1))
        }
    })
    
  
  # gets pen-level ASF and routine morbidity  
  pendissick <- sapply(1:in_sim_parms$npens, function(x){get_pen_mildclin(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results)})
  pennormsick <- sapply(1:in_sim_parms$npens, function(x){get_pen_normsick(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normsick = sim_normsick )})
  # barn-level ASF and routine morbidity
  premdissick <- sum(pendissick); premnormsick <- sum(pennormsick); 
  
  
  # gets barn-level prevalence of viremic or oral-positive pigs. The algorithm changes with the sample type  
  if(in_sim_parms$samp_ind_viablood_not_oral){
    premlindprev <- get_prem_lbloodprev(time_test_readtindex = time_test_readtindex, initer = initer, sample_type = sample_type, num_dissickobs=premdissick, transmission_results_premlevel =transmission_results_premlevel,pensizearray=pensizearray )
  }else{
    premlindprev <- get_prem_lindoralprev(time_test_readtindex = time_test_readtindex, initer = initer, sample_type = sample_type, num_dissickobs=premdissick, transmission_results_premlevel =transmission_results_premlevel,pensizearray=pensizearray )  
  }
  
  } else {
    # only the inputted pens are used for sampling here
    
    # gets ASF mortality
    dismort_vec <- sapply(test_time_arr, function(x){
      if(x>0){
        
        temp_pen_dismort <- sapply(1:in_sim_parms$npens, function(y){
          get_pen_dismort(time_test_readtindex = x, initer = initer, inpen = y, transmission_results = transmission_results,in_sim_parms = in_sim_parms)
        })
        return(sum(temp_pen_dismort[in_selectedpensvec]))}else{
          return(0)
        }
    })
    
    # gets routine mortality
    tfirstweekmeannormmortrate<-mean(apply(sim_normmort[[initer]][, 1:7],2,sum))/length(pensizearray)
    #summing up for selected subsets only
    normmort_vec <- sapply(test_time_arr, function(x){
      if(x>0){
        temp_pen_normmort <- sapply(1:in_sim_parms$npens, function(y){
          get_pen_normmort(time_test_readtindex = x, initer = initer, inpen=y, sim_normmort = sim_normmort,in_sim_parms =in_sim_parms )
        })
        return(sum(temp_pen_normmort[in_selectedpensvec]))}else{
          #taking mean of first one week as a poission rate for days prior to zero
          return(rpois(lambda = tfirstweekmeannormmortrate*length(in_selectedpensvec),n = 1))
        }
    })
    
    
    # gets pen-level ASF and routine morbidity
    pendissick <- sapply(1:in_sim_parms$npens, function(x){get_pen_mildclin(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results)})
    pennormsick <- sapply(1:in_sim_parms$npens, function(x){get_pen_normsick(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normsick = sim_normsick )})
    # barn-level ASF-related and routine morbidity
    premdissick <- sum(pendissick[in_selectedpensvec]); premnormsick <- sum(pennormsick[in_selectedpensvec]); 
    
    
    # get barn-level blood- or oral-positive pig prevalence (algorithm changes depending on sample_type)
    if(in_sim_parms$samp_ind_viablood_not_oral){
      premlindprev <- get_prem_lbloodprev_subsetpens(time_test_readtindex = time_test_readtindex, initer = initer, sample_type = sample_type, num_dissickobs=premdissick, transmission_results =transmission_results,in_pen_subset = in_selectedpensvec,pensizearray = pensizearray )
    }else{
      premlindprev <- get_prem_lindoralprev_subsetpens(time_test_readtindex = time_test_readtindex, initer = initer, sample_type = sample_type, num_dissickobs=premdissick, transmission_results =transmission_results,in_pen_subset = in_selectedpensvec,pensizearray = pensizearray   )  
    } 
    
  }
  
  
  # initialize PCR results
  temppremresults<-NA
  if(sample_type == 1){
    # Simulates PCR results with dead, sick, then apparently healthy sampling prioritization
    temppremresults <- perform_pcr_surv_dead_firstsick_livepooled(lprevaldbl = premlindprev, normdead = sum(normmort_vec), disdead = sum(dismort_vec),
                                                                  normsick = premnormsick, dissick = premdissick, tempsampsize = insamplesize,boolretnsamp = TRUE,inmaxdeadsamples = inmaxdeadsamples,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 2){
    # Simulates PCR results with sick, dead, then apparently healthy pig sampling prioritization
    temppremresults <- perform_pcr_surv_sickfirstdead_livepooled(lprevaldbl = premlindprev, normdead = sum(normmort_vec), disdead = sum(dismort_vec),
                                                                 normsick = premnormsick, dissick = premdissick, tempsampsize = insamplesize,boolretnsamp = TRUE,inmaxdeadsamples = inmaxdeadsamples,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 3){
    # Simulates PCR results with dead then random live pig sampling prioritization 
    temppremresults <- perform_pcr_surv_dead_firstsick_livepooled(lprevaldbl = premlindprev, normdead = sum(normmort_vec), disdead = sum(dismort_vec),
                                                              normsick = 0, dissick = 0, tempsampsize = insamplesize,boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 4){
    # Simulates PCR results with dead pig only sampling
    temppremresults <- perform_pcr_surv_dead_pooled(lprevaldbl = premlindprev, normdead = sum(normmort_vec), disdead = sum(dismort_vec),
                                                tempsampsize = insamplesize, boolretnsamp=TRUE,in_sim_parms = in_sim_parms)

    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 5){
    # Simulates PCR results with dead then sick pig sampling prioritization
    temppremresults <- perform_pcr_surv_deadfirst_sick_pooled(lprevaldbl = premlindprev, normdead = sum(normmort_vec), disdead = sum(dismort_vec),
                                                              normsick = premnormsick, dissick = premdissick, tempsampsize = insamplesize,boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
  } else if(sample_type == 6){
    # Simulates PCR results for sick then apparently healthy sampling prioritization
    temppremresults <- perform_pcr_surv_sickfirstdead_livepooled(lprevaldbl = premlindprev, normdead = 0, disdead = 0,
                                                                 normsick = premnormsick, dissick = premdissick, tempsampsize = insamplesize, boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR
    premresults<-temppremresults$is_pcr_detected
    
  } else if(sample_type == 7){
    # Simulates PCR results for random live pig sampling
    temppremresults <- perform_pcr_surv_sickfirstdead_livepooled(lprevaldbl = premlindprev, normdead = 0, disdead = 0,
                                                                normsick = 0, dissick = 0, tempsampsize = insamplesize, boolretnsamp = TRUE,in_sim_parms = in_sim_parms)
    # 1 or 0 whether detection occurred by PCR  
    premresults<-temppremresults$is_pcr_detected
  }
  

  # return desired output
  if(boolretnumsamp==FALSE){
    if(premresults > 0){
      return(1)
    } else{return(0)}
  }else{
    return(temppremresults)
  }
}


# This function performs mortality trigger suveillance at a given simulation timepoint and for a given simulation
# iteration
# Input: time_test_readtindex = simulation time index; initer = simulation iteration number; transmission_results = 
#   pen-level transmission model output; sim_normmort = simulated routine mortality; transmission_results_premlevel =
#   barn-level transmission model output; inmorttrig = mortality trigger; tnumalive = number of pigs alive at timepoint;
#   in_sim_parms = simulation parameters;
# Output: 1 if mortality at simulation timepoint exceeds trigger threshold, 0 if not
perform_mort_surv_iter<-function(time_test_readtindex = 81, initer = 5, 
                                 transmission_results=transmission_results,
                                 sim_normmort=sim_normmort,
                                 transmission_results_premlevel=transmission_results_premlevel,
                                 inmorttrig,
                                 tnumalive,in_sim_parms=in_sim_parms){
 
  # get pen- and barn-level ASF-related and routine mortality
  pendismort <- sapply(1:in_sim_parms$npens, function(x){get_pen_dismort(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results,in_sim_parms =in_sim_parms)})  
  pennormmort <- sapply(1:in_sim_parms$npens, function(x){get_pen_normmort(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normmort = sim_normmort,in_sim_parms = in_sim_parms  )})
  premdismort <- sum(pendismort)
  premnormmort <- sum(pennormmort)
  
  # initialize whether detection by mortality trigger
  mortres<-0
  # checks if ASF-related and routine mortality exceeds the trigger amount
  if(premdismort+premnormmort>inmorttrig*tnumalive){
    mortres<-1
  }
  return(mortres)
  
}


# Inputs: time_test_readtindex = simulation time index; initer = simulation iteration number; inpen = pen number;
#   transmission_results = pen-level transmission model output for each iteration; in_sim_parms = simulation parameters
# Outputs: disease mortality accrued over the day prior to time_test_readtindex
get_pen_dismort <- function(time_test_readtindex = 110, initer = 5, inpen = 5,transmission_results=transmission_results, in_sim_parms){
  if(time_test_readtindex < 5){
    # if time_test_readtindex is within a day of virus exposure, mortality of 0 is returned
    return(transmission_results[[initer]][[inpen]]$dead[time_test_readtindex])
  } else{
    # Daily disease mortality accrued over a day
    return(transmission_results[[initer]][[inpen]]$dead[time_test_readtindex] - transmission_results[[initer]][[inpen]]$dead[(time_test_readtindex - 1/in_sim_parms$readdt)])
  }
}


# Input: time_test_readtindex = simulation time; initer = simulation iteration; inpen = pen number; sim_normmort = simulated
#  routine mortality; in_sim_parms = simulation parameters
# Output: normmortality occurring over the day prior to time_test_readtindex
get_pen_normmort <- function(time_test_readtindex = 110, initer = 5, inpen = 5,sim_normmort=sim_normmort, in_sim_parms){
  
  if(time_test_readtindex == 1){
    return(sim_normmort[[initer]][inpen, 1])
  } else{
    return(sim_normmort[[initer]][inpen, ceiling((time_test_readtindex-1)*in_sim_parms$readdt)])
  }
}


# Input: time_test_readtindex = simulation time; initer = iteration number; inpen = pen number; transmission_results = 
#   simulated pen-level transmission model output for each iteration
# Output: number of pigs with clinical signs related to ASF at the given timepoint and for the given iteration and pen
get_pen_mildclin <- function(time_test_readtindex = 110, initer = 5, inpen = 5,transmission_results=transmission_results){
  return(transmission_results[[initer]][[inpen]]$clinical[time_test_readtindex])
}

# Input: time_test_readtindex = simulation time; initer = iteration number; inpen = pen number; sim_normsick = 
#   simulated routine morbidity
# Output: routine morbidity at the given time point and simulation iteration and pen
get_pen_normsick <- function(time_test_readtindex = 110, initer = 5, inpen = 5,sim_normsick=sim_normsick){
  
  return(sim_normsick[[initer]][inpen,time_test_readtindex])
  
}


# Function for prevalence of blood positive among the apparently healthy (no clinical signs) or living pigs
# Inputs: time_test_readtindex = simulation time; initer = iteration number; sample_type = sampling prioritization scheme;
#   num_dissickobs = number of pigs with clinical signs at the simulation time; transmission_results_premlevel = barn-level
#   transmission model output; pensizearry = number of pigs per pen
# Outputs: prevalence of blood positive pigs among apparently healthy pigs or living pigs depending on the sampling
#   prioritization scheme
get_prem_lbloodprev <- function(time_test_readtindex = 110, initer = 5, sample_type=1, num_dissickobs = 0, transmission_results_premlevel=transmission_results_premlevel,pensizearray){
  # Checks if the sampling prioritization includes sampling sick pigs
  if(sample_type == 1 || sample_type == 2 || sample_type == 5||sample_type==6){
    
    # gets number of viremic pigs
    temp_blood <- transmission_results_premlevel[[initer]]$bloodpos[time_test_readtindex]
    # gets number morbidity
    temp_mildclin <- num_dissickobs
    # gets number of living pigs (alive = starting pen size - cumulative dead)
    temp_alive <- sum(pensizearray) - transmission_results_premlevel[[initer]]$dead[time_test_readtindex]
    
    # returns prevalence of viremic pigs among the apparently healthy (not showing clinical signs) population
    if ((temp_blood - temp_mildclin) >= 0 || (temp_alive - temp_mildclin) > 0){
      return((temp_blood - temp_mildclin) / (temp_alive - temp_mildclin))
    } else{
      return(0)
    }
  } else if (sample_type == 3 || sample_type == 4 || sample_type == 7){ 
    #if sampling only done on live/dead or dead only
    
    # gets number of viremic pigs
    temp_blood <- transmission_results_premlevel[[initer]]$bloodpos[time_test_readtindex]
    # gets number of living pigs (alive = starting pen size - cumulative dead)
    temp_alive <- sum(pensizearray) - transmission_results_premlevel[[initer]]$dead[time_test_readtindex]
    
    # returns prevalence of viremic pigs among all living pigs
    if(temp_alive > 0){
      return(temp_blood / temp_alive)
    } else{
      return(0)
    }
  }
  
}


# Function for prevalence of blood positive among the apparently healthy (no clinical signs) or living pigs in selected pens
# Inputs: time_test_readtindex = simulation time; initer = iteration number; sample_type = sampling prioritization scheme;
#   num_dissickobs = number of pigs with clinical signs at the simulation time; transmission_results = pen-level
#   transmission model output; in_pen_subset = selected pens to target for sampling; pensizearry = number of pigs per pen
# Outputs: prevalence of blood positive pigs among apparently healthy pigs or living pigs depending on the sampling
#   prioritization scheme in the selected pens
get_prem_lbloodprev_subsetpens <- function(time_test_readtindex = 110, initer = 5, sample_type=1, num_dissickobs=0, transmission_results=transmission_results,in_pen_subset=NA,pensizearray){
  # If sampling includes sick pigs
  if(sample_type == 1 || sample_type == 2 || sample_type == 5||sample_type==6){

    # gets number of viremic pigs
    temp_blood <- sum(do.call("c",lapply(in_pen_subset,function(h){transmission_results[[initer]][[h]]$bloodpos[time_test_readtindex]})))
    # gets number morbidity
    temp_mildclin <- num_dissickobs
    # gets number of living pigs (alive = starting pen size - cumulative dead)
    temp_alive <- sum(pensizearray[in_pen_subset]) - sum(do.call("c",lapply(in_pen_subset,function(h){transmission_results[[initer]][[h]]$dead[time_test_readtindex]})))
    
    # returns prevalence of viremic pigs among the apparently healthy (not showing clinical signs) population
    if ((temp_blood - temp_mildclin) >= 0 || (temp_alive - temp_mildclin) > 0){
      return((temp_blood - temp_mildclin) / (temp_alive - temp_mildclin))
    } else{
      return(0)
    }
  } else if (sample_type == 3 || sample_type == 4 || sample_type == 7){
     #if sampling only done on live/dead or dead only
    
    # gets number of viremic pigs
    temp_blood <- sum(do.call("c",lapply(in_pen_subset,function(h){transmission_results[[initer]][[h]]$bloodpos[time_test_readtindex]})))
    # # gets number of living pigs (alive = starting pen size - cumulative dead)
    temp_alive <- sum(pensizearray[in_pen_subset]) - sum(do.call("c",lapply(in_pen_subset,function(h){transmission_results[[initer]][[h]]$dead[time_test_readtindex]})))
    
    # returns prevalence of viremic pigs among all living pigs
    if(temp_alive > 0){
      return(temp_blood / temp_alive)
    } else{
      return(0)
    }
  }
  
}



# Function for prevalence of infectious alive pigs at the premises-level. Note that infectiousness can be based on nasal and oral swabs in data parametrization if interested in oral sampling
# Inputs: time_test_readtindex = simulation time; initer = iteration number; sample_type = sampling prioritization scheme;
#   num_dissickobs = number of pigs with clinical signs at the simulation time; transmission_results_premlevel = barn-level
#   transmission model output; pensizearry = number of pigs per pen
# Outputs: prevalence of infectious pigs among apparently healthy pigs or living pigs depending on the sampling
#   prioritization scheme
get_prem_lindoralprev <- function(time_test_readtindex = 110, initer = 5, sample_type=1, num_dissickobs = 0, transmission_results_premlevel=transmission_results_premlevel,pensizearray){
  # If sampling includes sick pigs
  if(sample_type == 1 || sample_type == 2 || sample_type == 5||sample_type == 6){
    # get number of infectious pigs
    temp_blood <- transmission_results_premlevel[[initer]]$infrec[time_test_readtindex]+transmission_results_premlevel[[initer]]$infdead[time_test_readtindex]
    # get total number of pigs with clinical signs
    temp_mildclin <- num_dissickobs
    # get number of living pigs (alive = starting pen size - dead)
    temp_alive <- sum(pensizearray) - transmission_results_premlevel[[initer]]$dead[time_test_readtindex]
    
    if ((temp_blood - temp_mildclin) >= 0 || (temp_alive - temp_mildclin) > 0){
      # prevalence of infectiousness in apparently healthy (no clinical signs) pigs
      return((temp_blood - temp_mildclin) / (temp_alive - temp_mildclin))
    } else{
      return(0)
    }
  } else if (sample_type == 3 || sample_type == 4 || sample_type == 7){ 
    #if sampling only done on live/dead or dead only i.e., for purely random live regardless of sick or non sick
    
    # get blood-positive prevalence
    temp_blood <- transmission_results_premlevel[[initer]]$infrec[time_test_readtindex]+transmission_results_premlevel[[initer]]$infdead[time_test_readtindex]
    # get number of alive pigs (alive = starting pen size - dead)
    temp_alive <- sum(pensizearray) - transmission_results_premlevel[[initer]]$dead[time_test_readtindex]
    
    # return infection prevalence
    if(temp_alive > 0){
      return(temp_blood / temp_alive)
    } else{
      return(0)
    }
  }
  
}


# Function to get the infection prevalence based on a subset of pens rather than the premises. Either prevalence in 
# apparently healthy or among all live pigs depending on the protocol selected
# Inputs: time_test_readtindex = simulation time; initer = iteration number; sample_type = sampling prioritization scheme;
#   num_dissickobs = number of pigs with clinical signs at the simulation time; transmission_results = pen-level
#   transmission model output; in_pen_subset = selected pens to target for sampling; pensizearry = number of pigs per pen
# Outputs: prevalence of infectious pigs among apparently healthy pigs or living pigs depending on the sampling
#   prioritization scheme in the selected pens
get_prem_lindoralprev_subsetpens <- function(time_test_readtindex = 110, initer = 5, sample_type=1, num_dissickobs=0, transmission_results=transmission_results,in_pen_subset=NA,pensizearray){
  # If sampling includes sick pigs
  if(sample_type == 1 || sample_type == 2 || sample_type == 5||sample_type==6){
    # get number of infectious pigs in the selected pens
    temp_blood <- sum(do.call("c",lapply(in_pen_subset,function(h){transmission_results[[initer]][[h]]$infrec[time_test_readtindex]+transmission_results[[initer]][[h]]$infdead[time_test_readtindex]})))
    # get number of pigs with clinical signs in the selected pens
    temp_mildclin <- num_dissickobs
    # get number of living pigs (alive = starting pen size - dead)
    temp_alive <- sum(pensizearray[in_pen_subset]) - sum(do.call("c",lapply(in_pen_subset,function(h){transmission_results[[initer]][[h]]$dead[time_test_readtindex]})))
    
    if ((temp_blood - temp_mildclin) >= 0 || (temp_alive - temp_mildclin) > 0){
      # prevalence of infectious pigs in apparently healthy (no clinical signs) pigs
      return((temp_blood - temp_mildclin) / (temp_alive - temp_mildclin))
    } else{
      return(0)
    }
  } else if (sample_type == 3 || sample_type == 4 || sample_type == 7){ 
    #if sampling only done on live/dead or dead only
    
    # get number of infectious pigs in the selected pens
    temp_blood <- sum(do.call("c",lapply(in_pen_subset,function(h){transmission_results[[initer]][[h]]$infrec[time_test_readtindex]+transmission_results[[initer]][[h]]$infdead[time_test_readtindex]})))
    # get number of living pigs (alive = starting pen size - dead)
    temp_alive <- sum(pensizearray[in_pen_subset]) - sum(do.call("c",lapply(in_pen_subset,function(h){transmission_results[[initer]][[h]]$dead[time_test_readtindex]})))
    
    # return infection prevalence
    if(temp_alive > 0){
      return(temp_blood / temp_alive)
    } else{
      return(0)
    }
  }
  
}



# This function simulates PCR tests with dead, sick, then apparently healthy sampling prioritization. Samples from dead
# pigs pooled separately from sick and apparently healthy samples.
# Inputs: lprevaldbl = prevalence of viremic or infectious pigs; normdead = number of dead pigs due routine causes;
#   disdead = number of dead pigs due to ASF; normsick = number of pigs with clinical signs due to routine causes;
#   dissick = number of pigs with clinical signs due to ASF; tempsampsize = total number of samples; boolretnsamp =
#   T/F for whether to return number of sick/dead/apparently healthy pigs sampled; inmaxdeadsamples = ceiling for
#   the number of dead pigs that can be sampled; in_sim_parms = simulation parameters
# Outputs: list with the variables is_pcr_detected = 1 or 0 if detection occurred; N_pospcR_tubes = the total number
#   of pooled samples that tested positive; If boolretnsamp == TRUE, also returns outnumsamplesvec = the number of
#   samples taken from sick, dead, and live pigs
perform_pcr_surv_dead_firstsick_livepooled<-function(lprevaldbl=0.0,normdead=0,disdead=0,normsick=0,dissick=0,
                                                     tempsampsize=30,boolretnsamp=FALSE,inmaxdeadsamples=500, in_sim_parms){
  
  # number of pigs dead due to routine causes
  Snormal<-normdead
  # number of pigs dead due to ASF
  sdisease<-disdead
  
  
  ## your normal mortality may happen to have ASF, adjustment for that
  if(is.na(lprevaldbl)){
    disease_among_normal = 0
  } else{
    disease_among_normal = rbinom(1,Snormal, lprevaldbl)
  }
  sdisease = sdisease + disease_among_normal
  Snormal = Snormal - disease_among_normal
  totaldead_pigs<-sdisease+Snormal
  
  # total number of pigs with clinical signs (due to ASF or routine causes)
  totalsickpigs<-normsick+dissick
  
  
  # prioritize dead pigs for testing first
  deadpigstested<-min(totaldead_pigs,tempsampsize)
  # put cap on dead pigs to test
  deadpigstested<-min(inmaxdeadsamples,deadpigstested)
  
  # next, prioritize sick pigs for testing
  sigpigstest<-min(totalsickpigs,tempsampsize-deadpigstested)
  # lastly, prioritize apparently healthy pigs
  livepigstest<-tempsampsize- deadpigstested-sigpigstest
  
  # store number to be tested from each sub-population
  outnumsamplesvec<-c(deadpigstested, sigpigstest,livepigstest)

  
  # pool dead pig tissue samples in pools of 5
  N_tubes_deadpcr_on_day<-ceiling(deadpigstested/in_sim_parms$pcr_targetswabs)
  # initialize detection variable
  detected=FALSE
  # check if there is any dead pig testing to be done
  if (N_tubes_deadpcr_on_day>0){
    # initialize number of dead pigs in each pooled sample
    ith_n_arr_dead<-rep(0,  N_tubes_deadpcr_on_day)
    
    # divide dead pig samples among total number of pools
    ith_n_arr_dead[]<-floor(( deadpigstested) /  N_tubes_deadpcr_on_day)
    zeta = deadpigstested - (floor((deadpigstested) /  N_tubes_deadpcr_on_day) *  N_tubes_deadpcr_on_day)
    if(zeta>0){
      ith_n_arr_dead[1:zeta]<- ith_n_arr_dead[1:zeta]+1
    }
    
    # initialize array indicating if detection occurred for each pooled sample
    PCR_res_arr_dead<-rep(0,N_tubes_deadpcr_on_day)
    
    # loop through number of pooled samples
    for (i in 1 : N_tubes_deadpcr_on_day){
      # initialize variable for if positive
      pcrpositive=0
      
      # simulate the number of samples from dead pigs due to ASF in pooled sample i
      snumberdintube = wraphyper(Snormal, sdisease, ith_n_arr_dead[i])
      # remaining pool size consists of samples from dead pigs due to routine causes
      n_norm_in_tube = ith_n_arr_dead[i] - snumberdintube
      if (snumberdintube > 0 ){
        # if there is at least one sample in the pool from a pig that died due to ASF, simulate the result of PCR
        # testing based on the test sensitivity
        pcrpositive = rbinom(1, 1,in_sim_parms$PCRsens)
      }

      # update the total number of dead pigs due to routine and ASF-related causes available for sampling
      Snormal = Snormal - (n_norm_in_tube)
      if (Snormal < 0) {Snormal = 0}
      sdisease = sdisease - (snumberdintube)
      
      # if detected, update flag for overall detection and array indicating whether each pooled sample tested positive
      if( pcrpositive > 0){detected = TRUE
      PCR_res_arr_dead[i]<-1
      }
    }
    
  } else {PCR_res_arr_dead <-0}
  
  ############################################end dead##################################
  
  # Simulate PCR testing of blood samples (pooled separately from tissue samples from dead pigs)
  # Blood from sick and apparently healthy pigs pooled together
  
  # total number of pooled blood samples
  N_sicklivetubes_onaday<-ceiling(( sigpigstest+livepigstest)/in_sim_parms$pcr_targetswabs)
  if(N_sicklivetubes_onaday>0 ){
    # number of sick pigs in each pooled sample
    ith_n_sickarr<-rep(0, N_sicklivetubes_onaday)
    # number of apparently healthy pigs in each pooled sample
    ith_n_livearr<-rep(0, N_sicklivetubes_onaday)
    
    # total number of pigs with clinical signs due to ASF and non-ASF related causes
    snormsick<-normsick
    sdissick<-dissick
    
    
    # divide sick pig samples among the pools
    ith_n_sickarr[]<-floor((sigpigstest) /  N_sicklivetubes_onaday)
    szeta = sigpigstest - (floor((sigpigstest) /  N_sicklivetubes_onaday) *  N_sicklivetubes_onaday)
    if(szeta>0){
      ith_n_sickarr[(1):szeta]<-ith_n_sickarr[(1):szeta ]+1
    }
    
    # divide apparently healthy pig samples among the pools
    ith_n_livearr[]<-floor((livepigstest) /  N_sicklivetubes_onaday)
    lzeta = livepigstest - (floor((livepigstest) /  N_sicklivetubes_onaday) *  N_sicklivetubes_onaday)
    if(lzeta>0){
      ith_n_livearr[(N_sicklivetubes_onaday-lzeta+1):N_sicklivetubes_onaday]<-ith_n_livearr[(N_sicklivetubes_onaday-lzeta+1):N_sicklivetubes_onaday]+1
    }
    
    
    ##actual testing for sick and live
    # initialize PCR result for each pool
    PCR_res_sicklivearr<-rep(0, N_sicklivetubes_onaday)
    
    # loop through pooled samples
    for (i in 1 : N_sicklivetubes_onaday){
      # initialize PCR result for a pool
      pcrpositive=0
      
      # simulate the number of pigs with clinical signs due to ASF in pooled sample i
      snumberdintube_sick = wraphyper(snormsick, sdissick, ith_n_sickarr[i])
      # remaining samples are due to pigs with clinical signs due to routine causes
      n_norm_in_tube_sick = ith_n_sickarr[i] - snumberdintube_sick
      
      
      # simulate the number of apparently healthy pigs in pooled sample i that are viremic
      templbirdstobetested =ith_n_livearr[i] 
      if (templbirdstobetested>0){
        if(is.na(lprevaldbl)){
          nposlivebirds_pertube = 0
        } else {
          nposlivebirds_pertube = rbinom(1,templbirdstobetested, lprevaldbl)
        }
      }else { nposlivebirds_pertube=0}
      
      
      # if there is at least one ASF-positive pig in the sample, simulate the PCR result using the test sensitivity
      if (    snumberdintube_sick>0 || nposlivebirds_pertube>0){
        pcrpositive = rbinom(1, 1,in_sim_parms$PCRsens)
      }
      
      
      # update the number of pigs with clinical signs due to ASF and non-ASF related causes available for sampling
      snormsick = snormsick - (n_norm_in_tube_sick)
      if (snormsick< 0) {snormsick = 0}
      sdissick = sdissick - (snumberdintube_sick)
      
      
      # if detected, update variable for overall detection and array with indicator for detection for each pooled sample
      if( pcrpositive > 0){detected = TRUE
      PCR_res_sicklivearr[i]<-1
      }
      
    }
  } else {PCR_res_sicklivearr<-0}
  
  
  if (detected){resdetected<-1}else {resdetected<-0}
  
  # return desired output
  if(boolretnsamp==TRUE){
  return(list(is_pcr_detected=resdetected, N_pospcR_tubes=  sum(PCR_res_sicklivearr)+sum(PCR_res_arr_dead),outnumsamplesvec=outnumsamplesvec))
  }else{
  return(list(is_pcr_detected=resdetected, N_pospcR_tubes=  sum(PCR_res_sicklivearr)+sum(PCR_res_arr_dead)))
  }
}


# This function simulates PCR tests with sick, dead, then apparently healthy sampling prioritization. Samples from dead
# pigs pooled separately from sick and apparently healthy samples.
# Inputs: lprevaldbl = prevalence of viremic or infectious pigs; normdead = number of dead pigs due routine causes;
#   disdead = number of dead pigs due to ASF; normsick = number of pigs with clinical signs due to routine causes;
#   dissick = number of pigs with clinical signs due to ASF; tempsampsize = total number of samples; boolretnsamp =
#   T/F for whether to return number of sick/dead/apparently healthy pigs sampled; inmaxdeadsamples = ceiling for
#   the number of dead pigs that can be sampled; in_sim_parms = simulation parameters
# Outputs: list with the variables is_pcr_detected = 1 or 0 if detection occurred; N_pospcR_tubes = the total number
#   of pooled samples that tested positive; If boolretnsamp == TRUE, also returns outnumsamplesvec = the number of
#   samples taken from sick, dead, and live pigs
perform_pcr_surv_sickfirstdead_livepooled<-function(lprevaldbl=0.0,normdead=0,disdead=0,normsick=0,dissick=0,
                                                    tempsampsize=30,boolretnsamp=FALSE,inmaxdeadsamples=500, in_sim_parms){
  
  # number of dead pigs due to ASF
  Snormal<-normdead
  # number of dead pigs due to routine causes
  sdisease<-disdead
  
  
  ## your normal mortality may happen to have LPAI, adjustment for that
  if(is.na(lprevaldbl)){
    disease_among_normal = 0
  } else {
    disease_among_normal = rbinom(1,Snormal, lprevaldbl)
  }
  sdisease = sdisease + disease_among_normal
  Snormal = Snormal - disease_among_normal
  totaldead_pigs<-sdisease+Snormal
  
  # total number of pigs with clinical signs due to ASF-related and non-ASF-related causes
  totalsickpigs<-normsick+dissick
  
  
  # sick pigs prioritized first for sampling
  sigpigstest<-min(totalsickpigs,tempsampsize)
  # next, dead pigs selected for sampling
  deadpigstested<-min(totaldead_pigs,tempsampsize-sigpigstest)
  deadpigstested<-min(inmaxdeadsamples,deadpigstested)
  # lastly, number of apparently healthy pigs sampled determined
  livepigstest<-tempsampsize- deadpigstested-sigpigstest
  # stores number of samples taken from dead, sick, and apparently healthy pigs
  outnumsamplesvec<-c(deadpigstested, sigpigstest,livepigstest)

  ##### Dead pig testing ######
  # determine total number of pools for dead pig testing
  N_tubes_deadpcr_on_day<-ceiling(deadpigstested/in_sim_parms$pcr_targetswabs)
  detected=FALSE #key flag has to be false before any testing
  # check that there are pooled samples to be tested
  if (N_tubes_deadpcr_on_day>0){
    # initialize number of dead pigs in each pooled sample
    ith_n_arr_dead<-rep(0,  N_tubes_deadpcr_on_day)
    
    # divides dead pigs among the total number of samples
    ith_n_arr_dead[]<-floor(( deadpigstested) /  N_tubes_deadpcr_on_day)
    zeta = deadpigstested - (floor((deadpigstested) /  N_tubes_deadpcr_on_day) *  N_tubes_deadpcr_on_day)
    if(zeta>0){
      ith_n_arr_dead[1:zeta]<- ith_n_arr_dead[1:zeta]+1
    }
    
    # initialize PCR result for each pooled sample
    PCR_res_arr_dead<-rep(0,N_tubes_deadpcr_on_day)
    
    # loop through PCR pooled samples
    for (i in 1 : N_tubes_deadpcr_on_day){
      # initialize PCR result
      pcrpositive=0
      
      # simulate the number of dead pigs due to ASF pooled sample i
      snumberdintube = wraphyper(Snormal, sdisease, ith_n_arr_dead[i])
      # remaining dead pig samples due to non-ASF reasons
      n_norm_in_tube = ith_n_arr_dead[i] - snumberdintube
      # if there is at least one sample from a pig dead with ASF, PCR result simulated based on the test sensitivity
      if (snumberdintube > 0 ){
        pcrpositive = rbinom(1, 1,in_sim_parms$PCRsens)
      }

      # update mortality due to ASF and non-ASF reasons available for sampling
      Snormal = Snormal - (n_norm_in_tube)
      if (Snormal < 0) {Snormal = 0}
      sdisease = sdisease - (snumberdintube)
      
      # update variable for overall detection and array with PCR result for each pooled sample
      if( pcrpositive > 0){detected = TRUE
      PCR_res_arr_dead[i]<-1
      }
      
    }
    
  } else {PCR_res_arr_dead <-0}
  
  ############################################end dead##################################
  
  
  ##### Sick and apparently healthy testing #####
  # sick and apparently healthy samples pooled separately from dead pig samples
  
  # determine number of pooled samples
  N_sicklivetubes_onaday<-ceiling(( sigpigstest+livepigstest)/in_sim_parms$pcr_targetswabs)
  if(N_sicklivetubes_onaday>0 ){
    # initialize number of sick pigs in each pooled sample
    ith_n_sickarr<-rep(0, N_sicklivetubes_onaday)
    # initialize number of apparently healthy pigs in each pooled sample
    ith_n_livearr<-rep(0, N_sicklivetubes_onaday)
    
    # set the number of pigs with clinical signs due to non-ASF reasons
    snormsick<-normsick
    # set the number of pigs with clinical signs due to ASF
    sdissick<-dissick
    
    # fill the pooled samples with sick pig samples
    ith_n_sickarr[]<-floor((sigpigstest) /  N_sicklivetubes_onaday)
    szeta = sigpigstest - (floor((sigpigstest) /  N_sicklivetubes_onaday) *  N_sicklivetubes_onaday)
    if(szeta>0){
      ith_n_sickarr[(1):szeta]<-ith_n_sickarr[(1):szeta ]+1
    }
    
    # fill the pooled samples with apparently healthy samples
    ith_n_livearr[]<-floor((livepigstest) /  N_sicklivetubes_onaday)
    lzeta = livepigstest - (floor((livepigstest) /  N_sicklivetubes_onaday) *  N_sicklivetubes_onaday)
    if(lzeta>0){
      ith_n_livearr[(N_sicklivetubes_onaday-lzeta+1):N_sicklivetubes_onaday]<-ith_n_livearr[(N_sicklivetubes_onaday-lzeta+1):N_sicklivetubes_onaday]+1
    }
    
    
    ##actual testing for sick and live
    # initialize array with PCR results for each pooled sample
    PCR_res_sicklivearr<-rep(0, N_sicklivetubes_onaday)
    # loop through pooled samples
    for (i in 1 : N_sicklivetubes_onaday){
      # initialize PCR result
      pcrpositive=0
      
      # simulate number of pigs with clinical signs due to ASF in pooled sample i
      snumberdintube_sick = wraphyper(snormsick, sdissick, ith_n_sickarr[i])
      # remainder of samples in pooled sample are from pigs with clinical signs due to non-ASF reasons
      n_norm_in_tube_sick = ith_n_sickarr[i] - snumberdintube_sick
      
      
      # simulate number of apparently healthy pigs in sample that are ASF-positive
      templbirdstobetested =ith_n_livearr[i] 
      if (templbirdstobetested>0){
        if(is.na(lprevaldbl)){
          nposlivebirds_pertube = 0
        } else {
          nposlivebirds_pertube = rbinom(1,templbirdstobetested, lprevaldbl)
        }
        
      }else { nposlivebirds_pertube=0}
      
      # if at least one sample is included from a pig with ASF, simulates the PCR test result based on the test sensitivity
      if (snumberdintube_sick>0 || nposlivebirds_pertube>0){
        pcrpositive = rbinom(1, 1,in_sim_parms$PCRsens)
      }

     
      # updates the number of pigs with clinical signs with and without ASF available for sampling
      snormsick = snormsick - (n_norm_in_tube_sick)
      if (snormsick< 0) {snormsick = 0}
      sdissick = sdissick - (snumberdintube_sick)
      
      # updates PCR results
      if( pcrpositive > 0){detected = TRUE
      PCR_res_sicklivearr[i]<-1
      }
 
    }
  } else {PCR_res_sicklivearr<-0}
  
  # returns desired output
  if (detected){resdetected<-1}else {resdetected<-0}
  if(boolretnsamp==FALSE){
    return(list(is_pcr_detected=resdetected, N_pospcR_tubes=  sum(PCR_res_sicklivearr)+sum(PCR_res_arr_dead)))
  }else{
    return(list(is_pcr_detected=resdetected, N_pospcR_tubes=  sum(PCR_res_sicklivearr)+sum(PCR_res_arr_dead),outnumsamplesvec=outnumsamplesvec)) 
  }
}



# This function simulates PCR testing with samples from dead pigs only
# Inputs: lprevaldbl = prevalence of viremic or infectious pigs; normdead = number of dead pigs due routine causes;
#   disdead = number of dead pigs due to ASF; tempsampsize = total number of samples; boolretnsamp =
#   T/F for whether to return number of sick/dead/apparently healthy pigs sampled; in_sim_parms = simulation parameters
# Outputs: list containing the variables is_pcr_detected = 1 or 0 if detection occurred; N_pospcR_tubes = the total number
#   of pooled samples that tested positive; If boolretnsamp == TRUE, also returns outnumsamplesvec = the number of
#   samples taken from sick, dead, and live pigs
perform_pcr_surv_dead_pooled<-function(lprevaldbl=0.0, normdead=0,disdead=0,tempsampsize=30, boolretnsamp=FALSE, 
                                       in_sim_parms){
  
  # dead pigs due to routine causes
  Snormal<-normdead
  # dead pigs due to ASF
  sdisease<-disdead
  
  
  ## your normal mortality may happen to have ASF, adjustment for that
  if(is.na(lprevaldbl)){
    disease_among_normal = 0
  } else{
    disease_among_normal = rbinom(1,Snormal, lprevaldbl)
  }
  sdisease = sdisease + disease_among_normal
  Snormal = Snormal - disease_among_normal
  # total number of dead pigs
  totaldead_pigs<-sdisease+Snormal
  
  
  # total number of samples from dead pigs
  deadpigstested<-min(totaldead_pigs,tempsampsize)
  
  # stores number of samples from dead, sick, and apparently healthy pigs
  outnumsamplesvec<-c(deadpigstested, 0,0)
  
  
  # number of pooled samples
  N_tubes_deadpcr_on_day<-ceiling(deadpigstested/in_sim_parms$pcr_targetswabs)
  # initializes boolean variable that tracks whether detection occurred
  detected=FALSE 
  # checks if there is any testing to be done
  if (N_tubes_deadpcr_on_day>0){
    # initializes array with the number of dead pig samples per pool
    ith_n_arr_dead<-rep(0,  N_tubes_deadpcr_on_day)
    
    # fills the array with the number of dead pig samples per pool
    ith_n_arr_dead[]<-floor(( deadpigstested) /  N_tubes_deadpcr_on_day)
    zeta = deadpigstested - (floor((deadpigstested) /  N_tubes_deadpcr_on_day) *  N_tubes_deadpcr_on_day)
    if(zeta>0){
      ith_n_arr_dead[1:zeta]<- ith_n_arr_dead[1:zeta]+1
    }
    
    # initializes array with whether or not ASF was detected in each pooled sample
    PCR_res_arr_dead<-rep(0,N_tubes_deadpcr_on_day)
    
    # loopes through the number of pooled samples
    for (i in 1 : N_tubes_deadpcr_on_day){
      # initializes PCR result
      pcrpositive=0
      
      # simulates the number of samples from pigs with ASF included in the pooled sample
      snumberdintube = wraphyper(Snormal, sdisease, ith_n_arr_dead[i])
      # the rest of the samples come from pigs not infected with ASF
      n_norm_in_tube = ith_n_arr_dead[i] - snumberdintube
      
      # if there is at least one sample from a detectable pig with ASF, a PCR result is simulated based on the diagnostic senstivity
      if (snumberdintube > 0 ){
        pcrpositive = rbinom(1, 1,in_sim_parms$PCRsens)
      }

      
      # updates the number of ASF mortality and non-ASF mortality available for sampling
      Snormal = Snormal - (n_norm_in_tube)
      if (Snormal < 0) {Snormal = 0}
      sdisease = sdisease - (snumberdintube)
      
      
      # if detection occurs, updates variables
      if( pcrpositive > 0){detected = TRUE
      PCR_res_arr_dead[i]<-1
      }
      
    }
  } else {PCR_res_arr_dead <-0}
  
  
  # returns desired results
  if (detected){resdetected<-1}else {resdetected<-0}
  if(boolretnsamp==FALSE){
    return(list(is_pcr_detected=resdetected, N_pospcR_tubes=sum(PCR_res_arr_dead)))
  }else{
    return(list(is_pcr_detected=resdetected, N_pospcR_tubes=sum(PCR_res_arr_dead),outnumsamplesvec=outnumsamplesvec)) 
  }
}



# This function simulates PCR testing with dead then sick pig sampling prioritization
# Inputs: lprevaldbl = prevalence of viremic or infectious pigs; normdead = number of dead pigs due routine causes;
#   disdead = number of dead pigs due to ASF; normsick = number of pigs with clinical signs due to routine causes;
#   dissick = number of pigs with clinical signs due to ASF; tempsampsize = total number of samples; boolretnsamp =
#   T/F for whether to return number of sick/dead/apparently healthy pigs sampled; inmaxdeadsamples = ceiling for
#   the number of dead pigs that can be sampled; in_sim_parms = simulation parameters
# Outputs: list with the variables is_pcr_detected = 1 or 0 if detection occurred; N_pospcR_tubes = the total number
#   of pooled samples that tested positive; If boolretnsamp == TRUE, also returns outnumsamplesvec = the number of
#   samples taken from sick, dead, and live pigs
perform_pcr_surv_deadfirst_sick_pooled<-function(lprevaldbl=0.0,normdead=0,disdead=0,normsick=0,dissick=0,tempsampsize=30,
                                                 boolretnsamp=FALSE,inmaxdeadsamples=500, in_sim_parms){
  
  # routine mortality
  Snormal<-normdead
  # mortality related to ASF
  sdisease<-disdead
  
  
  ## your normal mortality may happen to have ASF, adjustment for that
  if(is.na(lprevaldbl)){
    disease_among_normal = 0
  } else{
    disease_among_normal = rbinom(1,Snormal, lprevaldbl)
  }
  sdisease = sdisease + disease_among_normal
  Snormal = Snormal - disease_among_normal
  
  # total number of available dead pigs
  totaldead_pigs<-sdisease+Snormal
  # total number of observed pigs with clinical signs
  totalsickpigs<-normsick+dissick
  
  
  # dead pigs are sampled first up to inmaxdeadsamples
  deadpigstested<-min(totaldead_pigs,tempsampsize)
  deadpigstested<-min(inmaxdeadsamples,deadpigstested)
  
  # sick pigs are priortized next for sampling
  sigpigstest<-min(totalsickpigs,tempsampsize-deadpigstested)

  # stores the number of samples from dead, sick, and apparently healthy pigs
  outnumsamplesvec<-c(deadpigstested, sigpigstest,0)

  
  # number of pooled samples from dead pigs (pooled separately from blood samples from sick pigs)
  N_tubes_deadpcr_on_day<-ceiling(deadpigstested/in_sim_parms$pcr_targetswabs)
  # initialize variable indicating whether detection occurred
  detected=FALSE 
  # checks if there are any dead pig samples to be tested
  if (N_tubes_deadpcr_on_day>0){
    
    # initialize array with the number of dead pig tissue samples in each pooled sample
    ith_n_arr_dead<-rep(0,  N_tubes_deadpcr_on_day)
    # determines number of dead pig tissue samples in each pooled sample
    ith_n_arr_dead[]<-floor(( deadpigstested) /  N_tubes_deadpcr_on_day)
    zeta = deadpigstested - (floor((deadpigstested) /  N_tubes_deadpcr_on_day) *  N_tubes_deadpcr_on_day)
    if(zeta>0){
      ith_n_arr_dead[1:zeta]<- ith_n_arr_dead[1:zeta]+1
    }
    
    # initialize array with indicator for detection for each pooled sample
    PCR_res_arr_dead<-rep(0,N_tubes_deadpcr_on_day)
    
    
    # loop through the pooled samples
    for (i in 1 : N_tubes_deadpcr_on_day){
      # initialize variable indicating detection
      pcrpositive=0
      
      # simulate the number of samples from detectable pigs with ASF
      snumberdintube = wraphyper(Snormal, sdisease, ith_n_arr_dead[i])
      # remaining samples in pool from pigs not detectable for ASF (uninfected or infected but not viremic)
      n_norm_in_tube = ith_n_arr_dead[i] - snumberdintube
      
      # if at least one sample from a detectable pig was included in the sample, the test result is simulated based on the
      # PCR diagnostic test sensitivity
      if (snumberdintube > 0 ){
        pcrpositive = rbinom(1, 1,in_sim_parms$PCRsens)
      }

      
      # updates the number of mortality due to ASF and not due to ASF available for sampling
      Snormal = Snormal - (n_norm_in_tube)
      if (Snormal < 0) {Snormal = 0}
      sdisease = sdisease - (snumberdintube)
      
      
      # if the sample was PCR positive, variables are updated indicating this
      if( pcrpositive > 0){detected = TRUE
      PCR_res_arr_dead[i]<-1
      }
      
    }
    
  } else {PCR_res_arr_dead <-0}
  

  
  # number of pooled samples for sick pig testing
  N_sicklivetubes_onaday<-ceiling(( sigpigstest)/in_sim_parms$pcr_targetswabs)
  if(N_sicklivetubes_onaday>0 ){
    
    # initializes number of blood samples from sick pigs in each pooled sample
    ith_n_sickarr<-rep(0, N_sicklivetubes_onaday)

    # number of sick pigs due to routine causes
    snormsick<-normsick
    # number of sick pigs due to ASF
    sdissick<-dissick
    
    # divides blood samples between the pooled samples
    ith_n_sickarr[]<-floor((sigpigstest) /  N_sicklivetubes_onaday)
    szeta = sigpigstest - (floor((sigpigstest) /  N_sicklivetubes_onaday) *  N_sicklivetubes_onaday)
    if(szeta>0){
      ith_n_sickarr[(1):szeta]<-ith_n_sickarr[(1):szeta ]+1
    }
    

    # initializes array indicating detection for each pooled sample
    PCR_res_sicklivearr<-rep(0, N_sicklivetubes_onaday)
    
    
    # loops through pooled samples
    for (i in 1 : N_sicklivetubes_onaday){
      # initializes indicator for detection
      pcrpositive=0
      
      # simulates the number of samples from sick pigs with ASF included in the pooled sample
      snumberdintube_sick = wraphyper(snormsick, sdissick, ith_n_sickarr[i])
      # remaining samples from pigs not detectable for ASF
      n_norm_in_tube_sick = ith_n_sickarr[i] - snumberdintube_sick
      
      
      # if there is at least one sample from a sick pig with ASF, the PCR result is simulated based on the diagnostic
      # test sensitivity
      if (snumberdintube_sick>0){
        pcrpositive = rbinom(1, 1,in_sim_parms$PCRsens)
      }
      

      # updates the number of pigs with clinical signs detectable and not detectable for ASF available for sampling
      snormsick = snormsick - (n_norm_in_tube_sick)
      if (snormsick< 0) {snormsick = 0}
      sdissick = sdissick - (snumberdintube_sick)
      
      
      # if PCR was detected, updates variables indicating this
      if( pcrpositive > 0){detected = TRUE
      PCR_res_sicklivearr[i]<-1
      }
      
    }
  } else {PCR_res_sicklivearr<-0}
  
  
  # stores if detection occurred by dead or sick pig sampling
  if (detected){resdetected<-1}else {resdetected<-0}
  
  
  # returns desired results
  if(boolretnsamp==TRUE){
    return(list(is_pcr_detected=resdetected, N_pospcR_tubes=  sum(PCR_res_sicklivearr)+sum(PCR_res_arr_dead),outnumsamplesvec=outnumsamplesvec))
  }else{
    return(list(is_pcr_detected=resdetected, N_pospcR_tubes=  sum(PCR_res_sicklivearr)+sum(PCR_res_arr_dead)))
  }
}
