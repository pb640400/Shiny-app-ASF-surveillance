# This file contains functions for simulating oral fluids PCR testing



# A function that returns an oral fluid PCR test sensitivity for a given pen-level prevalence
# Input: inprev = a pen-level prevalence; sens_frame = data from a csv file with pen level prevalences and associated
#   oral fluid sensitivities
# Output: a oral fluid sensitivity for the given pen-level prevalence. If prevalence not in sens_frame, estimated by
#   linear interpolation
get_sen_penprev<-function (inprev=.3,sens_frame){
  # if prevalence below first data point in sens_frame, sensitivity set to zero
  if (inprev<sens_frame[1,1]){
    return(0)
  }else if (inprev>=sens_frame[nrow(sens_frame),1]){
    # if prevalence is above last data point in sens_frame, sensitivity set to sensitivity of last data point in sens_frame
    return(sens_frame[nrow(sens_frame),2])
  }else{
    # if prevalence between first and last data point in sens_frame, sensitivity determined using linear interpolation between
    # the sensitivities of the prevalences in sens_frame that bracket inprev
    start_index<-max(which(sens_frame[,1]<=inprev))
    tout<-sens_frame[start_index,2]+ (inprev-sens_frame[start_index,1])/(sens_frame[start_index+1,1]-sens_frame[start_index,1])*((sens_frame[start_index+1,2]-sens_frame[start_index,2]))
    return(tout)
  }
  
}


# Input: time_test_readtindex = simulation time index; inpen = number id of pen in barn; transmission_results = the number
#   of pigs in each clinical state over time per pen for each simulation iteration; initer = simulation iteration number,
#   pensizearray = the number of pigs in each pen
# Output: the infection prevalence for a given simulation time, simulation iteration, and pen
get_pen_prev<-function(time_test_readtindex,inpen,initer,transmission_results=transmission_results,pensizearray=pensizearray){
  talive<-pensizearray[inpen]-transmission_results[[initer]][[inpen]][[7]][time_test_readtindex]
  if(talive>0){
    return((transmission_results[[initer]][[inpen]][[4]][time_test_readtindex]+transmission_results[[initer]][[inpen]][[5]][time_test_readtindex])/talive)
  } else {
    return(NA)
  }
}


# Simulates oral fluid PCR restult for a given pen-level prevalence
# Input: inprev = a pen-level prevalence; sens_frame = data containing pen-level prevalence and associated PCR sensitivities
#   for oral fluid testing
# Output: returns a 1 if oral fluid tests positive and 0 if oral fluid tests negative
perform_one_penropesamp<-function(inprev,sens_frame=sens_frame){
  
  if (is.na(inprev)){
    return(NA)
  } else {
    # determine test sensitivity based on pen-level prevalence
    tsens<-get_sen_penprev(inprev =inprev,sens_frame = sens_frame )
    # simulate PCR test based on tsens
    rbinom(n = 1,size = 1,prob =tsens )
  }
}



# Simulates oral fluid testing for a given simulation iteration and simulation time
# Inputs: time_test_readtindex = simulation time index; initer = simulation iteration number; npens_sampled = number of
#   pens to sample for oral fluids; sim_normmort = simulated routine mortality; sim_normsick = simulated routine morbidity;
#   transmission_results = the number of pigs in each clinical state over time per pen for each simulation iteration;
#   in_OFsampletype = number for oral fluids sampling scenario (1 is for random pen sampling, 2 prioritizes pens with
#   sick pigs; 3 prioritizes pens with dead pigs, 4 prioritizes pens with dead pigs first then pens with sick pigs);
#   in_selectedpensvec = array with pens to be tested by oral fluids; pensizearray = number of pigs per pen;
#   in_sim_parms = simulation parameters; sens_frame = dataframe with pen level prevalences and oral fluid PCR test senstivities
# Outputs: list with variables ropedet = 1 if detection by oral fluid occurred and 0 otherwise; sampled_pens = the pen
#   number of the sampled pens
sample_oral_iter_simtime<-function(time_test_readtindex = 110,initer = 5,npens_sampled, sim_normmort=sim_normmort, sim_normsick = sim_normsick, 
                                   transmission_results=transmission_results,in_OFsampletype=1,in_selectedpensvec=NA,pensizearray = pensizearray,in_sim_parms,sens_frame){
  
  # array with the infection prevalence for each pen at the given simulation time
  penprevalences<-sapply(1:in_sim_parms$npens,function(x){get_pen_prev(time_test_readtindex =time_test_readtindex,initer = initer,inpen = x,transmission_results = transmission_results,pensizearray = pensizearray)})
  
  
  # selects pens to be tested
 if (is.na(in_selectedpensvec[1])){
  sampled_pens <- get_sampledpens(time_test_readtindex = time_test_readtindex, initer = initer, npens_sampled=npens_sampled,
                               sim_normmort=sim_normmort, sim_normsick = sim_normsick, transmission_results = transmission_results,
                               OFsampletype = in_OFsampletype,in_sim_parms =in_sim_parms  )
 } else{
   sampled_pens<-in_selectedpensvec
 }
  
  # simulates PCR result for each pen to be tested
  penresults<-sapply(sampled_pens,function(x){perform_one_penropesamp(inprev = penprevalences[x],sens_frame =sens_frame )})
  
  
  if (sum(penresults,na.rm = TRUE)>0){
    return(list("ropedet"=1,"sampled_pens"=sampled_pens))
  }else{return(list("ropedet"=0, "sampled_pens"=sampled_pens))}
  
}


# This function selects pens to be sampled depending on the sampling prioritization protocol
# OFsampletype definitions: 1 == random pen samples; 2 == sick pig pen prioritization; 3 == dead pig pen prioritization; 4 == dead then sick pig pen priortization
# Input: time_test_readtindex = simulation time index; initer = simulation iteration number; npens_sampled = number of
#   pens to be tested; sim_normmort = simulated routine mortality; sim_normsick = simulated routine morbidity;
#   transmission_results = the number of pigs in each clinical state over time per pen for each simulation iteration;
#   OFsampletype = number indicating the sampling prioritization; in_sim_parms = simulation parameters
# Output: array with pen indices to be sampled
get_sampledpens <- function(time_test_readtindex = 110, initer = 5, npens_sampled=in_sim_parms$npens, sim_normmort=sim_normmort, sim_normsick=sim_normsick, 
                            transmission_results=transmission_results, OFsampletype = 1,in_sim_parms){
  
  if(OFsampletype == 1){
    #random pen sampling
    sampled_pens<-sample(x = 1:in_sim_parms$npens,replace = FALSE,size = npens_sampled)
  } else if(OFsampletype == 2){
    #sick pig pen prioritization
    
    #from Blood_detect_exp.R file
    # selects ASF morbidity and routine morbidity 
    pendissick <- sapply(1:in_sim_parms$npens, function(x){get_pen_mildclin(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results)})
    pennormsick <- sapply(1:in_sim_parms$npens, function(x){get_pen_normsick(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normsick = sim_normsick)})
    
    totalsick_pen <- pendissick + pennormsick
    penswithsick <- which(totalsick_pen > 0)
    
    #if more pens have sick pigs than the number of pens to sample, then selects the pens with the highest numbers of sick
    if(length(penswithsick) >= npens_sampled){
      sampled_pens <- order(totalsick_pen, decreasing = TRUE)[1:npens_sampled]
    }else{
      #if there are fewer pens with sick pigs than the number of pens we need to sample, the pens with sick pigs are selected
      # and the remaining pens are selected randomly. Need special condition to do all random if there are no pens with sick pigs.
      
      if(length(penswithsick) == 0){
        sample_pensnosick <- sample(x=c(1:in_sim_parms$npens),replace=FALSE,size = (npens_sampled - length(penswithsick)))
        sampled_pens <- c(sample_pensnosick)
        
      } else{
        pensnosick <- c(1:in_sim_parms$npens)[-penswithsick]
        sample_pensnosick <- sample(x=pensnosick,replace=FALSE,size = (npens_sampled - length(penswithsick)))
        
        sampled_pens <- c(penswithsick, sample_pensnosick)
      }

    }
  } else if(OFsampletype == 3){#dead pig pen prioritization
    
    #from Blood_detect_exp.R file
    # selects ASF mortality and routine mortality
    pendismort <- sapply(1:in_sim_parms$npens, function(x){get_pen_dismort(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results,in_sim_parms =in_sim_parms )})
    pennormmort <- sapply(1:in_sim_parms$npens, function(x){get_pen_normmort(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normmort = sim_normmort,in_sim_parms =in_sim_parms )})
    
    totalmort_pen <- pendismort + pennormmort
    penswithdead <- which(totalmort_pen > 0)
    
    #if more pens have dead pigs than the number of pens to sample, then select the pens with the highest numbers of dead
    if(length(penswithdead) >= npens_sampled){
      sampled_pens <- order(totalmort_pen, decreasing = TRUE)[1:npens_sampled]
    }else{
      #if more pens sampled then have dead pigs, then sample all pens with dead then random
      
      #if no pens have dead, do all random
      if(length(penswithdead) == 0){
        sample_pensnodead <- sample(x=c(1:in_sim_parms$npens),replace = FALSE,size = npens_sampled)
        sampled_pens <- c(sample_pensnodead)
      } else{
        pensnodead <- c(1:in_sim_parms$npens)[-penswithdead]
        sample_pensnodead <- sample(x=pensnodead,replace=FALSE,size = (npens_sampled - length(penswithdead)))
        
        sampled_pens <- c(penswithdead, sample_pensnodead)
      }
      
    }
  } else if(OFsampletype == 4){ #dead then sick pen prioritization
    
    #from Blood_detect_exp.R file
    # Selects ASF-related and routine mortality
    pendismort <- sapply(1:in_sim_parms$npens, function(x){get_pen_dismort(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results,in_sim_parms =in_sim_parms )})
    pennormmort <- sapply(1:in_sim_parms$npens, function(x){get_pen_normmort(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normmort = sim_normmort,in_sim_parms=in_sim_parms )})
    
    totalmort_pen <- pendismort + pennormmort
    penswithdead <- which(totalmort_pen > 0)
    
    #if more pens with dead pigs than pens to sample, select pens with highest numbers of dead
    if(length(penswithdead) >= npens_sampled){
      sampled_pens <- order(totalmort_pen, decreasing = TRUE)[1:npens_sampled]
    } else{
      
      # Get sick per pen
      pendissick <- sapply(1:in_sim_parms$npens, function(x){get_pen_mildclin(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,transmission_results = transmission_results)})
      pennormsick <- sapply(1:in_sim_parms$npens, function(x){get_pen_normsick(time_test_readtindex = time_test_readtindex, initer = initer, inpen = x,sim_normsick = sim_normsick)})
      
      totalsick_pen <- pendissick + pennormsick
      # identify pens with sick but no dead
      penswithsick <- which(totalsick_pen > 0 & totalmort_pen == 0)
      
      
      ndeadsamples <- length(penswithdead)
      #if dead + sick more than pens to sample, sample all dead then top off remaining with sick (highest first)
      if(length(penswithsick) >= (npens_sampled - ndeadsamples)){
        
        temp_totalsick_pen <- totalsick_pen
        temp_totalsick_pen[-penswithsick] <- 0
        
        sicksamples <- order(temp_totalsick_pen, decreasing = TRUE)[1:(npens_sampled - ndeadsamples)]
        
        sampled_pens <- c(penswithdead, sicksamples)
      } else{
        #if dead + sick is not more than pens to sample, then do all dead, all sick, and top off with random
        
        nsicksamples <- length(penswithsick)
        
        # need condition if we have no pens with sick or dead
        if((nsicksamples + ndeadsamples) == 0){
          healthysamples <- sample(x=c(1:in_sim_parms$npens),replace=FALSE,size = npens_sampled)
          sampled_pens <- c(healthysamples)
        } else{
          pensnodead_nosick <- c(1:in_sim_parms$npens)[-c(penswithdead, penswithsick)] #integer(0) if no penswithdead or penswithsick
          healthysamples <- sample(x=pensnodead_nosick,replace=FALSE,size = (npens_sampled - ndeadsamples - nsicksamples))
          
          sampled_pens <- c(penswithdead, penswithsick, healthysamples)
        }
        
      }
    }
    
  } #end dead>sick prioritization if statement

  return(sampled_pens)
}

  

