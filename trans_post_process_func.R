# This file contains helper functions for processing transmission model results in file "Transmission model visualization.R"
# Functions mainly have to do with when pens have at least one infectious pig. The remaining function identifies the 
# iteration where infection did not take off (defined as having mortality < 10)


# Function for getting the time until the first infectious pig in each pen ordered from shortest to longest time
# Inputs: initer = simulation iteration number; transmission_results = pen-level transmission model output for each
#   iteration
# Outputs: peninforder = pens from shortest to longest time until first infectious pig; sorted_peninftimes = time of
#   first infectious pig in the pen from shortest to longest
get_infectious_ordpens<-function(initer=1,transmission_results=transmission_results){
  # array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # Inputs: inpenno = pen number
  # Outputs: first time in tseq when the number of infectious pigs was greater than 0
  get_time_of_first_inbfectious<-function(inpenno=1){
    return(tseq[ (transmission_results[[initer]][[inpenno]][[4]]+transmission_results[[initer]][[inpenno]][[5]]>0)][1])
  }
  
  # array with the time of the first infectious pig in each pen
  pen_inftimes<-sapply(1:sim_parms$npens,get_time_of_first_inbfectious)
  # if a pen never had an infectious pig, set the time to the number of days run in the simulation
  pen_inftimes[is.na(pen_inftimes)]<-sim_parms$ndays
  # array with the pen numbers from shortest to longest time until the first infectious pig
  peninforder=order(pen_inftimes,1:sim_parms$npens)
  return(list( peninforder= peninforder,sorted_peninftimes=pen_inftimes[peninforder]))
}


# Inputs: initer = simulation iteration number; transmission_results = pen-level transmission model output for each
#   iteration
# Outputs: Number of pens with the number of infectious pigs exceeding 0 for each simulation time point
get_num_infectious_pens_iter<-function(initer=1,transmission_results=transmission_results){
  # array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # Inputs: pen number
  # Outputs: T/F if the total number of infectious pigs in the pen was larger than 0 at each time point
  get_if_infectious_pen<-function(inpenno=1){
   (transmission_results[[initer]][[inpenno]][[4]]+transmission_results[[initer]][[inpenno]][[5]]>0)
  }
  
  # get T/F for when infectious pigs exceeded 0 for each pen
  pen_infstatus_mat<-sapply(1:sim_parms$npens,get_if_infectious_pen)
  # for each time point, sum number of pens with infectious pigs at that time point
  Num_inf_penvec<-apply(pen_infstatus_mat,1,sum,na.rm=TRUE)
 
  return(Num_inf_penvec)
}


# Inputs: initer = iteration number; intimept = simulation time point; transmission_results = pen-level
#   transmission model output
# Outputs: vector with T/F for which iterations had greater than 0 infectious pigs for a given time point
get_which_infectious_pens_iter<-function(initer=1,intimept=10,transmission_results=transmission_results){
  # array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # Inputs: pen number
  # Outputs: T/F if the total number of infectious pigs in the pen was larger than 0 at each time point
  get_if_infectious_pen<-function(inpenno=1){
    (transmission_results[[initer]][[inpenno]][[4]]+transmission_results[[initer]][[inpenno]][[5]]>0)
  }
  
  # get T/F for when infectious pigs exceeded 0 for each pen
  pen_infstatus_mat<-sapply(1:sim_parms$npens,get_if_infectious_pen)
  # vector with T/F for which iterations had greater than 0 infectious pigs for a given time point
  outvector<-as.vector(pen_infstatus_mat[intimept,])
  # label the pens
  names(outvector)<-paste(rep("pen",sim_parms$npens),1:sim_parms$npens)
  return(outvector)
}


# Inputs: in_prem_transresults = barn-level transmission output
# Outputs: proportion of transmission runs that had cumulative mortality greater than 10
has_epdiemic_diedout<-function(in_prem_transresults=global_trans_list$transmission_results_premlevel){
  
  mean(sapply(1:sim_parms$num_iterations,function(h){if(in_prem_transresults[[h]][[1]][1]-in_prem_transresults[[h]][[1]][length(in_prem_transresults[[h]][[1]])]>=10){
    diedout<-FALSE
  }else{
    diedout<-TRUE
  }
  return(diedout)}))
  
}

