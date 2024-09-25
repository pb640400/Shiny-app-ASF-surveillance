# The functions in this file make tables giving the mean and 95% P.I. of the number of pigs with clinical signs (total
# and due to ASF alone) and number of dead pigs (total and due to ASF alone)


# Clinical signs data table
# Inputs: sim_parms = simulation parameters; transmission_results_premlevel = barn-level transmission model output;
#   sim_normsick = simulated routine morbidity; outdays = the days post virus exposure to store results
# Outputs: table with the mean and 90% P.I. of the total number of pigs with clinical signs and the number of pigs with
#   clinical signs related to ASF only
getdatatable_clin_prem <- function(sim_parms,transmission_results_premlevel=transmission_results_premlevel, sim_normsick=sim_normsick, outdays=NA ){
  # array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  
  # matrix with the total routine morbidity each day for each iteration (row iteration, column day)
  normsickmat <- t(sapply(1:length(sim_normsick), function(x){
    colSums(sim_normsick[[x]])
  }))
  # matrix with the number of pigs with clinical signs over time per iteration
  tmatrix_dissick <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]]$clinical})) 
  # matrix with the total number of pigs with clinical signs (due to ASF and not) over time post virus exposure for each iteration
  tmatrix_total <- t(sapply(1:sim_parms$num_iterations, function(h){
    normsickmat[h,] + tmatrix_dissick[h,]
  }))
  
  
  # determines the mean and 90% P.I. for the total number of pigs with clinical signs and the number of pigs with
  # ASF-related clinical signs
  df<-data.frame(
    mean_state_total=apply(tmatrix_total,2,mean,na.rm=TRUE),
    state_total_95_percent=apply(tmatrix_total,2,quantile,na.rm=TRUE,probs=0.95),
    state_total_5_percent=apply(tmatrix_total,2,quantile,na.rm=TRUE,probs=0.05),
    mean_state_dissick=apply(tmatrix_dissick,2,mean,na.rm=TRUE),
    state_dissick_95_percent=apply(tmatrix_dissick,2,quantile,na.rm=TRUE,probs=0.95),
    state_dissick_5_percent=apply(tmatrix_dissick,2,quantile,na.rm=TRUE,probs=0.05),
    'DPE'=tseq
  )
  cdf<-df
  # selects the days post virus exposure specified by outdays
  cdf<-cdf[ tseq[outdays/sim_parms$readdt+2],]
  # formats output
  cdf$'Pigs_with_clinical_signs'=paste(round(cdf$mean_state_total,2),"(",round(cdf$state_total_5_percent),"-",(cdf$state_total_95_percent),")",sep  ="")
  cdf$'Detectable_pigs_with_clinical_signs'=paste(round(cdf$mean_state_dissick,2),"(",round(cdf$state_dissick_5_percent),"-",round(cdf$state_dissick_95_percent),")",sep  ="")
  
  
  # selects columns of data frame to output
  toutdf=as.data.frame((cdf)[,c(7,8,9)])
  # names selection
  names(toutdf)[2]<-('Pigs with clinical signs,mean (90% prediction interval)')
  names(toutdf)[3]<-('ASF pigs with clinical signs,mean (90% prediction interval)')
  
  
  return((toutdf))
}


# Daily mortality data table
# Inputs: sim_parms = simulation parameters; transmission_results_premlevel = number of birds over time post exposure in
#   the different disease states at the barn level; sim_normdead = simulated routine mortality; outdays = days post exposure
#   selected for the output
# Outputs: table with the mean and 90% P.I. of the total number of dead pigs and the number of dead pigs related to ASF only
getdatatable_dead_prem <- function(sim_parms,transmission_results_premlevel=transmission_results_premlevel, sim_normdead=sim_normdead, outdays=NA ){
  
  # array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  
  # matrix with the total number of dead pigs across the pens over time for each iteration (row iteration, column day)
  normdeadmat <- t(sapply(1:length(sim_normdead), function(x){
    colSums(sim_normdead[[x]])
  }))
  
  
  # mortality is output as cumulative in the transmission model. This code snippet adjusts the output to get
  # daily mortality
  tmatrix_dismort <- do.call('rbind',lapply(1:sim_parms$num_iterations, function(h){
    zeropadmort <- c(0, transmission_results_premlevel[[h]]$dead[-(sim_parms$ndays+1)])
    # Get daily  mort instead of total mort
    temp_dailydismort <- transmission_results_premlevel[[h]]$dead - zeropadmort
    return(temp_dailydismort)
  }))
  
  
  # gets the cumulative mortality (routine + ASF mortality) over time for each iteration
  tmatrix_totalmort <- t(sapply(1:sim_parms$num_iterations, function(h){
   # browser()
    normdeadmat[h,] + tmatrix_dismort[h,]
  }))
  
  
  # calculates the mean and 90% P.I. over time for the ASF and cumulative mortality
  df<-data.frame(
    mean_state_totalmort=apply(tmatrix_totalmort,2,mean,na.rm=TRUE),
    state_totalmort_95_percent=apply(tmatrix_totalmort,2,quantile,na.rm=TRUE,probs=0.95),
    state_totalmort_5_percent=apply(tmatrix_totalmort,2,quantile,na.rm=TRUE,probs=0.05),
    mean_state_dismort=apply(tmatrix_dismort,2,mean,na.rm=TRUE),
    state_dismort_95_percent=apply(tmatrix_dismort,2,quantile,na.rm=TRUE,probs=0.95),
    state_dismort_5_percent=apply(tmatrix_dismort,2,quantile,na.rm=TRUE,probs=0.05),
    'DPE'=tseq
  )
  
  
  cdf<-df
  # selects the days to print
  cdf<-cdf[ tseq[outdays/sim_parms$readdt+2],]
  # formats the data
  cdf$'Dead_pigs'=paste(round(cdf$mean_state_totalmort,2),"(",round(cdf$state_totalmort_5_percent),"-",round(cdf$state_totalmort_95_percent),")",sep  ="")
  cdf$'Detectable_dead_pigs'=paste(round(cdf$mean_state_dismort,2),"(",round(cdf$state_dismort_5_percent),"-",round(cdf$state_dismort_95_percent),")",sep  ="")
  
  # selects the columns to output
  toutdf=as.data.frame((cdf)[,c(7,8,9)])
  names(toutdf)[2]<-('Daily dead pigs,mean (90% prediction interval)')
  names(toutdf)[3]<-('Daily dead pigs due to ASF,mean (90% prediction interval)')
  
  
  return((toutdf))
}



