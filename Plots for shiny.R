# This file contains code for generating plots and tables for the number of pigs in different disease states over time
# post exposure both for a single transmission run and for comparing two transmission runs


# Plots comparison of pigs in a disease state based on two sets of transmission model results
# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = transmission
#   model output at the barn level; past_transmission_results_premlevel = another set of transmission model out at the barn level
# Outputs: prints a plot displaying the mean number and 90% interval of pigs in a given disease state over time for
#   two sets of transmission results (results given on the same plot in a different color)
Comb_Visualize_geom_ribbon_state_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel,
                                                past_transmission_results_premlevel=transmission_results_premlevel ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # makes a matrix with the total number of pigs in the disease state at each simulation time point for each iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[instate]]}))
  tmatrix2 <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){past_transmission_results_premlevel[[h]][[instate]]}))
  
  # makes data frames with the mean and 90% P.I. of the number of pigs in the disease state at each time point
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  df2<-data.frame(
    mean_state=apply(tmatrix2,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  df$Scenario='current'
  df2$Scenario='past'
  # combines data frames from the two transmission runs
  cdf<-rbind(df,df2)
  names(cdf)[1]<-paste('Mean ',pstate_list_names[instate],sep="")
  names(cdf)[2]<-paste('perc_95th',pstate_list_names[instate],sep="")
  names(cdf)[3]<-paste('perc_5th',pstate_list_names[instate],sep="")
  cdf$dfill='blue'
  cdf$dfill=factor(cdf$dfill)
  
  # plots the mean number of pigs over time in the disease state with the 90% interval for each set of results 
  pl<-ggplot(data=cdf,aes_string(x=as.name(names(cdf)[4]),y=as.name(names(cdf)[1]),fill=names(cdf)[5]))
  pl<-pl+geom_line(size=1.2, aes(linetype=Scenario,color=Scenario))
  pl<-pl+labs(x="Time post exposure (Days)")
  pl<-pl+scale_color_manual( values=c('slateblue','orangered4'))
  pl<-pl+labs(color='Scenario')
  pl<-pl+theme(axis.title.x=element_text(margin = margin(t = 15)))
  pl<-pl+theme(axis.title.y=element_text(margin = margin(r = 15)))
  pl<-pl+theme(strip.text.y = element_text(size = 15))
  pl<-pl+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  pl<-pl+theme(panel.background = element_rect(fill="ivory"))
  pl<-pl+theme(axis.text=element_text(size=12),
               axis.title=element_text(size=14),
               legend.text=element_text(size=12),
               legend.title=element_text(size=14),
               legend.position = "top")
  dfill<-'dodgerblue2'
  pl<-pl+geom_ribbon(data = cdf, aes_string(x=as.name(names(cdf)[4]), ymin=as.name(names(cdf)[3]),ymax=as.name(names(cdf)[2]),fill=as.name(names(cdf)[5])),alpha=0.25)+
    scale_fill_manual(name='Scenario', values=c('slateblue','orangered4'),guide='legend')
  print(pl)
  
}


# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = transmission
#   model output at the barn level; past_transmission_results_premlevel = another set of transmission model out at the barn level
# Outputs: prints a plot of the number of infectious pigs over time post virus exposure for each transmission model run
  Comb_Visualize_geom_ribbon_infectious_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel,
                                                     past_transmission_results_premlevel=transmission_results_premlevel ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # matrices with the number of infectious pigs (infectious to recover and infectious to die) over time for each iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[4]]+
    transmission_results_premlevel[[h]][[5]]}))
  tmatrix2 <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){past_transmission_results_premlevel[[h]][[4]]+
    past_transmission_results_premlevel[[h]][[5]]}))
  
  # derives mean and 90% interval for the infectious pigs at each time point
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  df2<-data.frame(
    mean_state=apply(tmatrix2,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  df$Scenario='current'
  df2$Scenario='past'
  
  instate<-4
  # combines data frames from each transmission model run
  cdf<-rbind(df,df2)
  names(cdf)[1]<-'Infectious pigs'
  names(cdf)[2]<-paste('perc_95th',state_list_names[instate],sep="")
  names(cdf)[3]<-paste('perc_5th',state_list_names[instate],sep="")
  
  
  # plots mean and 90% interval of the number of infectious pigs over time post exposure for each transmission run
  pl<-ggplot(data=cdf,aes_string(x=as.name(names(cdf)[4]),y=as.name(names(cdf)[1]),fill=names(cdf)[5]))
  pl<-pl+geom_line(size=1.2, aes(linetype=Scenario,color=Scenario))
  pl<-pl+labs(x="Time post exposure (Days)")
  pl<-pl+scale_color_manual( values=c('slateblue','orangered4'))
  pl<-pl+labs(color='Scenario')
  pl<-pl+theme(axis.title.x=element_text(margin = margin(t = 15)))
  pl<-pl+theme(axis.title.y=element_text(margin = margin(r = 15)))
  pl<-pl+theme(strip.text.y = element_text(size = 15))
  pl<-pl+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  pl<-pl+theme(panel.background = element_rect(fill="ivory"))
  pl<-pl+theme(axis.text=element_text(size=12),
               axis.title=element_text(size=14),
               legend.text=element_text(size=12),
               legend.title=element_text(size=14),
               legend.position = "top")
  dfill<-'dodgerblue2'
  pl<-pl+geom_ribbon(data = cdf, aes_string(x=names(cdf)[4], ymin=names(cdf)[3],ymax=names(cdf)[2],fill=names(cdf)[5]),alpha=0.25)+
    scale_fill_manual(name='Scenario', values=c('slateblue','orangered4'),guide='legend')
  print(pl)
  
}


# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = transmission
#   model output at the barn level; outdays = selects days to output results
# Outputs: table with the mean and 90% P.I. number of infectious pigs for selected days post virus exposure 
getdatatable_infectious_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel, outdays=NA ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # matrix with the number of infectious pigs in the barn over time post virus exposure for each simulation iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[4]]+
    transmission_results_premlevel[[h]][[5]]}))
  
 
 # determines the mean and 90% P.I. for the number of infectious pigs over time 
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    'DPE'=tseq
  )
  
 # selects output from days given in outdays
  instate<-4
  cdf<-df
  cdf<-cdf[ tseq[outdays/sim_parms$readdt+2],]
  # formats output
  names(cdf)[1]<-'Infectious_pigs'
  names(cdf)[2]<-paste('perc_95th',state_list_names[instate],sep="")
  names(cdf)[3]<-paste('perc_5th',state_list_names[instate],sep="")
  cdf$'Infectious pigs'=paste(round(cdf$Infectious_pigs,1),"(",round(cdf$perc_5thinfrec),"-",round(cdf$perc_95thinfrec),")",sep  ="")

  # selects columns to output  
  toutdf=as.data.frame((cdf)[,c(4,5)])
  names(toutdf)[2]<-('Infectious pigs,mean (95% prediction interval)')

  return((toutdf))
}


# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = transmission
#   model output at the barn level; past_transmission_results_premlevel = another set of transmission model out at the barn 
#   level; outdays = selects days to output results
# Outputs: table with the mean and 90% P.I. for the number of infectious pigs over time for each of the transmission model
#   runs
Comb_getdatatable_infectious_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel,
                                            past_transmission_results_premlevel=transmission_results_premlevel,outdays=NA ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # matrices with the number infectious pigs over time post virus exposure for each simulation iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[4]]+
    transmission_results_premlevel[[h]][[5]]}))
  tmatrix2 <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){past_transmission_results_premlevel[[h]][[4]]+
    past_transmission_results_premlevel[[h]][[5]]}))
  
  # data frames with mean and 90% P.I. of the number of infectious pigs over time post virus exposure
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  df2<-data.frame(
    mean_state=apply(tmatrix2,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  # formats output
  df$'Infectious pigs current'=paste(round(df$mean_state,1),"(",round(df$state_5_percentile),"-",round(df$state_95_percentile),")",sep  ="")
  df2$'Infectious pigs previous'=paste(round(df2$mean_state,1),"(",round(df2$state_5_percentile),"-",round(df2$state_95_percentile),")",sep  ="")
  
  # combines data frames from the two transmission model runs
  cdf<-data.frame(DPE=df$times,`Infectious pigs current` = df$'Infectious pigs current',
                  `Infectious pigs previous`=df2$'Infectious pigs previous',row.names = NULL)
  names(cdf)[3]<-('Infectious pigs previous simulation mean (95% prediction interval)')
  names(cdf)[2]<-('Infectious pigs current simulation mean (95% prediction interval)')
  # selects days to output data
  cdf<-cdf[ tseq[outdays/sim_parms$readdt+2],]
  
  return(cdf)
}


# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = transmission
#   model output at the barn level; outdays = selects days to output results
# Outputs: table with the mean and 90% P.I. of the number of pigs disease state 'instate' at the times post virus exposure
# given by 'outdays'
get_datatbl_state_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel,
                                 outdays=NA ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # matrix with the number of pigs in the disease state indicated by 'instate' over time post virus exposure for each
  # iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[instate]]}))
  
  # data frame with the mean and 90% P.I. for the number of pigs over time post virus exposure in the disease state
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    DPE=tseq
  )
  
  # Selects data to output
  cdf<-df
  cdf<-cdf[ tseq[outdays/sim_parms$readdt+2],]
  # formats output
  cdf[[pstate_list_names[instate]]]=paste(round((cdf)[,1],2),"(",round((cdf)[,3]),"-",round((cdf)[,2]),")",sep  ="")
  
  
  return(cdf[,c(4,5)])
}


# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = transmission
#   model output at the barn level; past_transmission_results_premlevel = another set of transmission model out at the barn 
#   level; outdays = selects days to output results
# Outputs: table with the mean and 90% P.I. of the number of pigs in a selected disease state on the days selected in
#   outdays
Comb_tbl_state_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel,
                              past_transmission_results_premlevel=transmission_results_premlevel,outdays=NA ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # matrices with the number of pigs in the given disease state over time post virus exposure for each iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[instate]]}))
  tmatrix2 <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){past_transmission_results_premlevel[[h]][[instate]]}))
  
  # determines mean and 90% interval for the number pigs in the given disease state over time post virus exposure
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  df2<-data.frame(
    mean_state=apply(tmatrix2,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )

  # formats output
  cdf<-data.frame(DPE=df$times, paste(round(df$mean_state,1),"(",round(df$state_5_percentile),"-",round(df$state_95_percentile),")",sep  =""),paste(round(df2$mean_state,1),"(",df2$state_5_percentile,"-",df2$state_95_percentile,")",sep  =""),row.names = NULL)
  names(cdf)[3]<-paste(pstate_list_names[instate],'previous')
  names(cdf)[2]<-paste(pstate_list_names[instate],'current')
  
  # selects days to print 
  cdf<-cdf[ tseq[outdays/sim_parms$readdt+2],]
  return(cdf)
}


# Inputs: sens_frame = data frame with pen-level prevalences and associated oral fluid PCR sensitivities
# Outputs: plots the oral fluids diagnostic sensitivity for different pen-level infection prevalence values
get_oral_sens_plot<-function(sens_frame){
  # vector with the test sensitivities over pen level prevalence
  sensvec<-sapply(X = seq(0,1,length.out = 100),function(y){get_sen_penprev(inprev = y,sens_frame = sens_frame)})
  df<-data.frame("Within pen prevalence"=seq(0,1,length.out = 100),"Diagnostic Sensitvity"=sensvec)
  colnames(df)<-c("Within pen prevalence","Diagnostic Sensitvity")
  
  # plots diagnostic sensitivity over pen prevalence
  pl<-ggplot(data=df,aes(x=`Within pen prevalence`,y=`Diagnostic Sensitvity`))+geom_line(size=1.2,color='slateblue')
  pl<-pl+theme(axis.title.x=element_text(margin = margin(t = 15)))
  pl<-pl+theme(axis.title.y=element_text(margin = margin(r = 15)))
  pl<-pl+theme(strip.text.y = element_text(size = 15))
  pl<-pl+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  pl<-pl+theme(panel.background = element_rect(fill="ivory"))
  return(pl)
}


# Inputs: survscenstring = dummy variable; in_detcurveoutframe = data frame with the probability
#   of detection for different days post virus exposure and for different active surveillance protocols
# Outputs: prints a plot with the probability of detection over time post virus exposure for different active
#   surveillance protocols
plot_geom_detcurve<-function(survscenstring,in_detcurveoutframe){
  
  pl<-ggplot(data=in_detcurveoutframe,aes(x=testdaypostexposure,y=Detect.percent))
  pl<-pl+geom_line(size=1.2, aes(linetype=survscenario,color=survscenario))
  pl<-pl+labs(x="Time post exposure (Days)")
  pl<-pl+labs(y="Detection probability")
  pl<-pl+scale_color_brewer(type = 'qual',palette = 'Dark2')
  pl<-pl+scale_linetype_discrete()
  pl<-pl+labs(color='Scenario')
  pl<-pl+labs(linetype='Scenario')
  pl<-pl+theme(axis.title.x=element_text(margin = margin(t = 15)))
  pl<-pl+theme(axis.title.y=element_text(margin = margin(r = 15)))
  pl<-pl+theme(strip.text.y = element_text(size = 15))
  pl<-pl+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  pl<-pl+theme(panel.background = element_rect(fill="ivory"))
  pl<-pl+theme(axis.text=element_text(size=12),
               axis.title=element_text(size=14),
               legend.text=element_text(size=12),
               legend.title=element_text(size=14),
               legend.position = "top")
  
  return(pl)
}


# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = transmission
#   model output at the barn level; outdays = days to use for outputting results
# Outputs: a table with the number of latently infected pigs on the days post exposure given by 'outdays'
getdatatable_latent_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel,
                                   outdays=NA ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # matrix with the number of latently infected pigs (latently infected to die + latently infected to recover) over
  # time post virus exposure for each simulation iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[2]]+
    transmission_results_premlevel[[h]][[3]]}))

  # data frame with the mean and 90% P.I. for the number of latently infected pigs  
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    'DPE'=tseq
  )
  
  instate<-3
  cdf<-df
  # picks the results to store
  cdf<-cdf[ tseq[outdays/sim_parms$readdt+2],]
  names(cdf)[1]<-'Latent_pigs'
  names(cdf)[2]<-paste('perc_95th',state_list_names[instate],sep="")
  names(cdf)[3]<-paste('perc_5th',state_list_names[instate],sep="")
  # formats output
  cdf$'Latent pigs'=paste(round(cdf$Latent_pigs,1),"(",round(cdf$perc_95thlatdead),"-",round(cdf$perc_95thlatdead),")",sep  ="")
  
  # selects columns from data fram 'cdf'
  toutdf=as.data.frame((cdf)[,c(4,5)])
  # updates column names
  names(toutdf)[2]<-('Latent,mean (95% prediction interval)')
  
  return((toutdf))
}


# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = transmission
#   model output at the barn level; past_transmission_results_premlevel = another set of transmission model out at the barn 
#   level; outdays = selects days to output results
# Outputs: the number of latently infected pigs at the time points specified by 'outdays' for each transmission model run
Comb_getdatatable_latent_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel,
                                        past_transmission_results_premlevel=transmission_results_premlevel,outdays=NA ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  
  # matrices with the number of latently infected pigs over time post exposure for each simulation iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[2]]+
    transmission_results_premlevel[[h]][[3]]}))
  tmatrix2 <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){past_transmission_results_premlevel[[h]][[2]]+
    past_transmission_results_premlevel[[h]][[3]]}))
  
  # stores the mean and 90% P.I. for the number of latently infected pigs
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  df2<-data.frame(
    mean_state=apply(tmatrix2,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix2,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  df$'Latent pigs current'=paste(round(df$mean_state,1),"(",round(df$state_5_percentile),"-",round(df$state_95_percentile),")",sep  ="")
  df2$'Latent pigs previous'=paste(round(df2$mean_state,1),"(",round(df2$state_5_percentile),"-",round(df2$state_95_percentile),")",sep  ="")
  
  # combines output from each transmission model run into a single data set
  cdf<-data.frame(DPE=df$times,`Latent pigs current` = df$'Latent pigs current',
                  `Latent pigs previous`=df2$'Latent pigs previous',row.names = NULL)
  names(cdf)[3]<-('Latent pigs previous simulation mean (95% prediction interval)')
  names(cdf)[2]<-('Latent pigs current simulation mean (95% prediction interval)')
  
  # selects output to print
  cdf<-cdf[ tseq[outdays/sim_parms$readdt+2],]
  
  return(cdf)
}

