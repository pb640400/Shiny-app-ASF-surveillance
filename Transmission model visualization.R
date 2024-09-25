# This file contains functions for plotting transmission model output in different ways


library(ggplot2)
source('trans_post_process_func.R')


# define the disease states
pstate_list_names<-c("Susceptible","Latent_recovering","Latent_dying","Infectious_recovering","Infectious_dying","Recovered","Dead","Clinical","Blood_positive","Severe_clinical")


# Function prints a plot comparing any number of pens in a specified state and specified iteration
# Inputs: instate = number id for a disease state; initer = iteration number; in_penlist = array of pen number id's;
#   transmission_results = pen-level transmission output for each simulation iteration; in_facet = T/F to make separate
#   panel in the plot for each pen
# Outputs: prints a plot with the number of pigs over time post exposure in the state specified by instate, simulation
#   iteration specified by initer, and pens specified by in_penlist
plot_penstate_iter<-function(
    instate=1,
    initer=1,
    in_penlist=c(1,2),
    transmission_results=transmission_results, 
    in_facet=TRUE
){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  # Input: pen number
  # Output: tranmission model output for the desired iteration, pen, disease state
  get_df_forpen<-function(inpen){
    df<-data.frame(state_plotted=transmission_results[[initer]][[inpen]][[instate]],times=tseq,pen=inpen)
  }
  
  # get number of pigs over time in the desired disease state and simulation iteration for each pen in in_penlist
  cdf<-do.call("rbind",lapply(in_penlist,get_df_forpen))
  names(cdf)[1]<-pstate_list_names[instate]
  
  # plots number of pigs over time in the desired disease state
  pl<-ggplot(data=cdf,aes_string(x=names(cdf)[2],y=names(cdf)[1]))
  pl<-pl+geom_point(size=1.2, aes(color=factor(pen)))
  pl<-pl+labs(x="Time post exposure (Days)")
  pl<-pl+labs(y=paste(  strsplit(names(cdf)[1],split="_")[[1]],sep = " ",collapse = " "))
  if(in_facet==TRUE){
    pl<-pl+facet_grid(pen~.)
  }
  pl<-pl+labs(color='Pen')
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
  
  print(pl)
}



# Function returns matrix for a given state for the pens infected in a given order (e.g. 3rd) in each iteration
# Inputs: instate = number specifying disease state (1=sus; 2=latrec; 3=latdead; 4=infrec; 5=infdead; 6=recovered; 7=dead; 
#   8=clinical; 9=bloodpos; 10=severeclin); pen_num_postexp = order of when a pen got infectious pigs (e.g. if 1, user
#   wants output from the first pig with infectious pigs in each iteration); transmission results = number of pigs over
#   time in each disease state for each pen and each iteration
# Outputs: matrix with the number of pigs over time in the desired disease state and for the pen that had infectious pigs in a
#   certain order relative to the other pens (e.g., first or last) for each iteration
get_penpostexp_state <- function(  instate=1,
                                  pen_num_postexp=1,transmission_results=transmission_results ){

  # makes a matrix with the desired disease state output for a pen infected in a certain order for each iteration
  temp_penpostexp_state <- sapply(1:sim_parms$num_iterations, function(x){
    # gets order of pens from the shortest to the longest time until the pen first contains an infectious pig
    temp_orderlist <- get_infectious_ordpens(initer = x,transmission_results = transmission_results)
    # gets a disease state output for the pen that had infectious pigs in a certain order (e.g. first or second to have infectious pigs)
    temp_state <- (transmission_results[[x]][[(temp_orderlist$peninforder[pen_num_postexp])]][instate][[1]])
    return(temp_state)
  })
  
  #matrix rows are iteration, columns are time step
  return(t(temp_penpostexp_state))
}


# Inputs: initer = iteration number; transmission_results = pen-level transmission simulation output; instate = number
#   specifying a disease state; pen_orderlist = pen orders for when infected (e.g. c(1,5,20) is the 1st, 5th, and 20th 
#   pens with infectious pigs)
# Outputs: a plot with the number of pigs over time in the disease state state specified by instate for pens infected in
#   a certain order as chosen in pen_orderlist
single_iterstate_vispenorder<-function(initer=1,transmission_results = transmission_results,instate = 1,pen_orderlist=c(1,5,20)){
  # gets order of pens by when infectious pigs first appeared in the given iteration
  orderlist<-get_infectious_ordpens(initer = initer,transmission_results=transmission_results)
  # plots the number of pigs over time in the given state and iteration by when pens had infectious pigs
  plot_penstate_iter(instate = instate,initer = initer,in_penlist = orderlist[[1]][pen_orderlist],transmission_results = transmission_results)
}


# Inputs: instate = disease state number; transmission_results = pen-level transmission simulation output; pen_num_list =
#   array with the pens the user wants plotted by the order infected
# Outputs: plot with the number of pigs in instate over time (mean and 90% C.I.) for the pens infected in the order given
#   in pen_num_list (e.g., c(1,2,3) are the 1st, 2nd, and 3rd pens with infectious pigs in each iteration) across all simulation iterations
Visualize_geom_ribbon_state_facet<-function(instate = instate,transmission_results = transmission_results,pen_num_list=c(1,2,3)){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  # Inputs: the pen with infectious pigs appearing in a desired order relative to other pens
  # Outputs: data frame with the mean and 95% P.I. for the number of pigs over time in the pen infected in a certain order
  get_one_penstatedf<-function(in_pen_order){
    # matrix with the number of pigs over time for the inputted disease state for the pen with infectious pigs in a certain order
    # across the iterations
    tmatrix <- get_penpostexp_state(transmission_results = transmission_results,instate = instate,pen_num_postexp = in_pen_order)
    # stores mean and 90% interval over time for the output in tmatrix
    df<-data.frame(
      mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
      state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
      state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
      penorder=in_pen_order,
      times=tseq
    )
  }

  
  # get the mean and 90% C.I. of the number pigs in a given disease state for each pen infected in the order specified
  # in pen_num_list
  cdf<-do.call("rbind",lapply(pen_num_list,get_one_penstatedf))
  names(cdf)[1]<-paste(pstate_list_names[instate],sep="")
  names(cdf)[2]<-paste('perc_95th',pstate_list_names[instate],sep="")
  names(cdf)[3]<-paste('perc_5th',pstate_list_names[instate],sep="")
  names(cdf)[4]<-'Pen_infection_order'
  cdf$Pen_infection_order<-factor(cdf$Pen_infection_order)
  
  # plot the mean and 90% P.I. of the number of pigs in the disease state over time
  pl<-ggplot(data=cdf,aes_string(x=names(cdf)[5],y=names(cdf)[1]))
  pl<-pl+geom_line(size=1.2, aes(color=factor(Pen_infection_order)))
  pl<-pl+facet_grid(Pen_infection_order~.)
  pl<-pl+labs(x="Time post exposure (Days)")
  pl<-pl+labs(color='Infection order of pen')
  pl<-pl+labs(y=paste(  strsplit(names(cdf)[1],split="_")[[1]],sep = " ",collapse = " "))
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
  pl<-pl+geom_ribbon(data = cdf, aes_string(x=names(cdf)[5], ymin=names(cdf)[3],ymax=names(cdf)[2],fill= names(cdf)[4]),alpha=0.25)+
    scale_fill_discrete(name="Infection order of pen")
  print(pl)
    
}


# Inputs: sim_parms = simulation parameters; instate = disease state number; transmission_results_premlevel = barn-level
#   transmission model output
# Outputs: plot of the mean and 90% P.I. of the number of pigs over time in the barn in the disease state specified by 
#   instate across the simulation iterations
Visualize_geom_ribbon_state_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel ){
  
  #  array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  # matrix with the number of pigs in a disease state over time at the barn level for each simulation iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[instate]]}))
      
  # get mean and 95% P.I. from output in tmatrix   
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  # names column variables in df  
  cdf<-df
  names(cdf)[1]<-paste('Mean ',pstate_list_names[instate],sep="")
  names(cdf)[2]<-paste('perc_95th',pstate_list_names[instate],sep="")
  names(cdf)[3]<-paste('perc_5th',pstate_list_names[instate],sep="")
  cdf$dfill='blue'
  cdf$dfill=factor(cdf$dfill)
  
  # plots mean and 90% P.I. for the number of pigs over time in the disease state
  pl<-ggplot(data=cdf,aes_string(x=as.name(names(cdf)[4]),y=as.name(names(cdf)[1]),fill=names(cdf)[5]))
  pl<-pl+scale_linetype_identity(name="", guide='legend', labels=c("Median"))
  pl<-pl+geom_line(size=1.2, aes(linetype="solid"))
  pl<-pl+labs(x="Time post exposure (Days)")
  pl<-pl+labs(color='Infection order of pen')
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
    scale_fill_identity(name="", guide='legend',labels=c("90 percent prediction interval"))
  print(pl)
  
}


# Inputs: transmsision_results = the number of pigs in each disease state over time post exposure for each pen across 
#   the simulation iterations
# Outputs: plots the mean and 90% P.I. (as estimated across the simulation iterations) for the number of infectious pens 
#   over time post virus exposure
Visualize_geom_ribbon_ninfpens<-function(transmission_results=transmission_results ){
  
  # array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  # gets the number of pens with the number of infectious pigs exceeding 0 for each simulation time point across the iterations
  # (get_num_infectious_pen_iter defined in trans_post_process_func.R)
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){get_num_infectious_pens_iter(initer = h,transmission_results = transmission_results)}))
  
  
  # gets mean and 90% P.I. of output in tmatrix
  df<-data.frame(
    mean_infpens=apply(tmatrix,2,mean,na.rm=TRUE),
    Infectiouspens_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95,type=3),
    Infectiouspens_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05,type=3),
    times=tseq
  )
  # further processing of data frame
  cdf<-df
  names(cdf)[1]<-paste('mean_','infectious_pens',sep="")
  names(cdf)[2]<-paste('perc_95th','infectious_pens',sep="")
  names(cdf)[3]<-paste('perc_5th','infectious_pens',sep="")
  cdf$dfill='blue'
  cdf$dfill=factor(cdf$dfill)
  
  
  # plot mean and 90% P.I.
  pl<-ggplot(data=cdf,aes(x=times,y=mean_infectious_pens,fill='blue'))
  pl<-pl+scale_linetype_identity(name="", guide='legend', labels=c("Median"))
  pl<-pl+geom_line(size=1.2, aes(linetype="solid"))
  pl<-pl+labs(x="Time post exposure (Days)")
  pl<-pl+labs(color='Infection order of pen')
  pl<-pl+theme(axis.title.x=element_text(margin = margin(t = 15)))
  pl<-pl+theme(axis.title.y=element_text(margin = margin(r = 15)))
  pl<-pl+theme(strip.text.y = element_text(size = 15))
  pl<-pl+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  pl<-pl+theme(panel.background = element_rect(fill="ivory"))
  pl<-pl+labs(y="Number of infectious pens")
  pl<-pl+theme(axis.text=element_text(size=12),
               axis.title=element_text(size=14),
               legend.text=element_text(size=12),
               legend.title=element_text(size=14),
               legend.position = "top")
  dfill<-'dodgerblue2'
  pl<-pl+geom_ribbon(data = cdf, aes(x=times, ymin=perc_5thinfectious_pens,ymax=perc_95thinfectious_pens,fill='blue'),alpha=0.25)+
    scale_fill_identity(name="", guide='legend',labels=c("90 percent prediction interval"))
  print(pl)
  
}


# Inputs: initer = simulation iteration number; instate = number ID of disease state; transmission_results_premlevel =
#   the number of pigs in each disease state over time post virus exposure at the barn-level
# Outputs: plot with the number of pigs over time post exposure at the barn level in the chosen iteration and disease state
Visualize_iter_state_prem<-function(initer=initer,instate = instate,transmission_results_premlevel=transmission_results_premlevel ){
  # array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  # number of pigs over time in the barn for the desired iteration and disease state
  tmatrix <- transmission_results_premlevel[[initer]][[instate]]
  
  # stores tmatrix output in a data frame
  df<-data.frame(
    mean_state=tmatrix,
    times=tseq
  )
  cdf<-df
  names(cdf)[1]<-paste(state_list_names[instate],sep="")
  
  
  # plot number of pigs over time in the barn in the disease state
  pl<-ggplot(data=cdf,aes_string(x=names(cdf)[2],y=names(cdf)[1]))
  pl<-pl+geom_line(size=1.2, linetype="solid", col="brown")
  pl<-pl+labs(x="Time post exposure (Days)")

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

  print(pl)
  
}


# Inputs: sim_parms = simulation parameters; instate = disease state number ID; transmission_results_premlevel = number
#   of pigs over time in each disease state at the barn-level
# Outputs: plot of the number of infectious pigs over time post exposure (mean and 90% P.I.)
Visualize_geom_ribbon_infectious_prem<-function(sim_parms,instate = instate,transmission_results_premlevel=transmission_results_premlevel ){
  # array with numbers from 0 to sim_parm$ndays by sim_parms$readdt time step
  tseq<-seq(0,sim_parms$ndays/sim_parms$time_step,by=sim_parms$readdt/sim_parms$time_step)*sim_parms$time_step
  # matrix with the number of infectious pigs in the barn over time for each iteration
  tmatrix <- do.call('rbind',lapply(1:sim_parms$num_iterations,function(h){transmission_results_premlevel[[h]][[4]]+
    transmission_results_premlevel[[h]][[5]]}))
  
  
  # store mean and 90% P.I. in data frame
  df<-data.frame(
    mean_state=apply(tmatrix,2,mean,na.rm=TRUE),
    state_95_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.95),
    state_5_percentile=apply(tmatrix,2,quantile,na.rm=TRUE,probs=0.05),
    times=tseq
  )
  instate<-4
  cdf<-df
  names(cdf)[1]<-'Infectious_pigs'
  names(cdf)[2]<-paste('perc_95th',state_list_names[instate],sep="")
  names(cdf)[3]<-paste('perc_5th',state_list_names[instate],sep="")
  cdf$dfill='blue'
  cdf$dfill=factor(cdf$dfill)
  
  
  # plot mean and 90% P.I. for the number of infectious pigs over time
  pl<-ggplot(data=cdf,aes_string(x=names(cdf)[4],y=names(cdf)[1],fill=names(cdf)[5]))
  pl<-pl+scale_linetype_identity(name="", guide='legend', labels=c("Median"))
  pl<-pl+geom_line(size=1.2, aes(linetype="solid"))
  #pl<-pl+facet_grid(Pen_infection_order~.)
  pl<-pl+labs(x="Time post exposure (Days)")
  pl<-pl+labs(y="Infectious pigs")
  pl<-pl+labs(color='Infection order of pen')
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
    scale_fill_identity(name="", guide='legend',labels=c("90 percent prediction interval"))
  print(pl)
  
}
