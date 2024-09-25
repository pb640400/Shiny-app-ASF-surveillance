# This file contains functions for simulating the application of active surveillance protocols and storing the desired
# output in a data frame. There are two sets of functions, one considering a random movement day and one for running
# a detection curve (probability of detection over increasing days of virus exposure relative to movement)


# Random movement day functions

# This function, given transmission model output, simulates the application of active surveillance protocols and returns
# the results
# Inputs: updateprogress = function for updating progress bar in Shiny app; premovsurvlist = active surveillance protocols;
#   trans_shin_list = list with the cumulative barn-level transmission results, the pen-level transmission results, simulation parameters,
#   simulated routine mortality, and simulated routine morbidity; inscenariosvec = active surveillance protocols to run;
#   in_sim_parms = simulation parameters; pensizearray = number of pigs per pen; sens_frame = data frame with pen level
#   prevalences and associated oral fluid sensitivities; contact_rate_string = name of the contact rate scenario;
#   strain_scenario_string = name of ASFV strain scenario
# Outputs: a data frame with the simulated active surveillance model outputs as well as information on the 
#   active surveillance protocol scenarios and simulation parameters used
run_and_process_randday_premovsurv_shiny <- function(updateprogress,premovsurvlist=premovsurvlist,trans_shin_list,
                                                     inscenariosvec,  in_sim_parms = in_sim_parms,pensizearray,sens_frame,
                                                     contact_rate_string,strain_scenario_string){
  

  
  # simulates the desired active surveillance protocols given in 'inscenariosvec'
  temp_surv_res <- run_trans_and_premovesurv_scenarios_shiny(updateprogress=updateprogress,in_sim_parms = in_sim_parms,
                                                             inscenariosvec = inscenariosvec,trans_shin_list =trans_shin_list,
                                                             premovsurvlist = premovsurvlist,pensizearray=pensizearray,sens_frame=sens_frame)
  
  

  
  # processes active surveillance model output and information and stores it in a data frame
  tout<-postprocess_premovsurv(strain_beta_surv_res = temp_surv_res,sim_parms =in_sim_parms,premovsurvlist =premovsurvlist,
                               strain_scenario_string =strain_scenario_string,
                               contact_rate_string = contact_rate_string,inscenariosvec=inscenariosvec)
  
  return( tout)
}



# This function simulates the application of given active surveillance protocols and returns a list of the results for
# each simulated protocol
# Inputs: updateprogress = function for updating progress bar in shiny app; in_sim_parms = simulation parameters;
#   inscenariosvec = active surveillance protocols to run; trans_shin_list = list with the cumulative barn-level 
#   transmission results, the pen-level transmission results, simulation parameters, simulated routine mortality, 
#   and simulated routine morbidity; premovsurvlist = list with the active surveillance protocol details; pensizearry =
#   number of pigs per pen; sens_frame = data frame with pen level prevalences and associated oral fluid sensitivities
# Outputs: list containing the simulated outcomes from the active surveillance model 
run_trans_and_premovesurv_scenarios_shiny<-function(updateprogress,in_sim_parms = in_sim_parms,
                                                    inscenariosvec,
                                                    trans_shin_list=trans_shin_list,
                                                    premovsurvlist=premovsurvlist,
                                                    pensizearray,sens_frame
                                                    ){
  
  
  # Simulates the active surveillance protocols in 'inscenariosvec'. The function perform_all_iter_premovsurv is
  # defined in the file 'run_surv_scenario.R'
  surv_out_list<-lapply(inscenariosvec,function(k,updateprogress){
    updateprogress(value=((which(inscenariosvec==k)    -.5)/length(inscenariosvec)*100), detail="Simulating...")
    perform_all_iter_premovsurv(in_scenario = k,transmission_results=trans_shin_list$transmission_results,
                                sim_normmort=trans_shin_list$sim_normmort,
                                sim_normsick=trans_shin_list$sim_normsick,
                                transmission_results_premlevel=trans_shin_list$transmission_results_premlevel,
                                ignore_dieouts=sim_parms$ignoredieouts,
                                boolretnumsamp = sim_parms$boolretnumsamp,premovsurvlist=premovsurvlist,
                                pensizearray=pensizearray,in_sim_parms = in_sim_parms,sens_frame=sens_frame)
  },updateprogress=updateprogress)
  
  
  return(surv_out_list)
}



# This function gets transmission model output for a given ASFV strain scenario and transmission rate scenario. Routine
# mortality and morbidity are also simulated.
# Input: updateprogress = function for updating progress bar in Shiny app; in_sim_parms = simulation parameters;
#   in_strain = strain scenario; in_beta_choice = pre-defined transmission scenario; inbetavec = user-defined
#   within-pen transmission rate (if relevant); in_barn_pen_arraylist = list of arrays with 1) the room of each pen, 
#   2) whether a pen is an edge pen or not, 3) which pens are adjacent, 4) which rooms are adjacent;  inbetaBbvec =
#   user-defined between-pen transmission rate (if relevant)
# Output: list with the cumulative barn-level transmission results, the pen-level transmission results, simulation parameters,
#   simulated routine mortality, and simulated routine morbidity
get_single_trans_results_progressbar<-function(updateprogress,in_sim_parms = sim_parms,in_strain=1,in_beta_choice=1,inbetavec=NA,in_barn_pen_arraylist,inbetaBbvec=NA){
  # load transmission model function written in C. File type needed depends on computer type.
  dyn.load("transmission_asf_debug_gammaall.dll")
  #dyn.load("transmission_asf_debug_gammaall.so")
  
  stime<-Sys.time()

  # updates disease state parameters according to strain choice
  in_sim_parms <- sim_sate_parms_adjust(in_st_parmchoice = in_strain, in_sim_parms)

  # Checks whether transmission rate scenario is user defined or preset
  if(in_beta_choice==4){
    # transmission parameters updated according to user defined values
    in_sim_parms <- sim_beta_pars(in_beta_choice = 3, in_sim_parms)
    in_sim_parms$betapars=inbetavec
    in_sim_parms$betabetween=inbetaBbvec
  }else{
    # transmission parameters updated according to pre-defined scenarios
    in_sim_parms <- sim_beta_pars(in_beta_choice = in_beta_choice, in_sim_parms)
  }   
  
  # simulates a within pen transmission rate for each simulation iteration
  beta_arr_init<-rspert(n = in_sim_parms$num_iterations,x.min = in_sim_parms$betapars[1],x.mode =in_sim_parms$betapars[2],x.max = in_sim_parms$betapars[3] )
  # simulates a between room transmission rate for each simulation iteration
  beta_roomarr<-rspert(n = in_sim_parms$num_iterations,x.min = in_sim_parms$betaroom[1],x.mode =in_sim_parms$betaroom[2],x.max = in_sim_parms$betaroom[3] )
  # array repeating the barn size for each simulation iteration
  barnsize_arr <- rep(in_sim_parms$flock_size_pars, in_sim_parms$num_iterations)
  
  # when the transmission scenario number is 7, adjusts the within and between pen transmission rates for distance independent spread
  if(in_sim_parms$trans_type==7) {
    # simulates between pen transmission rate for each sumulation iteration
    beta_betweenarr_initial<-runif(n = in_sim_parms$num_iterations,min = in_sim_parms$betabetween[1],max = in_sim_parms$betabetween[3] )
    # adjusts between pen transmission rate for proportion of between pen spread due to distance independent pathways
    beta_betweenarr<-beta_betweenarr_initial*(1-in_sim_parms$between_penpeoplefrac)#for nose
    # determines transmission rate for spread related to distance independent pathways
    beta_peoplearr<-(beta_betweenarr_initial-beta_betweenarr)/(1-1/(in_sim_parms$npens/in_sim_parms$nrooms))
    # adjusts within pen transmission for transmission due to distance independent pathways
    beta_arr<-beta_arr_init-beta_peoplearr*1/(in_sim_parms$npens/in_sim_parms$nrooms)
  }
  
  

  # code to facilitate progress bar for transmission model runs
  splits_forbar=20
  temp_translist=list()
  temp_splititerarray=rep(floor(in_sim_parms$num_iterations/splits_forbar),splits_forbar)
  temp_split_iterstorun=list()
  temp_splititerarray[splits_forbar]=temp_splititerarray[splits_forbar]+(in_sim_parms$num_iterations-
                                                                        floor(in_sim_parms$num_iterations/splits_forbar)*splits_forbar)              
  
  counter=0   
  for(i in 1:splits_forbar){ 
    temp_split_iterstorun[[i]]=(counter+1):(counter+temp_splititerarray[i])
    counter=counter+  temp_splititerarray[i]
  }
  
  
  for(i in 1:splits_forbar){
    # simulates virus transmission. Call to transmission function depends on transmission scenario. Simulated data
    # stored as a list per iteration.
    if(in_sim_parms$trans_type==7){
      updateprogress(value=(i-1)*temp_splititerarray[1]/in_sim_parms$num_iterations*100, detail="Simulating...")
      temp_translist[[i]] <- lapply((temp_split_iterstorun[[i]]), function(x){get_single_sim(betawithin = beta_arr[x], betabetween = beta_betweenarr[x], barnsize_arr[x],betaroom =beta_roomarr[x],
                                                                                             sim_parms = in_sim_parms,betapeople =beta_peoplearr[x],in_barn_pen_arraylist = in_barn_pen_arraylist  )})
    }else{
      temp_translist[[i]] <- lapply(temp_split_iterstorun[[i]], function(x){get_single_sim(beta_arr[x], beta_betweenarr[x],
                                                                                           barnsize_arr[x],betaroom =beta_roomarr[x],sim_parms = in_sim_parms,in_barn_pen_arraylist = in_barn_pen_arraylist )})
    }
    
  }
  
  # reformats transmission model output
  transmission_results=do.call("c",temp_translist)

  
  ## the statement adds up the pen level results to get aggregate results at premises level for all iterations
  transmission_results_premlevel<-lapply(1:in_sim_parms$num_iterations,
                                         function(g){summarize_prem_oneiter(transmission_results[[g]])})
  
 
  # array with the number of pigs per pens
  tpensizearray<-rep(round(in_sim_parms$flock_size_pars/in_sim_parms$npens),in_sim_parms$npens)
  
  # simulates normal (routine) mortality. Function defined in file Blood_detect_exp.R
  sim_normmort <- sim_normmort_func(sim_parms  =in_sim_parms,pensizearray = tpensizearray )
  
  # simulates pigs with clinical signs that would occur normally during routine production. Function defined in file Blood_detect_exp.R
  sim_normsick <- sim_normsick_func(transmission_results = transmission_results, in_sim_parms = in_sim_parms,pensizearray = tpensizearray)
  

  # collects aggregate transmission results, pen-level transmission results, simulation parameters, simulated normal mortality,
  #   and simulated routine pigs with clinical signs into a list
  trans_outlist=(list(transmission_results_premlevel=transmission_results_premlevel,transmission_results=transmission_results,
                      in_sim_parms=in_sim_parms,sim_normmort=sim_normmort,sim_normsick=sim_normsick))

  # unloads transmission simulation function written in C
 dyn.unload("transmission_asf_debug_gammaall.dll")

  return(trans_outlist)
  
}



# This function processes active surveilance model information and output and returns it in a data frame
# Inputs: strain_beta_surv_res = list containing output from active surveillance simulation model; sim_parms = simulation
#   parameters; premovsurvlist = active surveillance protocol details; contact_rate_string = name of transmission model
#   scenario; strain_scenario_string = name of ASFV strain scenario; inscenariosvec = vector selecting active surveillance
#   protocols
# Outputs: data frame with active surveillance simulation model output and active surveillance protocol details
postprocess_premovsurv <- function(strain_beta_surv_res,sim_parms,premovsurvlist,contact_rate_string,strain_scenario_string,
                                   inscenariosvec){
  
  # sampling prioritizations for blood/tissue and oral fluids sampling
  bloodsurv_samptypelist=c(" dead_sick_live", "sick_dead_live", "dead_live","dead only", "dead_sick","sick_live","random live")
  OFsurv_samptypelist =c("random", "sick", "dead", "dead_sick")
      
  
  # fills data frame with simulated active surveillance model output, simulation parameter values, and active surveillance
  # protocol details
  tout<-as.data.frame( do.call("rbind",lapply(1:length(strain_beta_surv_res),function(k){
    
    toutin<- data.frame(
        "Detect percent"=strain_beta_surv_res[[k]]$detect_mvec,
        "Ind_pcr_dtect"=strain_beta_surv_res[[k]]$blood_res_mvec,
        "Rope_det_perc"=strain_beta_surv_res[[k]]$rope_res_mvec,
        "Mortdetperc"=strain_beta_surv_res[[k]]$mort_res_mvec,
        "mort_sum_iterations"=strain_beta_surv_res[[k]]$mort_summ_mat,
        "mortdetfrac"=strain_beta_surv_res[[k]]$mort_avg_mat,
        "meanmoveday"=strain_beta_surv_res[[k]]$move_day_mvec,
        "infnotdetectdel1"=strain_beta_surv_res[[k]]$inf_not_det[1],
        "infnotdetectdel2"=strain_beta_surv_res[[k]]$inf_not_det[2],
        "infnotdetectdel3"=strain_beta_surv_res[[k]]$inf_not_det[3],
        "infnotdetectdel4"=strain_beta_surv_res[[k]]$inf_not_det[4],
        "infnot2.5del1"=strain_beta_surv_res[[k]]$inf_not_det_vec_95int[1,1],
        "infnot2.5del2"=strain_beta_surv_res[[k]]$inf_not_det_vec_95int[1,2],
        "infnot2.5del3"=strain_beta_surv_res[[k]]$inf_not_det_vec_95int[1,3],
        "infnot2.5del4"=strain_beta_surv_res[[k]]$inf_not_det_vec_95int[1,4],
        "infnot97.5del1"=strain_beta_surv_res[[k]]$inf_not_det_vec_95int[2,1],
        "infnot97.5del2"=strain_beta_surv_res[[k]]$inf_not_det_vec_95int[2,2],
        "infnot97.5del3"=strain_beta_surv_res[[k]]$inf_not_det_vec_95int[2,3],
        "infnot97.5del4"=strain_beta_surv_res[[k]]$inf_not_det_vec_95int[2,4],
        "infnotdetectmeandel1"=strain_beta_surv_res[[k]]$inf_not_det_vec_mean[1],
        "infnotdetectmeandel2"=strain_beta_surv_res[[k]]$inf_not_det_vec_mean[2],
        "infnotdetectmeandel3"=strain_beta_surv_res[[k]]$inf_not_det_vec_mean[3],
        "infnotdetectmeandel4"=strain_beta_surv_res[[k]]$inf_not_det_vec_mean[4],
        
        "meandeadsamples"=strain_beta_surv_res[[k]]$meannumsampresvec[1],
        "meansicksamples"=strain_beta_surv_res[[k]]$meannumsampresvec[2],
        "meanlivesamples"=strain_beta_surv_res[[k]]$meannumsampresvec[3],row.names = NULL
        
      )
      
    
    toutin$samplesize<-paste(premovsurvlist[[inscenariosvec[k]]]$indsampsizes,collapse = ",")
    toutin$premisesconfig<-sim_parms$prem_config
    toutin$barnsize<-sim_parms$flock_size_pars
    toutin$nrooms<-sim_parms$nrooms
    toutin$peoplespread<-sim_parms$between_penpeoplefrac
    toutin$refrig<-premovsurvlist[[inscenariosvec[k]]]$refrigyesno
    toutin$survscenario<-premovsurvlist[[inscenariosvec[k]]]$scenarioname
    toutin$bloodsurvtype<-bloodsurv_samptypelist[[premovsurvlist [[inscenariosvec[k]]]$bloottargetype]]
    toutin$OFsurvtype<-OFsurv_samptypelist[[premovsurvlist[[inscenariosvec[k]]]$OFtargetype]]
    toutin
  })))
  
  tout$Contactrate<-contact_rate_string
  tout$Strain<-strain_scenario_string
  

  # rearranges output
  outframe<-tout[,c((ncol(tout)-0:4),1:(ncol(tout)-5))]
  outframe<-outframe[,-(17:28)]

  
  return(outframe)
}




###############################################################################################################################
###############################################################################################################################

# Detection curve function


# This function simulates active surveillance model results for given active surveillance protocols and days of virus
# exposure prior to movement
# Inputs: updateprogress = function for updating progress bar in Shiny app; in_sim_parms = simulation parameters;
#   min_day = closest day prior to movement when virus exposure can occur; max_day = furthest day prior to movement when
#   virus exposure can occur; by_day = interval by which to simulate active surveillance for days of virus exposure from min_day to max_day;
#   inscenariosvec = vector with the active surveillance protocols to run; trans_shin_list = list with the cumulative 
#   barn-level transmission results, the pen-level transmission results, simulation parameters, simulated routine mortality, 
#   and simulated routine morbidity; premovsurvlist = list with details on the active surveillance protocols; pensizearray =
#   number of pigs per pen; sens_frame = data frame with pen level prevalences and associated oral fluid test sensitivities
# Outputs: list of active surveillance model output for each protocol and desired time of virus exposure prior to movement
run_trans_and_premovesurv_detcurve_shiny<-function(updateprogress,in_sim_parms = in_sim_parms,min_day = 1, max_day = 25,
                                                   by_day=1, inscenariosvec, trans_shin_list=trans_shin_list, 
                                                   premovsurvlist=premovsurvlist, pensizearray,sens_frame){
  
  # array of days prior to movement when virus exposure occurs
  tmovedayseq<-seq(min_day,max_day,by=by_day)
  
  # iterates through active surveillance protocols to evaluate
  tdet_survout_list<-lapply(inscenariosvec,function(k,updateprogress){
    
    # iterates through array with days of virus exposure
    lapply(tmovedayseq,function(ttt,updateprogress){
      
      # makes copy of simulation parameters
      tinsimparms<- in_sim_parms
      # only allows virus exposure to occur ttt days prior to movement
      tinsimparms$move_day_bounds<-c(ttt,ttt)
      
      # update progress bar on app
      tcounter=(which(inscenariosvec==k)-0.5)*length(tmovedayseq)+which(tmovedayseq==ttt)
      updateprogress(value=100*tcounter/(length(inscenariosvec)*length(tmovedayseq)), detail="Simulating...")
      
      # simulate application of active surveillance protocol. The function 'perform_all_iter_premovsurv' is defined in the
      # the file 'run_surv_scenario.R'
      perform_all_iter_premovsurv(in_scenario = k,transmission_results=trans_shin_list$transmission_results,
                                  sim_normmort=trans_shin_list$sim_normmort,
                                  sim_normsick=trans_shin_list$sim_normsick,
                                  transmission_results_premlevel=trans_shin_list$transmission_results_premlevel,
                                  ignore_dieouts=tinsimparms$ignoredieouts,
                                  boolretnumsamp = tinsimparms$boolretnumsamp,premovsurvlist=premovsurvlist,pensizearray=pensizearray,in_sim_parms = tinsimparms,sens_frame=sens_frame
      )

    },updateprogress =updateprogress)
    
  },updateprogress =updateprogress)
  
  
  return(tdet_survout_list)
}



# This function stores detection curve active surveillance model output
# Inputs: det_curv_out_list = list containing active surveillance model output; in_sim_parms = simulation parameters;
#   premovsurvlist = active surveillance protocol details; contact_rate_string = name of transmission scenario;
#   strain_scenario_string = name of ASFV strain scenario; inscenariosvec = active surveillance protocols that were run;
#   min_day = time of virus exposure to consider that is closest to the time of movement; max_day = time of virus exposure
#   to consider that is furthest from the time of movement; by_day = time step to use between min_day and max_day
# Outputs: data frame containing active surveillance simulation model detection curve results, surveillance protocol 
#   details, and values for key simulation parameters
postprocess_detcurv_shiny <- function(det_curv_out_list,in_sim_parms,premovsurvlist,contact_rate_string,
                                      strain_scenario_string,inscenariosvec,min_day = 1, max_day = 25,by_day){
  
  # sampling prioritizations for blood/tissue and oral fluids sampling
  bloodsurv_samptypelist=c(" dead_sick_live", "sick_dead_live", "dead_live","dead only", "dead_sick","sick_live","random live")
  OFsurv_samptypelist =c("random", "sick", "dead", "dead_sick")

  #exposure as relates to testing instead of movement (e.g. if testing occurs one day prior to movement and 
  # exposure occurs one day prior to movement then exposure occurs at the same time as testing)
  test_day_post_exp=seq(min_day,max_day,by=by_day) - 1 
  
  # iterate through active surveillance simulation results
  outframe<-as.data.frame(do.call("rbind",lapply(1:length(det_curv_out_list), function(k){
    
    # iterate through days of virus exposure
    toutframe <- do.call("rbind",lapply(1:length(det_curv_out_list[[k]]), function(x){
      
      # stores active surveillance simulation model output, protocol details, and key simulation parameters
      data.frame(
        
        "Detect percent"=det_curv_out_list[[k]][[x]]$detect_mvec,
        "Ind_pcr_dtect"=det_curv_out_list[[k]][[x]]$blood_res_mvec,
        "Rope_det_perc"=det_curv_out_list[[k]][[x]]$rope_res_mvec,
        "Mortdetperc"=det_curv_out_list[[k]][[x]]$mort_res_mvec,
        "mort_sum_iterations"=det_curv_out_list[[k]][[x]]$mort_summ_mat,
        "mortdetfrac"=det_curv_out_list[[k]][[x]]$mort_avg_mat,
        "meanmoveday"=det_curv_out_list[[k]][[x]]$move_day_mvec,
        "infnotdetectdel1"=det_curv_out_list[[k]][[x]]$inf_not_det[1],
        "infnotdetectdel2"=det_curv_out_list[[k]][[x]]$inf_not_det[2],
        "infnotdetectdel3"=det_curv_out_list[[k]][[x]]$inf_not_det[3],
        "infnotdetectdel4"=det_curv_out_list[[k]][[x]]$inf_not_det[4],
        
        "meandeadsamples"=det_curv_out_list[[k]][[x]]$meannumsampresvec[1],
        "meansicksamples"=det_curv_out_list[[k]][[x]]$meannumsampresvec[2],
        "meanlivesamples"=det_curv_out_list[[k]][[x]]$meannumsampresvec[3],
        
        "samplesize"=paste(premovsurvlist[[inscenariosvec[k]]]$indsampsizes,collapse = ","),
        "premisesconfig"=in_sim_parms$prem_config,
        "barnsize"=in_sim_parms$flock_size_pars,
        "nrooms"=in_sim_parms$nrooms,
        "peoplespread"=in_sim_parms$between_penpeoplefrac,
        "testdaypostexposure"= test_day_post_exp[x],
        "refrig"=premovsurvlist[[inscenariosvec[k]]]$refrigyesno,
        "survscenario"=premovsurvlist[[inscenariosvec[k]]]$scenarioname,
        "bloodsurvtype"=bloodsurv_samptypelist[[premovsurvlist[[inscenariosvec[k]]]$bloottargetype]],
        "OFsurvtype"=OFsurv_samptypelist[[premovsurvlist[[inscenariosvec[k]]]$OFtargetype]],
        "contactrate"=contact_rate_string,
        "strain"=strain_scenario_string,
        "survscenarionum"=inscenariosvec[k],
        
        row.names = NULL
      )
    }))
    
  })))
  

 # rearranges data output
  outframe<-outframe[,c((ncol(outframe)-0:6),1:(ncol(outframe)-7))]
  outframe<-outframe[,c(27,8,1:7,9:26)]
 
  
  return(outframe)
}



# Function that runs the active surveillance simulation model and stores the desired output
# Inputs: updateprogress = function for updating progress bar in Shiny app; in_sim_parms = simulation parameters;
#   min_day = time of virus exposure to consider that is closest to the time of movement; max_day = time of virus exposure
#   to consider that is furthest from the time of movement; by_day = time step to use between min_day and max_day;
#   premovsurvlist = list of active surveillance protocol details; trans_shin_list = list with the cumulative 
#   barn-level transmission results, the pen-level transmission results, simulation parameters, simulated routine mortality, 
#   and simulated routine morbidity; inscenariosvec = active surveillance protocols ran; in_sim_parms = simulation parameters;
#   pensizearray = number of pigs per pen; sens_frame = data frame with pen-level prevalences and associated sensitivities
#   of testing by oral fluids; contact_rate_string = name of transmission scenario; strain_scenario_string = name of 
#   transmission scenario
# Outputs: data frame with output from the active surveillance simulation model, key simulation parameters, and
#   details of the active surveillance protocol
run_and_process_premovsurv_detcurve_shiny <- function(updateprogress =updateProgress,min_day=1, max_day=25, by_day=1,
                                                 premovsurvlist=premovsurvlist,trans_shin_list,inscenariosvec,  in_sim_parms = in_sim_parms,
                                                 pensizearray,sens_frame,contact_rate_string,strain_scenario_string){
  
  # simulates application of active surveillance considering days of virus exposure from min_day to max_day discretized
  # by by_day
  det_curv_out_single <- run_trans_and_premovesurv_detcurve_shiny(updateprogress =updateprogress,in_sim_parms = in_sim_parms, inscenariosvec = inscenariosvec,trans_shin_list = trans_shin_list,premovsurvlist = premovsurvlist,pensizearray =pensizearray,sens_frame = sens_frame, 
                                                            min_day = min_day, max_day = max_day,by_day=by_day)
  

  
  # collects in a data frame output from the active surveillance simulation model, key simulation parameters, and
  # details of the active surveillance protocol
  outdetframe<-postprocess_detcurv_shiny(det_curv_out_list = det_curv_out_single,in_sim_parms =in_sim_parms,
                            premovsurvlist =premovsurvlist,
                            contact_rate_string =contact_rate_string,strain_scenario_string = strain_scenario_string,inscenariosvec = inscenariosvec,min_day = min_day,max_day =  max_day,by_day=by_day)

  return(outdetframe)
  
}  
  
  
  