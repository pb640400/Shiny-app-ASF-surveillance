#This file contains functions for running surveillance scenarios
#assumes transmission runs have been performed and pen and premises level transmission are available




# This function simulates active surveillance for a given transmission simulation iteration and active surveillance
# protocol
# Inputs: initer = transmission model simulation iteration number; in_scenario = number specifying active surveillance
#   protocol; transmission_results = number of pigs in each disease state over time per pen for each iteration;
#   sim_normmort = simulated routine mortality; sim_normsick = simulated routine morbidity; transmission_results_premlevel =
#   total number of pigs (barn-level) in each disease state over time; in_sim_parms = simulation parameters; boolretnumsamp =
#   T/F whether to return average number of dead/sick/live pigs sampled as part of blood/tissue PCR sampling; premovsurvlist =
#   specifies active surveillance protocols; pensizearry = number of pigs per pen; sens_frame = datapoints with PCR sensivity
#   at different pen level prevalences for oral fluids testing; in_randmoveday = the day prior to movement when exposure
#   to ASFV occurs
# Outputs: list with variables blood_res_iter = 1 or 0 if detection by blood/tissue testing; rope_res_iter = 1 or 0
#   if detection by oral fluids testing; mort_res_iter = 1 or 0 if detection by mortality trigger; mortsum = 
#   number of times detection occurred by mortality trigger; rand_move_day_int = the time of exposure prior to movement;
#   numpsampvec = number of samples taken from dead/sick/live pigs for blood/tissue testing
perform_one_iterpremovesurv<-function(initer,in_scenario=1,transmission_results=transmission_results,
                                      sim_normmort=sim_normmort,
                                      sim_normsick=sim_normsick,
                                      transmission_results_premlevel=transmission_results_premlevel,
                                      in_sim_parms=in_sim_parms,
                                      boolretnumsamp=FALSE,
                                      premovsurvlist=premovsurvlist,pensizearray=pensizearray,sens_frame,in_randmoveday){
  
 
  # store day prior to movement when exposure to ASF occurs in new variable
  rand_move_day_int=in_randmoveday

  
  # Determines number of rope samples tested
  if(is.na(premovsurvlist[[in_scenario]]$ropetestdays[1])){
  temp_npensropes_sampled=NA
  }else {
  temp_npensropes_sampled=round(premovsurvlist[[in_scenario]]$ropepenpercs*in_sim_parms$npens/100)
    
  }
  
  # determines day post exposure when ropes tested
  if(is.na(premovsurvlist[[in_scenario]]$ropetestdays[1])){
    ropetesttimesarr<-NA
    
  }else{
    
    if ( length(premovsurvlist[[in_scenario]]$ropetestdays)>0){
      ropetesttimesarr<-rand_move_day_int-premovsurvlist[[in_scenario]]$ropetestdays
      ropetesttimesarr<- ropetesttimesarr[ropetesttimesarr>0]
      ropetesttimesarr<-  ropetesttimesarr[length(ropetesttimesarr):1]
    }
  }

  # determines day post exposure when blood/tissues tested with associated sample sizes
  if(!is.na(premovsurvlist[[in_scenario]]$temp_movevec[1])){
    if ( length(premovsurvlist[[in_scenario]]$temp_movevec)>0){
      pcrtesttimesarr<-rand_move_day_int-premovsurvlist[[in_scenario]]$temp_movevec
      pcrsamplesizesarr<-premovsurvlist[[in_scenario]]$indsampsizes
      pcrsamplesizesarr<-pcrsamplesizesarr[pcrtesttimesarr>0]
      pcrtesttimesarr<- pcrtesttimesarr[pcrtesttimesarr>0]
      
      #reordering earlier first
      if(length(pcrtesttimesarr)>1){
      pcrtesttimesarr<-pcrtesttimesarr[length(pcrtesttimesarr):1]
      pcrsamplesizesarr<-pcrsamplesizesarr[length(pcrtesttimesarr):1]  
      }  
      if(length(pcrtesttimesarr) == 0){
        pcrtesttimesarr <- NA
        pcrsamplesizesarr <- premovsurvlist[[in_scenario]]$indsampsizes[1]
      }
      
    }
  } else {
    pcrtesttimesarr <- NA
  }
  

  # determines days post exposure when mortality trigger applied
  if(!is.na(premovsurvlist[[in_scenario]]$morttrigpriorday[1])){
    if ( length(premovsurvlist[[in_scenario]]$morttrigpriorday)>0){
      if(premovsurvlist[[in_scenario]]$ismorttrigconsecutive){
        morttesttimesarr<-rand_move_day_int-(1:premovsurvlist[[in_scenario]]$morttrigpriorday)
      } else {
        morttesttimesarr<-rand_move_day_int-premovsurvlist[[in_scenario]]$morttrigpriorday
      }
      morttesttimesarr<-  morttesttimesarr[ morttesttimesarr>0]
      morttesttimesarr<-  morttesttimesarr[length( morttesttimesarr):1]
      
      if(length(morttesttimesarr) == 0){morttesttimesarr <- NA}
    }
  } else {
    morttesttimesarr <- NA
  }
  

  
  # Simulates mild sick pigs to be missed/unobserved
  # sim_normsick matrix and transmission_results changed directly when oral fluid and individual pig blood samples are taken. Done this way to keep consistency
  # between sampling types without having to pass new arrays into the functions.
  testtimes <- unique(c(ropetesttimesarr, pcrtesttimesarr))
  for(i in 1:length(testtimes)){
    
    if(!is.na(testtimes[i])){
      temp_time <- testtimes[i]/in_sim_parms$readdt+1
      # get normal and disease sick for each pen at time period (functions from Blood_detect_exp.R)
      pennormsick <- sapply(1:in_sim_parms$npens, function(x){get_pen_normsick(time_test_readtindex = temp_time, initer = initer, inpen = x,sim_normsick = sim_normsick)})
      pendissick <- sapply(1:in_sim_parms$npens, function(x){get_pen_mildclin(time_test_readtindex = temp_time, initer = initer, inpen = x,transmission_results = transmission_results)})
      
      
      # Sim sick pigs to be missed
      # sim number sick pigs noticed
      num_sickobs <- rbinom(1, size = sum(pendissick, pennormsick), prob = in_sim_parms$propsickobs)
      # divide total sick observed between actual sick and normal sick
      num_dissickobs <- rhyper(1, m=sum(pendissick), n=sum(pennormsick), k=num_sickobs)
      num_normsickobs <- num_sickobs - num_dissickobs
      # divide number asf sick and normal sick observed between pens
      pendissick <- rmvhyper(nn=1, n=pendissick, k=num_dissickobs)
      pennormsick <- rmvhyper(nn=1, n=pennormsick, k=num_normsickobs)
      
      # Directly change sim_normsick and transmission_results
      sim_normsick[[initer]][,temp_time] <- pennormsick
      for(j in 1:in_sim_parms$npens){
        transmission_results[[initer]][[j]]$clinical[temp_time] <- pendissick[j]
      }
    }
    
  }
  
 
  # Checks if there is rope testing
  if(is.na(premovsurvlist[[in_scenario]]$ropetestdays[1])){
    # if no rope testing, whether detection by oral fluids occurred and the pens sampled are set to NA
    roperes<-NA; ropesampled_pens_arr<-NA
  }else{
    
    # initialize matrix for pens sampled for rope testing  
    ropesampled_pens_arr<-matrix(NA,nrow=length(ropetesttimesarr),ncol = temp_npensropes_sampled)
      
    # initialize rope test results
    roperes<-list(NA)
    
    # loop through rope testing times
    for(g in 1:length(ropetesttimesarr)){
      # translate rope testing time in day post exposure to simulation time point
      if(is.na(ropetesttimesarr[g])){
        ttime_test_readtindex<-1
      }else{
        ttime_test_readtindex<-(ropetesttimesarr[g]/in_sim_parms$readdt+1)
      }
        
      trepeatflag<-FALSE
        # can repeat testing of pens previously tested if not first test time, if repeatropepens is not NA and valuefor ifrepeatpens =TRUE
      if (g!=1 & !is.na(premovsurvlist[[in_scenario]]$If_repeat_ropepens)){
        if(premovsurvlist[[in_scenario]]$If_repeat_ropepens==TRUE) {
          trepeatflag<-TRUE
        }
      }
      
      # simulates oral fluid sampling  
      if (!trepeatflag){
        tempout<-sample_oral_iter_simtime(time_test_readtindex = (ttime_test_readtindex),initer = initer,npens_sampled=temp_npensropes_sampled,sim_normmort=sim_normmort,sim_normsick=sim_normsick,
                                          transmission_results=transmission_results,in_OFsampletype=premovsurvlist[[in_scenario]]$OFtargetype,pensizearray=pensizearray,in_sim_parms = in_sim_parms,sens_frame=sens_frame)
      }else {
        # this is where previous ropepens sampled have to be used
        tempout<-sample_oral_iter_simtime(time_test_readtindex = (ttime_test_readtindex),initer = initer,npens_sampled=temp_npensropes_sampled,sim_normmort=sim_normmort,sim_normsick=sim_normsick,
                                          transmission_results=transmission_results,in_OFsampletype=premovsurvlist[[in_scenario]]$OFtargetype,in_selectedpensvec = ropesampled_pens_arr[g-1,] ,pensizearray=pensizearray,in_sim_parms =in_sim_parms,sens_frame=sens_frame  )
      }
      ropesampled_pens_arr[g,]<-tempout$sampled_pens
      roperes[[g]]<-tempout
    }
      
    #vector with 0 or 1 for detection on each sampling day
    roperes <- sapply(1:length(roperes), function(x){roperes[[x]]$ropedet})
  }# end if loop for whether to do rope sampling
  

  
  # checks if there are legitimate blood/tissue sample times (e.g. samples after exposure)
  if(is.na(pcrtesttimesarr[1])){
    # if no blood/tissue PCR test times, pigs will be sampled at first time period which will give negative test result
    pcrtesttimesarr<-1
    pcrreslist<-lapply(1:length(pcrtesttimesarr),function(g){
      sample_blood_iter_simtime(initer = initer,time_test_readtindex = 1,
                                sample_type =premovsurvlist[[in_scenario]]$bloottargetype,
                                insamplesize = 1,transmission_results=transmission_results,
                                sim_normmort=sim_normmort,
                                sim_normsick=sim_normsick,
                                transmission_results_premlevel=transmission_results_premlevel,boolretnumsamp = boolretnumsamp,pensizearray=pensizearray,in_sim_parms =in_sim_parms )
    })
    pcrreslist[[1]]$is_pcr_detected<-0
    pcrreslist[[1]]$N_pospcR_tubes<-0
    pcrreslist[[1]]$outnumsamplesvec[]<-0
  }else{
    # converts PCR test times from day post exposure to simulation time point
    pcrreslist<-lapply(1:length(pcrtesttimesarr),function(g){
      if(is.na(pcrtesttimesarr[g])){
        ttime_test_readtindex<-1
      }else{
        ttime_test_readtindex<-(pcrtesttimesarr[g]/in_sim_parms$readdt+1)
      }
    
    
      # conditions when we prioritize pcr pens based on rope sampling
      Tflag_prioritizepens<-FALSE
      if(!is.na(premovsurvlist[[in_scenario]]$if_premovpen_same_rope)){
        if(premovsurvlist[[in_scenario]]$if_premovpen_same_rope==TRUE){
          if (!is.na(ropesampled_pens_arr[1,1])){
            if(!is.na(ropetesttimesarr[1])){
              
              if (pcrtesttimesarr[g]>=min(ropetesttimesarr)){# ropesampling has been performed
                Tflag_prioritizepens<-TRUE
              } 
            }
          } 
        }
      }
   
   
      # checks if we should test same pens by PCR that were rope tested       
      if (Tflag_prioritizepens){
        # checks whether PCR sampling includes samples across multiple days that were regrigerated
        if((premovsurvlist[[in_scenario]]$refrigyesno==FALSE)){
          
          # indicates pens tested by oral fluids sampling that should also be sampled for blood/tissue PCR
          tselectedpens<-ropesampled_pens_arr[which.max(ropetesttimesarr[ropetesttimesarr<=pcrtesttimesarr[g]]),]
          # simulates blood/tissue PCR testing
          tres<-sample_blood_iter_simtime(initer = initer,time_test_readtindex = ttime_test_readtindex,
                                          sample_type =premovsurvlist[[in_scenario]]$bloottargetype,
                                          insamplesize = pcrsamplesizesarr[g],transmission_results=transmission_results,
                                          sim_normmort=sim_normmort,
                                          sim_normsick=sim_normsick,
                                          transmission_results_premlevel=transmission_results_premlevel,boolretnumsamp = boolretnumsamp,
                                          in_selectedpensvec = tselectedpens,pensizearray=pensizearray,in_sim_parms = in_sim_parms
          )
          
        }else{
          tres<-sample_blood_iter_simtime_refri(initer = initer,time_test_readtindex = ttime_test_readtindex,
                                                sample_type =premovsurvlist[[in_scenario]]$bloottargetype,
                                                insamplesize = pcrsamplesizesarr[g],transmission_results=transmission_results,
                                                sim_normmort=sim_normmort,
                                                sim_normsick=sim_normsick,
                                                transmission_results_premlevel=transmission_results_premlevel,instordays = premovsurvlist[[in_scenario]]$refdaysaccum,inmaxdeadsamples =premovsurvlist[[in_scenario]]$refrimaxdead,
                                                boolretnumsamp = boolretnumsamp,
                                                in_selectedpensvec = tselectedpens,
                                                pensizearray=pensizearray,in_sim_parms = in_sim_parms) 
          
        }
        
      } else { 
        # this is the main scenario where pcr pens are not prioritized according to ropes
        if((premovsurvlist[[in_scenario]]$refrigyesno==FALSE)){
          # simulates blood/tissue PCR testing
          tres<-sample_blood_iter_simtime(initer = initer,time_test_readtindex = ttime_test_readtindex,
                                          sample_type =premovsurvlist[[in_scenario]]$bloottargetype,
                                          insamplesize = pcrsamplesizesarr[g],transmission_results=transmission_results,
                                          sim_normmort=sim_normmort,
                                          sim_normsick=sim_normsick,
                                          transmission_results_premlevel=transmission_results_premlevel,boolretnumsamp = boolretnumsamp,pensizearray=pensizearray,in_sim_parms = in_sim_parms)
          
        }else{
          # simulates blood/tissue PCR testing when samples are from multiple days
          tres<-sample_blood_iter_simtime_refri(initer = initer,time_test_readtindex = ttime_test_readtindex,
                                                sample_type =premovsurvlist[[in_scenario]]$bloottargetype,
                                                insamplesize = pcrsamplesizesarr[g],transmission_results=transmission_results,
                                                sim_normmort=sim_normmort,
                                                sim_normsick=sim_normsick,
                                                transmission_results_premlevel=transmission_results_premlevel,instordays = premovsurvlist[[in_scenario]]$refdaysaccum,inmaxdeadsamples =premovsurvlist[[in_scenario]]$refrimaxdead,boolretnumsamp = boolretnumsamp,pensizearray=pensizearray,in_sim_parms=in_sim_parms) 
          
        } 
        
      }
    
      
    })# brackets for sapply for PCR results list
  }
  
  # 1 or 0 for whether detected or not each pcr test day
  pcrres<-do.call("c",lapply(1:length(pcrreslist),function(tt){pcrreslist[[tt]][[1]] }))
  # number dead, sick, and random live pigs sampled
  numpsampvec<-NA
  if(boolretnumsamp==TRUE){
    numpsampmat<- do.call("rbind",lapply(1:length(pcrreslist),function(tt){pcrreslist[[tt]][[3]] }))
    numpsampvec<-apply( numpsampmat,2,sum,na.rm=TRUE)
  }
  
 
  # simulates application of mortality trigger returning 1 if detected and 0 if not
  mortres<-sapply(1:length(morttesttimesarr),function(g){
    if(is.na(morttesttimesarr[g])){
      return(0)
    }else{
      return(perform_mort_surv_iter(time_test_readtindex =(morttesttimesarr[g]/in_sim_parms$readdt+1),initer = initer,transmission_results = transmission_results,sim_normmort = sim_normmort,transmission_results_premlevel = transmission_results_premlevel,in_sim_parms = in_sim_parms,
                             inmorttrig = premovsurvlist[[in_scenario]]$morttrigfrac,
                             tnumalive = (sum(pensizearray)-transmission_results_premlevel[[initer]]$dead[(morttesttimesarr[g]/in_sim_parms$readdt+1)]))) 
    }
 
  })
  
  # overall indicator if detection occurred by blood/tissue sampling
  if(sum(pcrres)>0){
    blood_res_iter<-1
  }else{
    blood_res_iter<-0
  }
  
  # overall indicator if detection occurred by rope sampling
  if(!is.na(roperes[1])){
    if(sum(roperes)>0){
      rope_res_iter<-1
    }else{
      rope_res_iter<-0
    }
  }else{rope_res_iter<-NA}
  
  # overall indicator if detection occurred by mortality trigger
  if(sum(mortres)>0){
    mort_res_iter<-1
  }else{
    mort_res_iter<-0
  }
  
  
  return(list(blood_res_iter=blood_res_iter,rope_res_iter=rope_res_iter,mort_res_iter=mort_res_iter,mortsum=sum(mortres),
              rand_move_day_int=rand_move_day_int,numpsampvec=numpsampvec))
  
}



# This function simulates an active surveillance protocol for each transmission model iteration
# Inputs: in_scenario = number indicating the active surveillance protocol; transmission_results = number of pigs in each
#   disease state over time at the pen level; sim_normmort = simulated routine mortality; sim_normsick = simulated routine
#   morbidty; transmission_results_premlevel = total number of pigs in each disease state over time; in_sim_parms = simulation
#   parameters; ignore_dieouts = boolean variable for whether transmission iterations with maximum mortality less than
#   the threshold value should be ignored; boolretnumsamp = boolean variable for whether the average number of dead/sick/random live
#   pigs sampled for blood/tissue testing should be included in the output; premovsurvlist = active surveillance protocol inputs;
#   pensizearray = array with the number of pigs per pen; sens_frame = datapoints with PCR sensitivity per pen level
#   prevalence for oral fluids testing
# Outputs: list containing detectmvec = overall probability of detection; blood_res_mvec = probability of detection for
#   blood/tissue sampling; rope_res_mvec = probability of detection for oral fluids testing; mort_res_mvec = probability
#   of detection for mortality trigger; mort_summ_mat = number of days mortality trigger exceeded per iteration;
#   mort_avg_mat = average number of days mortality trigger exceeded per iteration; move_day_mvec = average day of
#   exposure prior to movement; inf_not_det = number of infectious but not detected pigs for different days post movement;
#   inf_not_det_vec_95int = 95% interval for infectious but not detected pigs; inf_not_det_vec_mean = mean number of
#   infectious but not detected pigs; meannumsampresvec = mean number of dead/sick/random live pigs sampled as part
#   of blood/tissue PCR testing
perform_all_iter_premovsurv<-function(in_scenario=1,transmission_results=transmission_results,
                                        sim_normmort=sim_normmort,
                                        sim_normsick=sim_normsick,
                                        transmission_results_premlevel=transmission_results_premlevel,
                                        in_sim_parms=in_sim_parms,
                                        ignore_dieouts=TRUE,boolretnumsamp=TRUE,
                                      premovsurvlist=premovsurvlist,pensizearray,sens_frame){
  
  
  # Inputs: in_ttweights = an array with the likelihood of exposure for each day prior to movement of interest;
  #   tnumsims = the number of transmission runs meeting die-out criteria (if relevant); in_sim_parms = simulation parameters
  # Outputs: array with the exposure day to use for each simulation of active surveillance 
  get_rand_move_arr<-function(in_ttweights,tnumsims, in_sim_parms){
    # divides up exposure times between the iterations to run dependent on the probability of exposure on a given day prior to movement
    titerdaysarr<-floor(in_ttweights/sum(in_ttweights)*tnumsims)
    extraiter=tnumsims-sum(titerdaysarr)
    for(k in 1:extraiter){
      ttt_day=sample(x = 1:length(titerdaysarr),size = 1,prob = in_ttweights)
      titerdaysarr[ ttt_day]=titerdaysarr[ ttt_day]+1
    }
    return(do.call("c",lapply(1:length(titerdaysarr),function(x){rep((in_sim_parms$move_day_bounds[1]:in_sim_parms$move_day_bounds[2])[x],titerdaysarr[x])}) ))
  }
  
  
  
  # Checks whether transmission run iterations where the infection dies-out should be ignored
  if(ignore_dieouts){
    
    # If infection die outs are ignored, calculates array indicating which transmission runs do not meet the threshold
    # for total number dead used to define a die-out
    is_dieout <- sapply(1:in_sim_parms$num_iterations, function(g){
      if(max(transmission_results_premlevel[[g]]$dead) <= in_sim_parms$dieout_mortthresh){
        return(1)
      } else {return(0)}
    })
    
    # array with the iterations where die-out condition is not met
    nodieoutindex <- which(is_dieout == 0)
  }else{
    # if infection die-outs are not ignored, all simulation iterations are kept
    nodieoutindex <-1:in_sim_parms$num_iterations
  }
    
   
  # Checks whether there is at least one day of PMIP
  if(in_sim_parms$PMIPdur<1){
    # if there is no PMIP, each day between move_day_bound[1] and move_day_bound[2] has an equal likelihood of virus exposure
    tweights<-rep(1,(in_sim_parms$move_day_bound[2]-in_sim_parms$move_day_bound[1]+1)) 
    # gets an array with the day of virus exposure prior to movement to use for each simulation of active surveillance
    randmovedayarr=get_rand_move_arr(in_ttweights =tweights,tnumsims =length(nodieoutindex),in_sim_parms = in_sim_parms  )
    
  }else{
    # if there is a PMIP, the days between move_day_bound[1] and move_day_bound[2] that are part of the PMIP have a
    # reduced likelihood of virus exposure dependent on PMIPdur
    tweights<-rep(1,(in_sim_parms$move_day_bound[2]-in_sim_parms$move_day_bound[1]+1)) 
    #to adjust for non zero start days. 
    ttindeces<-((1:in_sim_parms$PMIPdur)-in_sim_parms$move_day_bound[1]+1)
    ttindeces<-ttindeces[ttindeces>0]
    tweights[ttindeces]<-(1-in_sim_parms$PMIPeff)
    # gets an array with the day of virus exposure prior to movement to use for each simulation of active surveillance
    # according to its relative likelihood of exposure
    randmovedayarr=get_rand_move_arr(in_ttweights =tweights,tnumsims =length(nodieoutindex),in_sim_parms = in_sim_parms  ) 
    
  }
    
    
  # performs surveillance on iterations where spread did not die out
  ttesreslist <- lapply(1:length(nodieoutindex), function(g){perform_one_iterpremovesurv(initer=nodieoutindex[g],in_scenario=in_scenario,transmission_results=transmission_results,
                                                                                         sim_normmort=sim_normmort,
                                                                                         sim_normsick=sim_normsick,
                                                                                         transmission_results_premlevel=transmission_results_premlevel,boolretnumsamp = boolretnumsamp,premovsurvlist=premovsurvlist,
                                                                                         pensizearray=pensizearray,in_sim_parms = in_sim_parms,sens_frame =sens_frame,in_randmoveday = randmovedayarr[g]  )})
  
  
  ##### Processes active surveillance simulation results #####
    
  #delays for infectious pigs HARD CODE
  tdelarr<-c(0,2,3,4)
 
  
  # Gets the number of infectious but undetected pigs for each simulation iteration at each time point in tdelarr
  # relative to movement (time of movement, 2 days post movement, etc.)
  inf_not_det_mat <-  t( sapply(1:length(ttesreslist), function(h){
    #if not detected return number infectious pigs at time of movement, otherwise return NA
    # bug corrected here. ttrreslist does not include all iterations but no dieout iterations.Therefore when calling
    # accessing simulation results we need to use nodieoutindex
    # because we run an 80 day simulation there is no issue with respect to the additional delays
    inf_not_det_<-sapply(tdelarr, function(k){
    # checks if rope testing was part of active surveillance protocol
    if(is.na(premovsurvlist[[in_scenario]]$ropetestdays[1])){
      # checks if detection occured by blood/tissue testing or mortality trigger
      if(ttesreslist[[h]]$blood_res_iter + ttesreslist[[h]]$mort_res_iter == 0){
        # checks if infectious + latently infected pigs should be recorded or just infectious
        if(in_sim_parms$ret_infected_vs_infectious){
          #if the above is true return infected+infectious pigs
          num_inf <- transmission_results_premlevel[[nodieoutindex[h]]]$infdead[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ] +
            transmission_results_premlevel[[nodieoutindex[h]]]$infrec[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ]+
            transmission_results_premlevel[[nodieoutindex[h]]]$latdead[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ] +
            transmission_results_premlevel[[nodieoutindex[h]]]$latrec[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ]
        }else{
          num_inf <- transmission_results_premlevel[[nodieoutindex[h]]]$infdead[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ] +
            transmission_results_premlevel[[nodieoutindex[h]]]$infrec[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ]
        }
        
        return(num_inf)
      } else{
        return(NA)}
        
    }else{
      # if oral fluids testing performed, check if detection occurred by blood/tissue PCR testing, mortality trigger
      # of oral fluids PCR testing
      if(ttesreslist[[h]]$blood_res_iter + ttesreslist[[h]]$rope_res_iter + ttesreslist[[h]]$mort_res_iter == 0){
        if(in_sim_parms$ret_infected_vs_infectious){
          #if the above is true return infected+infectious pigs
          num_inf <- transmission_results_premlevel[[nodieoutindex[h]]]$infdead[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ] +
            transmission_results_premlevel[[nodieoutindex[h]]]$infrec[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ]+
            transmission_results_premlevel[[nodieoutindex[h]]]$latdead[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ] +
            transmission_results_premlevel[[nodieoutindex[h]]]$latrec[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ]
        }else{
          num_inf <- transmission_results_premlevel[[nodieoutindex[h]]]$infdead[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ] +
            transmission_results_premlevel[[nodieoutindex[h]]]$infrec[ (ttesreslist[[h]]$rand_move_day_int+k)/in_sim_parms$readdt + 1 ]
        }
        return(num_inf)
      } else{
        return(NA)
      }
    }
  })
  
  }))
   
  
  # format infectious undetected pigs as mean (95% prediction interval)
  inf_not_det_vec_mean <- apply(inf_not_det_mat,2,mean, na.rm=TRUE)
  inf_not_det_vec_95int <-apply(inf_not_det_mat,2, quantile,probs = c(0.025, 0.975), na.rm=TRUE) 
  inf_not_det_format <- sapply(1:length(tdelarr),function(kk){
    paste(round(inf_not_det_vec_mean[kk],0), "(", round(inf_not_det_vec_95int[1,kk],0), ",", round(inf_not_det_vec_95int[2,kk],0), ")",sep="")})
  

  # Records whether detected by blood/tissue PCR testing, oral fluids PCR testing, or mortality trigger per iteration
  blood_res_mat<-do.call("c",lapply(1:length(ttesreslist),function(h){ttesreslist[[h]]$blood_res_iter}))
  rope_res_mat<-do.call("c",lapply(1:length(ttesreslist),function(h){ttesreslist[[h]]$rope_res_iter}))
  mort_res_mat<-do.call("c",lapply(1:length(ttesreslist),function(h){ttesreslist[[h]]$mort_res_iter}))
  # Records day of exposure prior to movement
  move_day_mat<-do.call("c",lapply(1:length(ttesreslist),function(h){ttesreslist[[h]]$rand_move_day_int}))
 
  meannumsampresvec<-NA
  if(boolretnumsamp==TRUE){
    # records number of dead, sick, or random live pigs sampled per iteration if this data is desired
    meannumsampresvec<-apply(do.call("rbind",lapply(1:length(ttesreslist),function(h){ttesreslist[[h]]$numpsampvec})),2,mean,na.rm=TRUE)
  }
  
  
  # collects number of times mortality trigger exceeded per iteration
  mortsum_mat<-do.call("c",lapply(1:length(ttesreslist),function(h){ttesreslist[[h]]$mortsum}))
  # probability detection by blood/tissue PCR
  blood_res_mvec<-mean(blood_res_mat)
  # probability detection by oral fluids sampling
  if(is.na(premovsurvlist[[in_scenario]]$ropetestdays[1])){
    rope_res_mvec<-NA
  }else{
    rope_res_mvec<-mean(rope_res_mat,na.rm=TRUE)
  }
  # probability detection by mortality trigger
  mort_res_mvec<-mean(mort_res_mat,na.rm=TRUE)
  mort_summ_mat<-sum(mortsum_mat)
  # average number of times mortality trigger exceeded per iteration
  mort_avg_mat<-ifelse(length(premovsurvlist[[in_scenario]]$morttrigpriorday) > 0, mort_summ_mat/(length(ttesreslist)*length(premovsurvlist[[in_scenario]]$morttrigpriorday)),0)
  # overall probability of detection
  if(is.na(premovsurvlist[[in_scenario]]$ropetestdays[1])){
   detectedvec<- blood_res_mat+mort_res_mat
 }else{
   detectedvec<- blood_res_mat+mort_res_mat+rope_res_mat
 }
  detect_mvec<-mean(detectedvec>0)
  # average day of exposure prior to movement
  move_day_mvec<-mean(move_day_mat)
  
  
  return(list(detect_mvec=detect_mvec,blood_res_mvec=blood_res_mvec,rope_res_mvec=rope_res_mvec, mort_res_mvec= mort_res_mvec,
              mort_summ_mat=mort_summ_mat, mort_avg_mat=mort_avg_mat, move_day_mvec = move_day_mvec, inf_not_det=inf_not_det_format,
              inf_not_det_vec_95int=inf_not_det_vec_95int,
              inf_not_det_vec_mean =inf_not_det_vec_mean ,
              meannumsampresvec=meannumsampresvec))
}



##### These functions actually used??

robus_surv_loop<-function(in_scenario=1){
  out<-tryCatch({
    
    perform_all_iter_test_timeseq(in_scenario)
    
  }, error = function(err) {
    print(paste ("for surv scenario no", in_scenario,err))
    return(-1)
    
  }, finally = {
    
    
    
  })
  return(out)
}


perform_all_iter_test_timeseq<-function(in_scenario=1,transmission_results=transmission_results,
                                        sim_normmort=sim_normmort,
                                        transmission_results_premlevel=transmission_results_premlevel,
                                        ignore_dieouts=TRUE,in_sim_parms=in_sim_parms,sens_frame){
  #browser()
  #testtimesarr<-seq(1,12)/sim_parms$readdt + 1 #+1 to account for 1st array element being time zero
  testtimesarr<-seq(in_sim_parms$nsurvdays[1],in_sim_parms$nsurvdays[2])/in_sim_parms$readdt + 1
  
  if(ignore_dieouts){
    
    is_dieout <- sapply(1:in_sim_parms$num_iterations, function(g){
      if(max(transmission_results_premlevel[[g]]$dead) <= in_sim_parms$dieout_mortthresh){
        return(1)
      } else {return(0)}
    })
    # iterations where dieout condition is not met
    nodieoutindex <- which(is_dieout == 0)
    # performs surveillance on iterations where spread did not die out
    ttesreslist <- lapply(1:length(nodieoutindex), function(g){per_oneiter_test_timesequance(initer = nodieoutindex[g],in_scenario = in_scenario,testtimesarr = testtimesarr,transmission_results=transmission_results,
                                                                                             sim_normmort=sim_normmort,
                                                                                             transmission_results_premlevel=transmission_results_premlevel,in_sim_parms = in_sim_parms,sens_frame = sens_frame)})
    
  } else {
    ttesreslist<-lapply(1:in_sim_parms$num_iterations,function(g){per_oneiter_test_timesequance(initer = g,in_scenario = in_scenario,testtimesarr = testtimesarr,transmission_results=transmission_results,
                                                                                                sim_normmort=sim_normmort,
                                                                                                transmission_results_premlevel=transmission_results_premlevel,in_sim_parms = in_sim_parms,sens_frame=sens_frame )})
  }
  
  
  blood_res_mat<-do.call("rbind",lapply(1:length(ttesreslist),function(h){ttesreslist[[h]]$blood_res_iter}))
  rope_res_mat<-do.call("rbind",lapply(1:length(ttesreslist),function(h){ttesreslist[[h]]$rope_res_iter}))
  blood_res_mvec<-apply(blood_res_mat,2,mean)
  rope_res_mvec<-apply(rope_res_mat,2,mean)
  return(list(blood_res_mvec=blood_res_mvec,rope_res_mvec=rope_res_mvec))
}



per_oneiter_test_timesequance<-function(initer,in_scenario=5,testtimesarr,transmission_results=transmission_results,
                                        sim_normmort=sim_normmort,
                                        transmission_results_premlevel=transmission_results_premlevel,in_sim_parms=in_sim_parms,sens_frame){
  
  temp_npensropes_sampled=round(scendefframe$Percent_pens_sampled_rope[in_scenario]*in_sim_parms$npens/100)
  
  blood_res_iter<-sapply(testtimesarr,function(x){sample_blood_iter_simtime(initer = initer,time_test_readtindex = x,
                                                                            sample_type =scendefframe$blood_surv_type[in_scenario],
                                                                            insamplesize = scendefframe$Blood_surv_nsamples[in_scenario],transmission_results=transmission_results,
                                                                            sim_normmort=sim_normmort,
                                                                            transmission_results_premlevel=transmission_results_premlevel )})
  rope_res_iter<-sapply(testtimesarr,function(x){sample_oral_iter_simtime(time_test_readtindex = x,initer = initer,npens_sampled=temp_npensropes_sampled,transmission_results=transmission_results,pensizearray=pensizearray,in_sim_parms = in_sim_parms,sens_frame=sens_frame )})
  return(list(blood_res_iter=blood_res_iter,rope_res_iter=rope_res_iter))
}


