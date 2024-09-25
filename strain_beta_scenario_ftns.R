# This file declares functions for changing the disease state and contact rate parameters according to preset scenarios

# This function updates the disease state parameters given an ASF strain scenario
# Inputs: in_st_parmchoice = the ASF strain scenario; sim_parms = simulation parameters
# Outputs: simulation parameters with updated disease state parameters
sim_sate_parms_adjust<-function(in_st_parmchoice=state_parm_set_list$olsen, sim_parms){
  out_sim_parms<-sim_parms
  if(in_st_parmchoice==state_parm_set_list$olsen){
    out_sim_parms$is_georgia_strain=FALSE
    out_sim_parms$mort_lat_pars = c(12.39518,0.3682945)
    out_sim_parms$recov_lat = c(12.39518,0.3682945)#
    out_sim_parms$mort_inf_pars = c(2.516493,0.9012873)
    out_sim_parms$recov_inf = c(2.516493,0.9012873)
    out_sim_parms$p_mort=1
    out_sim_parms$bpospar=c(-0.92, 0.71)
    out_sim_parms$sero_pars = c(8.5646158,0.6127723)
    out_sim_parms$clinpar=c(30,1)# artificially high so that pigs are clinical till death
    out_sim_parms$p_sero=1 #Fraction of pigs that develop mild clinical signs following exposure
  }else if(in_st_parmchoice==state_parm_set_list$hu) {
    out_sim_parms$is_georgia_strain=TRUE # for simulating a different mean time to clinical signs per simulation run
    out_sim_parms$georgia_clinmean = c(4.985, 1.2463) #mean and standard deviation for normal distribution for mean time to clincial signs Georgia strain
    out_sim_parms$georgia_clinmean_bounds = c(3, 7.5) #min/max bounds for normal distributed mean time to clinical signs Georgia strain
    out_sim_parms$clinpar=c(30,1)# artificially high so that they are clinical till death
    out_sim_parms$p_mort=1
    out_sim_parms$mort_inf_pars = c(22.7, 0.403) #gamma shape and scale from Hu et al.
    out_sim_parms$recov_inf = c(22.7, 0.403)
    out_sim_parms$mort_lat_pars = c(17.25, 0.31) #gamma shape and scale derived from Guinat et al. (time to oral positive)
    out_sim_parms$recov_lat = c(17.25, 0.31)
    out_sim_parms$bpospar=c(-0.6, 3.38) #mean and variance normal distribution from Guitat et al. (pcr blood positive - oral positive; see georgia_clinsigns.R for details)
    out_sim_parms$p_sero=1 #Fraction of pigs that develop mild clinical signs following exposure
  }else if(in_st_parmchoice==state_parm_set_list$olsensens){
    #sensitivity analysis mean 4.5 infectious period
    out_sim_parms$is_georgia_strain=FALSE
    out_sim_parms$mort_lat_pars = c(12.39518,0.3682945)
    out_sim_parms$recov_lat = c(12.39518,0.3682945)#
    out_sim_parms$mort_inf_pars = c(2.516493,1.788203)
    out_sim_parms$recov_inf = c(2.516493,1.788203)
    out_sim_parms$p_mort=1
    out_sim_parms$bpospar=c(-0.92, 0.71)
    out_sim_parms$sero_pars = c(8.5646158,0.6127723)  
    out_sim_parms$clinpar=c(30,1)# artificially high so that they are clinical till death
    out_sim_parms$p_sero=1 #Fraction of pigs that develop mild clinical signs following exposure
  }else if(in_st_parmchoice==state_parm_set_list$Modvir_severeclin){
    #moderately virulent with mild clinical parameters replaced by severe clinical parameters
    out_sim_parms$is_georgia_strain=FALSE
    out_sim_parms$mort_lat_pars = c(13.299,0.3384482) #c(shape, scale) gamma
    out_sim_parms$recov_lat = c(13.299,0.3384482)
    out_sim_parms$mort_inf_pars = c(9.632,0.862) #c(shape, scale) gamma
    out_sim_parms$recov_inf =  c(55.42012,0.7950162)
    out_sim_parms$p_mort=0.4 #0.1
    out_sim_parms$bpospar=c(-0.82, 0.82)
    out_sim_parms$sero_pars = c(41.969,0.259) #shape and scale parameter for time to onset of severe clinical signs
    out_sim_parms$clinpar=c(1.027,6.515)#duration of severe clinical signs
    out_sim_parms$p_sero=0.688 #Fraction of pigs that develop SEVERE clinical signs following exposure
  } else{
    #moderately virulent
    out_sim_parms$is_georgia_strain=FALSE
    out_sim_parms$mort_lat_pars = c(13.299,0.3384482) #c(shape, scale) gamma
    out_sim_parms$recov_lat = c(13.299,0.3384482)
    out_sim_parms$mort_inf_pars = c(9.632,0.862) #c(shape, scale) gamma
    out_sim_parms$recov_inf =  c(55.42012,0.7950162)
    out_sim_parms$p_mort=0.4
    out_sim_parms$bpospar=c(-0.82, 0.82)
    out_sim_parms$sero_pars = c(26.257,0.214)
    out_sim_parms$clinpar=c(3.418,3.2)
    out_sim_parms$p_sero=1 #Fraction of pigs that develop mild clinical signs following exposure
  }
  return(out_sim_parms)
}



# This function updates the disease transmission parameters given a transmission scenario
# Inputs: in_beta_choice = the ASF transmission scenario; sim_parms = simulation parameters
# Outputs: simulation parameters with updated transmission parameters
sim_beta_pars <- function(in_beta_choice, sim_parms){
  out_sim_parms <- sim_parms
  if(in_beta_choice == beta_set_list$georgia_slow){
    #from Guinat et al. (2016)
    out_sim_parms$betapars=c(0.3,0.6,1.0)
    out_sim_parms$betabetween=c(0.1,0.3,0.5)
    out_sim_parms$betaroom=c(0.001,0.05,0.15)
    
  } else if(in_beta_choice == beta_set_list$georgia_fast){
    #from Hu et al. (2017)
    out_sim_parms$betapars=c(0.96,2.62,5.61)
    out_sim_parms$betabetween=c(0.31,0.99,1.98)
    out_sim_parms$betaroom=c(0.001,0.05,0.15)
  } else if(in_beta_choice == beta_set_list$Moderately_virulent){
    #within-pen from Malladi et al.; between pen from Guinat et al. (2016)
    out_sim_parms$betapars=c(1.00,1.64,2.74)
    out_sim_parms$betabetween=c(0.1,0.3,0.5)
    out_sim_parms$betaroom=c(0.001,0.05,0.15)
  }
  return(out_sim_parms)
}


