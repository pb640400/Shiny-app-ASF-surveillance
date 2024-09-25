# This file contains a function for re-formatting the active surveillance protocol details read in from a CSV file
# for use in the simulation models. It also contains a function for checking that the variables in an output of
# get_premov_survlist have the correct format.


# Input: data frame with active surveillance protocol details
# Output: active surveillance protocol details formatted in list form for simulation runs 
get_premov_survlist<-function(premovsurframe){
  # initialize list to store active surveillance details
  premovsurvlist<-as.list(rep(NA,nrow(premovsurframe)))
  # loop through active surveillance protocols
  for (k in 1:nrow(premovsurframe)){
    # store active surveillance protocol name
    scenarioname<-premovsurframe$Scenario_Label[k]
    
    # store vector with days prior to movement to perform blood and/or tissue PCR testing
    temp_movevec<-as.numeric(strsplit(as.character(premovsurframe$Premove.survdays[k]),split = ",")[[1]])
    if(!is.numeric(temp_movevec)){
      temp_movevec<-NA 
    }else if(length(temp_movevec)==0){
      temp_movevec=NA
    }
    
    # store vector with days prior to movement to perform oral fluids sampling
    if(!is.na(premovsurframe$Oral.surv.days[k]) ){
      if(premovsurframe$Oral.surv.days[k]!="" & tolower(premovsurframe$Oral.surv.days[k])!="na") { 
        ropetestdays<-as.numeric(strsplit(as.character(premovsurframe$Oral.surv.days[k]),split = ",")[[1]])
      }else{
        ropetestdays<-NA
      }
    }else{
      ropetestdays<-NA
    }
    
    # stores percent of total pens to sample for oral fluids
    if(!is.na(premovsurframe$Ropesamppenprev[k])){
      ropepenpercs<-as.numeric(strsplit(as.character(premovsurframe$Ropesamppenprev[k]),split = ",")[[1]])
    }else{
      ropepenpercs<-NA
    }
    
    # stores  days prior to movement to mortality trigger is applied
    if(!is.na(premovsurframe$Mortriggdaysprior[k])){
      morttrigpriorday=as.numeric(strsplit(as.character(premovsurframe$Mortriggdaysprior[k]),split=",")[[1]])
    }else{
      morttrigpriorday=NA
    }
    
    # stores mortality trigger threshold
    if(!is.na(premovsurframe$Mortrigper1000[k])){
      morttrigfrac=premovsurframe$Mortrigper1000[k]/1000
    }else{
      morttrigfrac=NA
    }
    
    # T/F is mortality threshold should be applied consecutively from day given by morttrigpriorday
    if(!is.na(premovsurframe$Ismorttrigconsec[k])){
      ismorttrigconsecutive = premovsurframe$Ismorttrigconsec[k]
    }else{
      ismorttrigconsecutive = NA
    }
    
    # T/F if mortality should be refrigerated (stored) across multiple days for testing
    refrigyesno=premovsurframe$BoolRefrigerationyesno[k]
    
    # maximum number of stored mortality to be tested
    if(!is.na(premovsurframe$Refregerationmaxdead[k])){
      refrimaxdead=premovsurframe$Refregerationmaxdead[k]
    }else{
      refrimaxdead=NA
    }
    
    # number of days for which mortality should be accumulated
    if(!is.na(premovsurframe$Refdaysaccumulate[k])){
      refdaysaccum=premovsurframe$Refdaysaccumulate[k]
    }else{
      refdaysaccum=NA
    }
    
    # T/F for if pens sampled for oral fluids testing should be re-tested on the next oral fluid sampling day
    if(!is.null(premovsurframe$If_repeat_ropepens[k])){
      if(!is.na(premovsurframe$If_repeat_ropepens[k])){
        If_repeat_ropepens=as.logical(premovsurframe$If_repeat_ropepens[k])
      }else{
        If_repeat_ropepens=NA
      }
    }else { If_repeat_ropepens=NA}
    
    # T/F for if the same pens sampled for oral fluids should be tested for blood/tissue sampling
    if(!is.null(premovsurframe$if_premovpen_same_rope[k])){
      #browser()
      if(!is.na(premovsurframe$if_premovpen_same_rope[k])){
        if(!is.na(ropetestdays)){
          if_premovpen_same_rope=as.logical(premovsurframe$if_premovpen_same_rope[k])
        }else{if_premovpen_same_rope=NA}
      }else{
        if_premovpen_same_rope=NA
      }
    } else {if_premovpen_same_rope=NA}
    
    # sampling prioritization for blood/tissue sampling
    bloottargetype=premovsurframe$Bloodtargetsampletype[k]
    # sampling prioritization for oral fluids sampling
    OFtargetype = premovsurframe$Ofsampletype[k]
    
    # number of blood and/or tissue samples to be taken
    indsampsizes<-as.numeric(strsplit(as.character(premovsurframe$Premov.sample.size[k]),split = ",")[[1]])
    
    # assemble all results
    premovsurvlist[[k]]=list(scenarioname=scenarioname,temp_movevec=temp_movevec,
                             ropetestdays=ropetestdays,indsampsizes=indsampsizes,ropepenpercs=ropepenpercs,
                             morttrigpriorday=morttrigpriorday,
                             morttrigfrac=morttrigfrac,
                             ismorttrigconsecutive=ismorttrigconsecutive,
                             refrigyesno=refrigyesno,
                             refrimaxdead=refrimaxdead,
                             bloottargetype=bloottargetype,
                             OFtargetype=OFtargetype,
                             refdaysaccum=refdaysaccum,
                             If_repeat_ropepens=If_repeat_ropepens,
                             if_premovpen_same_rope=if_premovpen_same_rope)
    
  }#end for
  # return list of results
  return(premovsurvlist)
}


# Inputs: x = a number; tol = tolerance for distance x must be from a whole number to be a whole number
# Outpus: T/F for if x is a whole number
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
} 


# Inputs: templis = list of formatted active surveillance protocol inputs
# Outputs: T/F if variables in templist have correct form
check_premovsurvlist<-function(templist){
  
  tflag=TRUE
  # loops through active surveillance protocols in templist
  for  (g in 1:length(templist)){
    # various checks are performed that active surveillance model inputs are formatted correctly
    
    if(!is.na(templist[[g]]$temp_movevec[1])){ 
      if(!is.vector(templist[[g]]$temp_movevec)){
        tflag=FALSE 
      }
      if(length(templist[[g]]$temp_movevec)!=
         length(templist[[g]]$indsampsizes)){
        tflag=FALSE
      }
    }
    if(!is.logical(templist[[g]]$if_premovpen_same_rope)){
      tflag=FALSE 
    } 
    if(!is.logical(templist[[g]]$If_repeat_ropepens)){
      tflag=FALSE 
    }
    
    if(!is.wholenumber(templist[[g]]$bloottargetype)){
      tflag=FALSE 
    }else{
      if(templist[[g]]$bloottargetype>7 |templist[[g]]$bloottargetype<1){
        tflag=FALSE 
      }
      
    }
    
    if(!is.wholenumber(templist[[g]]$OFtargetype)){
      tflag=FALSE 
    }else{
      if(templist[[g]]$OFtargetype>4 |templist[[g]]$OFtargetype<1){
        tflag=FALSE 
      }
      
    }
    
    if(!is.na(templist[[g]]$ropetestdays)){
      if(is.numeric(templist[[g]]$ropetestdays)){
        if(is.na(templist[[g]]$ropepenpercs)){
          tflag=FALSE 
        } else{
          if(length(templist[[g]]$ropepenpercs)==0){
            tflag=FALSE 
          }
        } 
      } 
      
    }
    
    
  }#ending for loop 
  return(tflag)
}