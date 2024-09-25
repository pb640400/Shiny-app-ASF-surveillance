# This file contains a function for deriving objects necessary to run the transmission model with sim_parms$trans_type == 7.
# The objects include a matrix and array indicating adjacent pens, an array indicating edge pens, and an array indicating
# rooms that can transmit infection to each other

# Inputs: sim_parms = simulation parameters
# Outputs: list containing the objects pen_distmult_arr = matrix indicating the pen pairs that are adjacent; is_edge_pen_arr = array indicating
#  the pens that are at the end of a row; pen_distmult_linvec = linearized version of pen_distmult_arr; room_distmult_linvec =
#  array indicating room transmission partners; room_of_pen_arr = array with the room of each pen
get_barn_pen_arraylist<-function(sim_parms=sim_parms){
   
   # initialize arrays indicating the room of each pen and whether or not a pen is at the end of a row ("edge pen")
   room_of_pen_arr<-rep(NA,sim_parms$npens)
   is_edge_pen_arr<-rep(0,sim_parms$npens)

   
   # assigns each pen a room
   for (i in 1:sim_parms$nrooms){
      for (k in 1:(sim_parms$npens/sim_parms$nrooms)){
         room_of_pen_arr[k+((i-1)*(sim_parms$npens/sim_parms$nrooms))]<-i
      }
   }

   
   # initializes array indicating adjacent pens
   pen_distmult_arr<-matrix(1,nrow =sim_parms$npen,ncol = sim_parms$npen)
   # initializes array indicating rooms that can transmit infection to each other
   room_distmult_arr<-matrix(1,nrow =sim_parms$nrooms,ncol = sim_parms$nrooms)


   # Function that returns 1) a matrix indicating adjacent pairs of pens ("pen_distmult_arr") where if 
   # pen_distmult_arr[i,j] == 1 pen i and j are adjacent while if pen_distmult_arr[i,j] == 0 they are not and
   # 2) an array indicating which pen are at the end of a row ("is_edge_pen_arr")
   get_adjacentpen_dmat <- function(){
      
      # barn with a single row per room
      if(sim_parms$prem_config==1 ){
         # initializes matrix for indicating adjacent pens
         pen_distmult_arr<-matrix(0,nrow =sim_parms$npen,ncol = sim_parms$npen) 
         # loop through the number of rooms
         for (i in 1:sim_parms$nrooms){
            # picks pen indices for the pens in room i
            temp_pensinroomarr<-((i-1)*(sim_parms$npens/sim_parms$nrooms)+1):((i)*(sim_parms$npens/sim_parms$nrooms))
            # first and last pens in the row are edge pens
            is_edge_pen_arr[temp_pensinroomarr[1]]<-1
            is_edge_pen_arr[temp_pensinroomarr[length(temp_pensinroomarr)]]<-1
            # all pens in the room are assigned the same row
            temp_rowofpenarr<-rep(NA,length(temp_pensinroomarr))
            temp_rowofpenarr[]<-1
            # loops through pairs of pens in the room. Pens with a difference of 1 between the indices are adjacent
            for (g in  ((i-1)*(sim_parms$npens/sim_parms$nrooms)+1):((i)*(sim_parms$npens/sim_parms$nrooms))){
               for (h in  ((i-1)*(sim_parms$npens/sim_parms$nrooms)+1):((i)*(sim_parms$npens/sim_parms$nrooms))){
                  if (g != h){
                     if(temp_rowofpenarr[which(temp_pensinroomarr==g)] == temp_rowofpenarr[which(temp_pensinroomarr==h)]){
                        if(abs(g-h)==1){
                           pen_distmult_arr[g,h] <- 1
                        }
                        
                     }  
                     
                  }
                  
               }#end h loop
               
            }# end g loop
            
         }# end i loop through rooms
         
      }
      
      
      ## barn where each room contains 2 rows separated by an aisle
      if(sim_parms$prem_config==2 ){
         # initializes matrix for indicating adjacent pens for each pen
         pen_distmult_arr<-matrix(0,nrow =sim_parms$npen,ncol = sim_parms$npen) 
         # loops through rooms in barn
         for (i in 1:sim_parms$nrooms){
            # selects pen indices of the pens in room i
            temp_pensinroomarr<-((i-1)*(sim_parms$npens/sim_parms$nrooms)+1):((i)*(sim_parms$npens/sim_parms$nrooms))
            
            # initializes array for the row of each pen in the room
            temp_rowofpenarr<-rep(NA,length(temp_pensinroomarr))
            # divides pens in rooms evenly into two rows where the row of the matrix defines the pens in the same row
            temp_peninroominmatrix<-(matrix(temp_pensinroomarr,nrow = 2,byrow = TRUE))
            
            # loops through rows in room. First and last pen in each row are edge pens.
            for (zzz in 1:nrow(temp_peninroominmatrix)){
               is_edge_pen_arr[temp_peninroominmatrix[zzz,1]]<-1 
               is_edge_pen_arr[temp_peninroominmatrix[zzz,ncol(temp_peninroominmatrix)]]<-1 
            }
            
            # assigns each pen a row in temp_rowofpenarr
            for (z in 1:2){
               for (y in 1:(length(temp_pensinroomarr) / 2)){
                  temp_rowofpenarr[y+((z-1)*(length(temp_pensinroomarr) / 2))]<-z
               }
            }
            
            
            # loops through pen pairs in the same room. Pens in the same row that are one index apart are adjacent.
            for (g in  ((i-1)*(sim_parms$npens/sim_parms$nrooms)+1):((i)*(sim_parms$npens/sim_parms$nrooms))){
               for (h in  ((i-1)*(sim_parms$npens/sim_parms$nrooms)+1):((i)*(sim_parms$npens/sim_parms$nrooms))){
                  if (g != h){
                     
                     if(temp_rowofpenarr[which(temp_pensinroomarr==g)] == temp_rowofpenarr[which(temp_pensinroomarr==h)]){
                        if (abs(g-h)==1){
                           pen_distmult_arr[g,h] <- 1
                        }
                     }  
                     
                  }
                  
               }#end h loop
               
            }# end g loop
            
         }# end i loop
         
      }
      
      
      
      return(list(pen_distmult_arr,is_edge_pen_arr))
   }#end distancematrix function
   

   # stores matrix indicating adjacent pen pairs
   pen_distmult_arr<-get_adjacentpen_dmat()[[1]]
   # stores array indicating edge pens
   is_edge_pen_arr<-get_adjacentpen_dmat()[[2]]
   # linearizes matrix indicating adjacent pen pairs
   pen_distmult_linvec<-Linearize_row_major(pen_distmult_arr)
   # linearizes matrix indicating rooms that can transmit infection to one another
   room_distmult_linvec<-Linearize_row_major(room_distmult_arr)
   
   outlist=list(pen_distmult_arr=pen_distmult_arr,
                is_edge_pen_arr=is_edge_pen_arr,pen_distmult_linvec=pen_distmult_linvec,
                room_distmult_linvec=room_distmult_linvec,room_of_pen_arr=room_of_pen_arr)
   return(outlist)
   
}