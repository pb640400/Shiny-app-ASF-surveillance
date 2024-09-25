# This file contains helper functions for linearizing an array, making a matrix from an array input, simulating from
# a beta-PERT distribution, and simulating from a hypergeometric distribution


# The function linerizes a 2 d array in a row major fashion. This is useful
# to reconstruct the array in C by point to the head of each row
# so that we can use same indices order to refer to the array in C
Linearize_row_major<-function(inp_arr){
  return(do.call("c",lapply(1:nrow(inp_arr),function(g){inp_arr[g,]})))
}


# create a matrix from a 2-D array
get_r2darr_fromrowmaj<-function(in_arr,in_rows,in_cols){
  return(t(array(in_arr,dim=c(in_cols,in_rows))))
}


# function for simulating from a beta-PERT distribution
rspert <- function( n, x.min, x.max, x.mode, lambda = 4 ){
  #browser()
  if( x.min > x.max || x.mode > x.max || x.mode < x.min ) stop( "invalid parameters" );
  
  x.range <- x.max - x.min;
  if( x.range == 0 ) return( rep( x.min, n ));
  
  mu <- ( x.min + x.max + lambda * x.mode ) / ( lambda + 2 );
  
  # special case if mu == mode
  if( mu == x.mode ){
    v <- ( lambda / 2 ) + 1
  }
  else {
    v <- (( mu - x.min ) * ( 2 * x.mode - x.min - x.max )) /
      (( x.mode - mu ) * ( x.max - x.min ));
  }
  
  w <- ( v * ( x.max - mu )) / ( mu - x.min );
  return ( rbeta( n, v, w ) * x.range + x.min );
}


# function to catch exception when total mortality is 0 when simulating from a hypergeometric distribution
wraphyper<-function(normal,  disease,  nswabs ){
  if((normal+disease)==0){
    hypergeometric=0
    
  }else{
    if (normal + disease >= nswabs ){
      hypergeometric = rhyper(1,disease,normal,nswabs )
    } else {
      hypergeometric = rhyper(1,disease ,normal ,(normal + disease))
    }
  }
  return(  hypergeometric)
  
}