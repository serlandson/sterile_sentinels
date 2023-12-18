#### subset code
# community table in rows=samples format
# sample data and community samples must be in the same order.
sample_subset <- function(sample_data, samp_type=NULL, which_time=NULL){
  if(!is.null(samp_type) && is.null(which_time)){
    sampsub <- sample_data[sample_data[ ,"type"] == samp_type, ]
    return(sampsub)
  }
  if(is.null(samp_type) && !is.null(which_time)){
    sampsub <- sample_data[sample_data[ ,"time"] == which_time, ]
    return(sampsub)
  }
  if(!is.null(samp_type) %% !is.null(which_time)){
    sampsub <- sample_data[sample_data[ ,"type"] == samp_type & sample_data[ ,"time"] == which_time, ]
    return(sampsub)
  }
  
}

comm_subset <- function(x, sample_data, samp_type=NULL, which_time=NULL){
  if(!is.null(samp_type) && is.null(which_time)){
    testcomm <- x[sample_data[ ,"type"] == samp_type, ]
    testcomm <- testcomm[ ,colSums(testcomm)>3]
    return(testcomm)
  }
  if(is.null(samp_type) && !is.null(which_time)){
    testcomm <- x[sample_data[ ,"time"] == which_time, ]
    testcomm <- testcomm[ ,colSums(testcomm)>3]
    return(testcomm)
  }
  if(!is.null(samp_type) %% !is.null(which_time)){
    testcomm <- x[sample_data[ ,"type"] == samp_type & sample_data[ ,"time"] == which_time, ]
    testcomm <- testcomm[ ,colSums(testcomm)>3]
    return(testcomm)
  }
  
}
