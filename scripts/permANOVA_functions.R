bulkANOVAs <- function(commdat, sampdata){
  # initiate rda summary df
  bulkrdamodelsummary <- data.frame(matrix(nrow = 6, ncol = 7))
  colnames(bulkrdamodelsummary) <- c("time","week", "type", "total_df", "rotation_R2",
                                     "F", "rotation_p_value")
  bulkrdamodelsummary$time <- c(0,1,2,3,4,5)
  bulkrdamodelsummary$week <- c(0,1,2,4,8,12)
  bulkrdamodelsummary$type <- c("bulk", "bulk", "bulk", "bulk", "bulk", "bulk")
  
  # loop through each time point
  for(t in 0:5){
    timecomm <- comm_subset(commdat, sampdata, "bulk", t)
    timesamp <- sample_subset(sampdata, "bulk", t)
    timeclr <- decostand(timecomm+1, method = "clr")
    
    #adonis2 - permanova
    timeanova <- adonis2(timeclr ~ rotation, data = timesamp, method = "euclid")
    
    #add data to summary table
    trow <- t+1
    bulkrdamodelsummary[trow,4] <- timeanova[3,1] # total df
    bulkrdamodelsummary[trow,5] <- timeanova[1,3] # rotation r2
    bulkrdamodelsummary[trow,6] <- timeanova[1,4] # F
    bulkrdamodelsummary[trow,7] <- timeanova[1,5] # rotation p-value
   
  }
  
  return(bulkrdamodelsummary)
}


bagANOVAs <- function(commdat, sampdata){
  #initiate rda summary df
  bagrdamodelsummary <- data.frame(matrix(nrow = 5, ncol = 7))
  colnames(bagrdamodelsummary) <- c("time","week", "type", "total_df" ,"rotation_R2",
                                    "F", "rotation_p_value")
  bagrdamodelsummary$time <- c(1,2,3,4,5)
  bagrdamodelsummary$week <- c(1,2,4,8,12)
  bagrdamodelsummary$type <- c("bag", "bag", "bag", "bag", "bag")
 
  for(t in 1:5){
    timecomm <- comm_subset(commdat, sampdata, "bag", t)
    timesamp <- sample_subset(sampdata, "bag", t)
    timeclr <- decostand(timecomm+1, method = "clr")
    
    #adonis2 - permanova
    timeanova <- adonis2(timeclr ~ rotation, data = timesamp, method = "euclid")
    
    #add data to summary table
    trow <- t
    bagrdamodelsummary[trow,4] <- timeanova[3,1] # total df
    bagrdamodelsummary[trow,5] <- timeanova[1,3] # rotation r2
    bagrdamodelsummary[trow,6] <- timeanova[1,4] # F
    bagrdamodelsummary[trow,7] <- timeanova[1,5] # rotation p-value
  }
  return(bagrdamodelsummary)
} 
