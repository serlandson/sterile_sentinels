
source("../OneDrive - USDA/files to save/sterile sentinels/manuscript/scripts/subset_samples_community_table.R")

summarydf <- data.frame()

for(t in 1:5){
  timecommbulk <- comm_subset(comm, samples, samp_type = "bulk", which_time = t)
  timesampbulk <- sample_subset(samples, samp_type = "bulk", which_time = t)
  timetaxbulk <- taxtrait[taxtrait$species %in% colnames(timecommbulk), ]
  timetaxbulk <- timetaxbulk[!is.na(timetaxbulk$genome_size), ]
  timecommbulk <- timecommbulk[, colnames(timecommbulk) %in% timetaxbulk$species]
  
  
  timecommbag <- comm_subset(comm, samples, samp_type = "bag", which_time = t)
  timesampbag <- sample_subset(samples, samp_type = "bag", which_time = t)
  timetaxbag <- taxtrait[taxtrait$species %in% colnames(timecommbag), ]
  timetaxbag <- timetaxbag[!is.na(timetaxbag$genome_size), ]
  timecommbag <- timecommbag[, colnames(timecommbag) %in% timetaxbag$species]
  
  timebagrelabund <- colSums(timecommbag)/sum(colSums(timecommbag))
  timebulkrelabund <- colSums(timecommbulk)/sum(colSums(timecommbulk))
  
  newdf <- data.frame(timetaxbag$species, timetaxbag$genome_size, timebagrelabund, 
                      rep("sentinels", length(timetaxbag$species)), rep(t, length(timetaxbag$species)))
  colnames(newdf) <- c("species", "genome_size","relative_abundance", "type", "time")
  
  newdf2 <- data.frame(timetaxbulk$species, timetaxbulk$genome_size, timebulkrelabund, 
                       rep("bulk", length(timetaxbulk$species)),rep(t, length(timetaxbulk$species)))
  colnames(newdf2) <- c("species", "genome_size", "relative_abundance","type", "time")
  
  subdf <- rbind(newdf, newdf2)
  summarydf <- rbind(summarydf, subdf)
}

write.csv(summarydf, "genome size ttest bacteria 2020.csv", quote =F)
write.csv(summarydf, "genome size ttest bacteria 2021.csv", quote =F)

##### plot
plotdf <- read.csv("filtered differential abundance data tables/trait differential abundance tables/genome size ttest bacteria 2020.csv")

plotdf <- summarydf
plotdf2 <- summarydf
hist(plotdf2$genome_size)

plotdf <- genomesize2021
plotdf2$time <- factor(plotdf$time, labels = c(1,2,3,4,5))

ttestslist <- list()

for(t in 1:5){
  ttestslist[[t]] <- wilcox.test(plotdf[plotdf$time==t, ]$genome_size ~ plotdf[plotdf$time==t, ]$type, exact = F)
}
ttestslist

weighted.mean(plotdf$genome_size[plotdf$type=="sentinels"], plotdf$relative_abundance[plotdf$type=="sentinels"])
weighted.mean(plotdf$genome_size[plotdf$type=="bulk"], plotdf$relative_abundance[plotdf$type=="bulk"])

mean(plotdf$genome_size[plotdf$type=="sentinels"])
mean(plotdf$genome_size[plotdf$type=="bulk"])

library(ggplot2)

plotdf$time <- factor(plotdf$time, labels = c(1,2,4,8,12))
genome2020 <- ggplot(plotdf2)+
  geom_violin(aes(x= type,
                  y = genome_size))+
  #geom_text(aes(x = "sentinels", y = 150, label = "wilcoxon t.test pvalue = 0.0009", fontface="plain"), size = 3)+
  #geom_text(aes(x = "sentinels", y = -1,label = "weighted mean doubling time = 1.9 h"), size = 3)+
  #geom_text(aes(x = "bulk", y = -1, label = "weighted mean doubling time = 9.9 h"),size = 3)+
  #theme(axis.title.y = element_blank())+
  ylab("potential genome size (bp)")+
  xlab("sample type")+
  facet_wrap(~time,ncol = 5)+
  ggtitle("Bacteria 2020")

library(patchwork)
genomesize <- genome2020+genome2021+plot_layout(guides = "collect")
genomesize
ggsave("genome size.tiff", plot = genomesize, 
       device = "tiff", units = "in", height = 4, width = 10.5)

genomesize2020 <- plotdf2
genomesize2021 <- plotdf

sample_rotation_subset <- function(sample_data, samp_type=NULL, which_time=NULL){
  if(!is.null(samp_type) && is.null(which_time)){
    sampsub <- sample_data[sample_data[ ,"rotation"] == samp_type, ]
    return(sampsub)
  }
  if(is.null(samp_type) && !is.null(which_time)){
    sampsub <- sample_data[sample_data[ ,"time"] == which_time, ]
    return(sampsub)
  }
  if(!is.null(samp_type) %% !is.null(which_time)){
    sampsub <- sample_data[sample_data[ ,"rotation"] == samp_type & sample_data[ ,"time"] == which_time, ]
    return(sampsub)
  }
  
}

comm_rotation_subset <- function(x, sample_data, samp_type=NULL, which_time=NULL){
  if(!is.null(samp_type) && is.null(which_time)){
    testcomm <- x[sample_data[ ,"rotation"] == samp_type, ]
    testcomm <- testcomm[ ,colSums(testcomm)>0]
    return(testcomm)
  }
  if(is.null(samp_type) && !is.null(which_time)){
    testcomm <- x[sample_data[ ,"time"] == which_time, ]
    testcomm <- testcomm[ ,colSums(testcomm)>0]
    return(testcomm)
  }
  if(!is.null(samp_type) %% !is.null(which_time)){
    testcomm <- x[sample_data[ ,"rotation"] == samp_type & sample_data[ ,"time"] == which_time, ]
    testcomm <- testcomm[ ,colSums(testcomm)>0]
    return(testcomm)
  }
  
}


summarydf <- data.frame()

bulk <- comm_subset(comm, samples, samp_type = "bulk", which_time = NULL)
bulksamp <- sample_subset(samples, samp_type = "bulk")
bag <-  comm_subset(comm, samples, samp_type = "bag", which_time = NULL)
bagsamp <-  sample_subset(samples, samp_type = "bag")

for(t in 1:5){
  timecommCS <- comm_rotation_subset(bulk, bulksamp, samp_type = "CS", which_time = t)
  timesampCS <- sample_rotation_subset(bulksamp, samp_type = "CS", which_time = t)
  timetaxCS <- taxtrait[taxtrait$species %in% colnames(timecommCS), ]
  timetaxCS <- timetaxCS[!is.na(timetaxCS$genome_size), ]
  timecommCS <- timecommCS[, colnames(timecommCS) %in% timetaxCS$species]
  
  timecommCSSwP <- comm_rotation_subset(bulk, bulksamp, samp_type = "CSSwP", which_time = t)
  timesampCSSwP <- sample_rotation_subset(bulksamp, samp_type = "CSSwP", which_time = t)
  timetaxCSSwP <- taxtrait[taxtrait$species %in% colnames(timecommCSSwP), ]
  timetaxCSSwP <- timetaxCSSwP[!is.na(timetaxCSSwP$genome_size), ]
  timecommCSSwP <- timecommCSSwP[, colnames(timecommCSSwP) %in% timetaxCSSwP$species]
  
  timeCSSwPrelabund <- colSums(timecommCSSwP)/sum(colSums(timecommCSSwP))
  timeCSrelabund <- colSums(timecommCS)/sum(colSums(timecommCS))
  
  
  newdf <- data.frame(timetaxCSSwP$species, timetaxCSSwP$genome_size, timeCSSwPrelabund, 
                      rep("CSSwP", length(timetaxCSSwP$species)), rep(t, length(timetaxCSSwP$species)))
  colnames(newdf) <- c("species", "genome_size","relative_abundance", "rotation", "time")
  
  newdf2 <- data.frame(timetaxCS$species, timetaxCS$genome_size, timeCSrelabund, 
                       rep("CS", length(timetaxCS$species)),rep(t, length(timetaxCS$species)))
  colnames(newdf2) <- c("species", "genome_size", "relative_abundance","rotation", "time")
  
  subdf <- rbind(newdf, newdf2)
  summarydf <- rbind(summarydf, subdf)
}


write.csv(summarydf, "genome size ttest ROTATION BULK bacteria 2020.csv", quote =F)

write.csv(summarydf, "genome size ttest ROTATION BULK bacteria 2021.csv", quote =F)
write.csv(summarydf, "genome size ttest ROTATION SENTINEL bacteria 2021.csv", quote =F)

##### plot
plotdf <- summarydf
hist(plotdf$genome_size)

rotationttestslist <- list()

for(t in 1:5){
  rotationttestslist[[t]] <- wilcox.test(plotdf[plotdf$time==t, ]$genome_size ~ plotdf[plotdf$time==t, ]$rotation, exact = F)
}


rotationttestslist


weighted.mean(plotdf$genome_size[plotdf$rotation=="CS"], plotdf$relative_abundance[plotdf$rotation=="CS"])
weighted.mean(plotdf$genome_size[plotdf$rotation=="CSSwP"], plotdf$relative_abundance[plotdf$rotation=="CSSwP"])

mean(plotdf$genome_size[plotdf$rotation=="CS"])
mean(plotdf$genome_size[plotdf$rotation=="CSSwP"])

library(ggplot2)

plotdf$week <- factor(plotdf$time, labels = c(1,2,4,8,12))
genomebulk <- ggplot(plotdf)+
  geom_violin(aes(x= rotation,
                  y = genome_size))+
  #geom_text(aes(x = "sentinels", y = 150, label = "wilcoxon t.test pvalue = 0.0009", fontface="plain"), size = 3)+
  #geom_text(aes(x = "sentinels", y = -1,label = "weighted mean doubling time = 1.9 h"), size = 3)+
  #geom_text(aes(x = "bulk", y = -1, label = "weighted mean doubling time = 9.9 h"),size = 3)+
  ylab("bacteria genome size")+
  xlab("sample type")+
  facet_wrap(~week,ncol = 5)+
  ggtitle("Bacteria Bulk 2020")



genomebulk2021 <- genomebulk
genomebag2021 <- genomebag

genomebag+genomebulk

