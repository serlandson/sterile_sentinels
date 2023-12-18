setwd("/Users/sonya.erlandson/OneDrive - USDA/sterile sentinels/sr analysis/")

library(data.table)
library(tidyverse)
# use the abundance files without threshold - some of the threshold files are missing.
paths <- dir("emu_sr_bacteria_2020/", pattern = ".tsv$", full.names = TRUE) # get path to files
names(paths) <- basename(paths) #name with filename
alldata <- map_dfr(paths, fread,  .id = "filename") #gather all files into one dataframe
alldata$barcode_id <- str_extract(alldata$filename, "barcode[0-9]{2}") # create a column of barcode ids
samples <- read.csv("../emu statistical analyses/data tables/precursor data tables/sterile_sentinels_sample_data_2020.csv") # load sample data

#summarize data.table relative abundance by tax_id and sample id?
# need two data frames 1) tax_id and taxonomy and 2) tax_id abundance for each barcode
#1
colnames(alldata)
str(alldata)
colnames(alldata)[15] <- "estimated_counts"
taxonomy <- aggregate(estimated_counts ~ ., data = alldata[ ,-c(1,3,11,12,13,14,16)], FUN = sum)

alldata <- alldata[alldata$tax_id != "unassigned", ]
taxonomy <- taxonomy[taxonomy$tax_id != "unassigned", ]

#2
table <- pivot_wider(alldata, id_cols = "species", names_from = "barcode_id", values_from = c("estimated_counts"), values_fill = 0, names_sort = TRUE)
colnames(table)

#get other controls
samples$bacteria_barcode <- gsub("BC", "barcode", samples$Bac.BC.Well)

bag_bulk <- samples$bacteria_barcode[c(which(samples$type == "bag"), which(samples$type == "bulk"))] # get list of bag and bulk smaples
controls <- samples[!samples$bacteria_barcode %in% bag_bulk, ] # get list of control samples
control.comm <- as.data.frame(table[ ,-1]) #get table
control.comm <- control.comm[ ,colnames(control.comm) %in% c(controls$bacteria_barcode)]# subset table to controls
control.comm$species <- table$species # get species names
control.comm$size <- rowSums(table[,-1]) # get size of taxa (overall abundance)
str(control.comm)
control.comm <- control.comm[rowSums(control.comm[,-c(6,7)])>0 , ] # remove taxa with zero abundance in control samples

write.csv(control.comm, "emu_sr_16S_control_communities_2020.csv", quote = FALSE)

# filter counts table
table <- as.data.frame(table) # get into base r dataframe
filtcountstab <- table[!table$species %in% control.comm$species, ] # remove control taxa
row.names(filtcountstab) <- filtcountstab$species # make row names
filtcountstab <- filtcountstab[ ,-1] # remove species column
filtcountstab <- filtcountstab[ ,colnames(filtcountstab) %in% bag_bulk] # remove control samples
filtcountstab <- filtcountstab[rowSums(filtcountstab)>4, ] # remove any empty taxa

# filter by prevalence - keep taxa in 3 or more samples
filt.tab <- as.data.frame(t(filtcountstab)) # transpose count table
str(filt.tab)
filt.tab <- as.data.frame(ceiling(filt.tab))
# get presence absence matrix
bacpa <- filt.tab
bacpa[bacpa>0] <- 1

summary(colSums(bacpa)>=2)
which(colSums(bacpa)<=1)

hist(colSums(bacpa))
filt.tab <- filt.tab[ ,colSums(bacpa)>=3]

# filter taxonomy
taxonomyfilt <- taxonomy[taxonomy$species %in% colnames(filt.tab), ] # filter taxonomy table
summary(taxonomyfilt$species %in% colnames(filt.tab))# check that all match
taxonomyfilt <- taxonomyfilt[match(colnames(filt.tab), taxonomyfilt$species), ]# make sure taxonomy and table are in same order


write.csv(taxonomyfilt, "data tables/emu_16S_taxonomy_2020.csv", quote = FALSE, row.names = FALSE )
write.csv(filt.tab, "data tables/emu_16S_counts_2020.csv", quote = FALSE, row.names = TRUE)


# remove barcode 37 and barcode67
filt.tab <- filt.tab[!row.names(filt.tab) %in% c("barcode37", "barcode67"), ]
summary(colSums(filt.tab)>0)
filt.tab <- filt.tab[, colSums(filt.tab)>0]

write.csv(taxonomyfilt, "data tables/emu_16S_taxonomy_2021.csv", quote = FALSE, row.names = FALSE )
write.csv(filt.tab, "data tables/emu_16S_counts_2021.csv", quote = FALSE, row.names = TRUE)

#####################################################################################
# EMU Output file summary
# get emu output files and put into a summary table
paths <- dir("emu_err_out", pattern = "out$", full.names = TRUE) # get path to files
names(paths) <- basename(paths) #name with filename
alldata <- map_dfr(paths, fread,  .id = "filename") #gather all files into one dataframe
newnames <- str_extract(alldata$filename, "[0-9]+.out")
newnames <-str_remove(newnames, ".out")
alldata$barcode <- newnames

alldatasub <- alldata[ , c(6,2,5)]
colnames(alldatasub) <- c("barcode", "read_status", "count")
finaldata <- pivot_wider(alldatasub, id_cols = "barcode", names_from = "read_status", values_from = "count")
write.table(finaldata, "results/Rdatatables/emu_16S_read_counts.txt", quote = FALSE, row.names = FALSE)

