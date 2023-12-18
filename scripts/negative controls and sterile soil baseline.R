library(tidyverse)
setwd(dir = "/Users/sonya.erlandson/OneDrive - USDA/")

bacmm2020 <- read.csv("sterile sentinels/minimap2 statistical analyses/data tables/16S_control_community_2020.csv")
bacmm2021<- read.csv("sterile sentinels/minimap2 statistical analyses/data tables/16S_control_community_2021.csv")
eukmm2020 <- read.csv("sterile sentinels/minimap2 statistical analyses/data tables/its_control_mock_communities_2020.csv")
eukmm2021 <- read.csv("sterile sentinels/minimap2 statistical analyses/data tables/its_control_mock_communities_2021.csv")

bacemu2020 <- read.csv("sterile sentinels/emu statistical analyses/data tables/emu_16S_control_communities_2020.csv")
bacemu2021 <- read.csv("sterile sentinels/emu statistical analyses/data tables/emu_16S_control_communities_2021.csv")
itsemu2020 <- read.csv("sterile sentinels/emu statistical analyses/data tables/emu_unite_control_communities_2020.csv")
itsemu2021 <- read.csv("sterile sentinels/emu statistical analyses/data tables/emu_unite_control_communities_2021.csv")



#bacmm2020$taxa <- str_split(bacmm2020$taxa, pattern = "; ", n=6, simplify = TRUE)[ ,6]
#bacmm2021$taxa <- str_split(bacmm2021$taxa, pattern = "; ", n=6, simplify = TRUE)[ ,6]

bacmm2020$method <- rep("16S minimap2 silva 2020", length(bacmm2020$sterile.soil))
bacmm2021$method <- rep("16S minimap2 silva 2021", length(bacmm2021$extn.nc))
eukmm2020$method <- rep("ITS minimap2 unite 2020", length(eukmm2020$sterile.soil))
eukmm2021$method <- rep("ITS minimap2 unite 2021", length(eukmm2021$sterile.soil))
bacemu2020$method <- rep("16S emu default 2020", length(bacemu2020$sterile.soil))
bacemu2021$method <- rep("16S emu default 2021", length(bacemu2021$sterile.soil))
itsemu2020$method <- rep("ITS emu unite 2020", length(itsemu2020$sterile.soil))
itsemu2021$method <- rep("ITS emu unite 2021", length(itsemu2021$extn.nc))

bacmm2020 <- bacmm2020[order(bacmm2020$taxa, decreasing = TRUE), ]
bacmm2021 <- bacmm2021[order(bacmm2021$taxa, decreasing = TRUE), ]

eukmm2020 <- eukmm2020[order(eukmm2020$taxa, decreasing = TRUE), ]
eukmm2021 <- eukmm2021[order(eukmm2021$taxa, decreasing = TRUE), ]

bacemu2020 <- bacemu2020[order(bacemu2020$taxa, decreasing = TRUE), ]
bacemu2021 <- bacemu2021[order(bacemu2021$taxa, decreasing = TRUE), ]

itsemu2020 <- itsemu2020[order(itsemu2020$taxa, decreasing = TRUE), ]
itsemu2021 <- itsemu2021[order(itsemu2021$taxa, decreasing = TRUE), ]

#get genus
bacemu2020$genus <- sapply(strsplit(bacemu2020$taxa, split = " "), "[[", 1)
bacemu2021$genus <- sapply(strsplit(bacemu2021$taxa, split = " "), "[[", 1)

itsemu2020$genus <- sapply(strsplit(itsemu2020$taxa, split = "_"), "[[", 1)
itsemu2021$genus <- sapply(strsplit(itsemu2021$taxa, split = "_"), "[[", 1)


cumsum(bacmm2020$sterile.soil)

sterilesoils <- rbind(bacmm2020[ ,c(which(colnames(bacmm2020) == "sterile.soil"), 
              which(colnames(bacmm2020) == "taxa"),
              which(colnames(bacmm2020) == "method"),
              which(colnames(bacmm2020) == "size")),], 
bacemu2020[ ,c(which(colnames(bacemu2020) == "sterile.soil"), 
              which(colnames(bacemu2020) == "taxa"),
              which(colnames(bacemu2020) == "method"),
              which(colnames(bacemu2020) == "size"))],
bacemu2021[ ,c(which(colnames(bacemu2021) == "sterile.soil"), 
              which(colnames(bacemu2021) == "taxa"),
              which(colnames(bacemu2021) == "method"),
              which(colnames(bacemu2021) == "size"))],
itsemu2020[ ,c(which(colnames(itsemu2020) == "sterile.soil"), 
              which(colnames(itsemu2020) == "taxa"),
              which(colnames(itsemu2020) == "method"),
              which(colnames(itsemu2020) == "size"))])




extn <- rbind(bacmm2020[ ,c(which(colnames(bacmm2020) == "extn.nc"), 
                                    which(colnames(bacmm2020) == "taxa"),
                                    which(colnames(bacmm2020) == "method"),
                                    which(colnames(bacmm2020) == "size"))], 
                      bacmm2021[ ,c(which(colnames(bacmm2021) == "extn.nc"), 
                                    which(colnames(bacmm2021) == "taxa"),
                                    which(colnames(bacmm2021) == "method"),
                                    which(colnames(bacmm2021) == "size"))], 
                      eukmm2020[ ,c(which(colnames(eukmm2020) == "extn.nc"), 
                                    which(colnames(eukmm2020) == "taxa"),
                                    which(colnames(eukmm2020) == "method"),
                                    which(colnames(eukmm2020) == "size"))],
                      eukmm2021[ ,c(which(colnames(eukmm2021) == "extn.nc"), 
                                    which(colnames(eukmm2021) == "taxa"),
                                    which(colnames(eukmm2021) == "method"),
                                    which(colnames(eukmm2021) == "size"))],
                      bacemu2021[ ,c(which(colnames(bacemu2021) == "extn.nc"), 
                                     which(colnames(bacemu2021) == "taxa"),
                                     which(colnames(bacemu2021) == "method"),
                                     which(colnames(bacemu2021) == "size"))],
                      itsemu2020[ ,c(which(colnames(itsemu2020) == "extn.nc"), 
                                     which(colnames(itsemu2020) == "taxa"),
                                     which(colnames(itsemu2020) == "method"),
                                     which(colnames(itsemu2020) == "size"))],
                      itsemu2021[ ,c(which(colnames(itsemu2021) == "extn.nc"), 
                                     which(colnames(itsemu2021) == "taxa"),
                                     which(colnames(itsemu2021) == "method"),
                                     which(colnames(itsemu2021) == "size"))])

#itsemu2020- the pcr nc is actually pcr2 not pcr1

pcr.nc <- rbind(bacmm2020[ ,c(which(colnames(bacmm2020) == "pcr1.nc"), 
                                    which(colnames(bacmm2020) == "taxa"),
                                    which(colnames(bacmm2020) == "method"),
                                    which(colnames(bacmm2020) == "size"))], 
                      bacmm2021[ ,c(which(colnames(bacmm2021) == "pcr1.nc"), 
                                    which(colnames(bacmm2021) == "taxa"),
                                    which(colnames(bacmm2021) == "method"),
                                    which(colnames(bacmm2021) == "size"))], 
                      bacemu2021[ ,c(which(colnames(bacemu2021) == "pcr1.nc"), 
                                     which(colnames(bacemu2021) == "taxa"),
                                     which(colnames(bacemu2021) == "method"),
                                     which(colnames(bacemu2021) == "size"))],
                      itsemu2020[ ,c(which(colnames(itsemu2020) == "pcr1.nc"), 
                                     which(colnames(itsemu2020) == "taxa"),
                                     which(colnames(itsemu2020) == "method"),
                                     which(colnames(itsemu2020) == "size"))])




sterilesoils <- sterilesoils[sterilesoils$sterile.soil>0,]
extn <- extn[extn$extn.nc>0,]
pcr.nc <- pcr.nc[pcr.nc$pcr1.nc>0,]


# make labels centered in cumulative count bar
sterilesoils <- sterilesoils %>%
  group_by(method) %>%
  mutate(label_y = cumsum(sterile.soil)- 0.5*sterile.soil)

extn <- extn %>%
  group_by(method) %>%
  mutate(label_y = cumsum(extn.nc)-0.5*extn.nc)

pcr.nc <- pcr.nc %>%
  group_by(method) %>%
  mutate(label_y = cumsum(pcr1.nc)-0.5*pcr1.nc)


ggplot(sterilesoils)+
    geom_bar(aes(y=sterile.soil, x = method, fill = taxa), stat = "identity", position = "stack", show.legend = FALSE)+
  #scale_fill_manual(values = newcols)+
  geom_text(aes(y=label_y, x = method, label = taxa), stat = "identity", colour = "black", size = 2.5, check_overlap = T)+
  xlab("amplicon + bioinformatics method + database + year")+
  ylab("taxa read count in control")+theme_minimal()+ggtitle("Sterile Soil")

ggplot(extn)+
  geom_bar(aes(y=extn.nc, x = method, fill = taxa), stat = "identity", position = "stack", show.legend = FALSE)+
  #scale_fill_manual(values = newcols)+
  geom_text(aes(y=label_y, x = method, label = taxa), stat = "identity", colour = "black", size = 2.5)+
  xlab("amplicon + bioinformatics method + database + year")+
  ylab("taxa read count in control")+theme_minimal()+ggtitle("Extraction Blank")

ggplot(pcr.nc)+
  geom_bar(aes(y=pcr1.nc, x = method, fill = taxa), stat = "identity", position = "stack", show.legend = FALSE)+
  #scale_fill_manual(values = newcols)+
  geom_text(aes(y=label_y, x = method, label = taxa), stat = "identity", colour = "black", size = 2.5)+
  xlab("amplicon + bioinformatics method + database + year")+
  ylab("taxa read count in control")+theme_minimal()+ggtitle("PCR negative control")


# what we are actually concerned about is: Are taxa present in the controls also present in the samples? 

#sterile soil
sslarge <- sterilesoils[sterilesoils$sizesterilesoils$size, ]
ssnotindata <- sterilesoils[sterilesoils$sterile.soil<=sterilesoils$size, ]



