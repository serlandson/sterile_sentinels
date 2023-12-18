setwd("/Users/sonya.erlandson/OneDrive - USDA/files to save/sterile sentinels/sr analysis/")
library(tidyverse)
library(RColorBrewer)
library(ggtext)

#species
allrot <- read.csv("differential abundance data tables/all_taxa_Rotation_diff_abund_no_filter.csv")
allrot <- separate(allrot, amplicon, into = c("sample", "amplicon", "year"), sep = "_")
allrot <- separate(allrot, type, into = c("type", "week"), sep = "_")

fungi <- allrot[allrot$amplicon=="fungi", ]
bacteria <- allrot[allrot$amplicon=="bac", ]

length(unique(fungi$species))
length(unique(bacteria$species))

# chose taxa found in multiple samples
fungibigeffect <- fungi[abs(fungi$effect)>1, ] # effect size larger than 1
#how many taxa?
# table of differentially abundant taxa 
# type | week | rotation/effect | year 
csfungi <- fungibigeffect[fungibigeffect$effect<0, ]
csswpfungi <- fungibigeffect[fungibigeffect$effect>0, ]


# count the number of taxa differentially abundant in each week/type/year.
table(csfungi[csfungi$type=="bulk" & csfungi$year == 2020, ]$week)
table(csfungi[csfungi$type=="sentinel" & csfungi$year == 2020, ]$week)

table(csswpfungi[csswpfungi$type=="bulk" & csswpfungi$year == 2020, ]$week)
table(csswpfungi[csswpfungi$type=="sentinel" & csswpfungi$year == 2020, ]$week)

table(csfungi[csfungi$type=="bulk" & csfungi$year == 2021, ]$week)
table(csfungi[csfungi$type=="sentinel" & csfungi$year == 2021, ]$week)

table(csswpfungi[csswpfungi$type=="bulk" & csswpfungi$year == 2021, ]$week)
table(csswpfungi[csswpfungi$type=="sentinel" & csswpfungi$year == 2021, ]$week)


length(unique(csfungi$species)) # 103
length(unique(csfungi$species[csfungi$type=="bulk"])) #73
length(unique(csfungi$species[csfungi$type=="bulk" & csfungi$year == 2020])) # 56
length(unique(csfungi$species[csfungi$type=="bulk" & csfungi$year == 2021])) # 23
length(unique(csfungi$species[csfungi$type=="sentinel"])) #49
length(unique(csfungi$species[csfungi$type=="sentinel" & csfungi$year == 2020])) # 26
length(unique(csfungi$species[csfungi$type=="sentinel" & csfungi$year == 2021])) # 25

length(unique(csswpfungi$species)) # 99
length(unique(csswpfungi$species[csswpfungi$type=="bulk"]))
length(unique(csswpfungi$species[csswpfungi$type=="bulk" & csswpfungi$year == 2020]))#47
length(unique(csswpfungi$species[csswpfungi$type=="bulk" & csswpfungi$year == 2021]))#36
length(unique(csswpfungi$species[csswpfungi$type=="sentinel"])) 
length(unique(csswpfungi$species[csswpfungi$type=="sentinel" & csswpfungi$year == 2020]))#29
length(unique(csswpfungi$species[csswpfungi$type=="sentinel" & csswpfungi$year == 2021]))#22


dups <- as.data.frame(table(fungibigeffect$species))
dups <- dups[dups$Freq>4,]
fungibigeffect <- fungibigeffect[fungibigeffect$species %in% dups$Var1, ]

selectedspp <- unique(fungibigeffect$species) # 17 taxa
#subset to just these 17 taxa with large effect sizes
fungisubset <- fungi[fungi$species %in% selectedspp, ]


# plot colors
library(RColorBrewer)
cols <- brewer.pal(9, name = "Blues")
cols <- cols[4:9]

fungisubset <- fungisubset[order(fungisubset$effect), ]
sppord <- unique(fungisubset$species)
fungisubset$species <- factor(fungisubset$species, levels = sppord)
fungisubset$year <- as.factor(fungisubset$year)
fungisubset$week <- factor(fungisubset$week, labels = c(0,1,2,4,8,12))

fungilabels <- sapply(strsplit(sppord, split = "_"), "[[", 1)
fungilabels <- paste0("*",fungilabels, "* ", "sp.")

# plot
fungidiffabund <- ggplot(fungisubset)+
  geom_point(aes(y = species, 
                 x = effect, 
                 color = week,
                 size = rab.all, 
                 shape = year,
                 alpha = abs(effect)>1), size = 4)+
  theme_light()+
  scale_color_manual(values = cols, 
                     name = c("sampling week"),
                     labels = c("0","1","2","4","8","12"))+
  scale_shape_manual(values = c(17,19,2,1),
                     name = c("year"))+
  geom_vline(xintercept=-1,linetype = 2, alpha = 0.5 )+
  geom_vline(xintercept= 1,linetype = 2, alpha = 0.5 )+
  labs(x = "effect size (CS vs CSSwP)",
       y = "differentially abundant taxa")+
  guides(alpha = "none")+
  scale_y_discrete(labels = fungilabels)+
  facet_grid(~type)+
  theme(axis.text.y = element_markdown())
fungidiffabund
ggsave("fungi_rotation_diff_abund.tiff", plot = fungidiffabund, width = 8, height = 6, units = "in")

##### bacteria

# chose taxa found in multiple samples
bacbigeffect <- bacteria[abs(bacteria$effect)>1, ]

csbac <- bacbigeffect[bacbigeffect$effect<0, ]
csswpbac <- bacbigeffect[bacbigeffect$effect>0, ]

table(csbac[csbac$type=="bulk" & csbac$year == 2020, ]$week)
table(csbac[csbac$type=="sentinel" & csbac$year == 2020, ]$week)

table(csswpbac[csswpbac$type=="bulk" & csswpbac$year == 2020, ]$week)
table(csswpbac[csswpbac$type=="sentinel" & csswpbac$year == 2020, ]$week)

table(csbac[csbac$type=="bulk" & csbac$year == 2021, ]$week)
table(csbac[csbac$type=="sentinel" & csbac$year == 2021, ]$week)

table(csswpbac[csswpbac$type=="bulk" & csswpbac$year == 2021, ]$week)
table(csswpbac[csswpbac$type=="sentinel" & csswpbac$year == 2021, ]$week)

length(unique(csbac$species)) # 172
length(unique(csbac$species[csbac$type=="bulk"])) #
length(unique(csbac$species[csbac$type=="bulk" & csbac$year == 2020])) # 33
length(unique(csbac$species[csbac$type=="bulk" & csbac$year == 2021])) # 56
length(unique(csbac$species[csbac$type=="sentinel"])) #
length(unique(csbac$species[csbac$type=="sentinel" & csbac$year == 2020])) # 37
length(unique(csbac$species[csbac$type=="sentinel" & csbac$year == 2021])) # 56

length(unique(csswpbac$species)) # 175
length(unique(csswpbac$species[csswpbac$type=="bulk"]))
length(unique(csswpbac$species[csswpbac$type=="bulk" & csswpbac$year == 2020]))# 47
length(unique(csswpbac$species[csswpbac$type=="bulk" & csswpbac$year == 2021]))# 20
length(unique(csswpbac$species[csswpbac$type=="sentinel"])) #118
length(unique(csswpbac$species[csswpbac$type=="sentinel" & csswpbac$year == 2020]))# 78
length(unique(csswpbac$species[csswpbac$type=="sentinel" & csswpbac$year == 2021]))# 47

# narrow down by selecting taxa that are more abundant or are significant in more than one week
dups <- as.data.frame(table(bacbigeffect$species))
dups <- dups[dups$Freq>2,]
bacbigeffect <- bacbigeffect[bacbigeffect$species %in% dups$Var1, ]

selectedspp <- unique(bacbigeffect$species) # 21 taxa
#subset to just these 17 taxa with large effect sizes
bacsubset <- bacteria[bacteria$species %in% selectedspp, ]

# plot colors
cols <- brewer.pal(9, name = "Blues")
cols <- cols[4:9]

bacsubset <- bacsubset[order(bacsubset$effect), ]
sppord <- unique(bacsubset$species)
bacsubset$species <- factor(bacsubset$species, levels = sppord)
bacsubset$year <- as.factor(bacsubset$year)
bacsubset$week <- factor(bacsubset$week, labels = c(0,1,2,4,8,12))

baclabels <- sapply(strsplit(sppord, split = " "), "[[", 1)
baclabels <- paste0("*", baclabels, "* ", "sp.")

# plot
bacdiffabund <- ggplot(bacsubset)+
  geom_point(aes(y = species, 
                 x = effect, 
                 color = week,
                 size = rab.all, 
                 shape = year,
                 alpha = abs(effect)>1), size = 4)+
  theme_light()+
  scale_color_manual(values = cols, 
                     name = c("sampling week"),
                     labels = c("0","1","2","4","8","12"))+
  scale_shape_manual(values = c(17,19,2,1),
                     name = c("year"))+
  geom_vline(xintercept=-1,linetype = 2, alpha = 0.5 )+
  geom_vline(xintercept= 1,linetype = 2, alpha = 0.5 )+
  labs(x = "effect size (CS vs CSSwP)",
       y = "differentially abundant taxa")+
  guides(alpha = "none")+
  scale_y_discrete(labels = baclabels)+
  facet_grid(~type)+
  theme(axis.text.y = element_markdown())

bacdiffabund

ggsave("bacteria_rotation_diff_abund.tiff", plot = bacdiffabund, width = 8, height = 6, units = "in")


# compare number of differentially abundant taxa between 
# sentinels and bulk soil.
counts <- read.csv("../OneDrive - USDA/files to save/sterile sentinels/manuscript/data/counts of differentially abundant taxa between rotation.csv")

fungicounts <- counts[grep("fungi", counts$amplicon_year), ]
fungicounts <- fungicounts[fungicounts$week!=0,]
baccounts <- counts[grep("bacteria", counts$amplicon_year), ]
baccounts <- baccounts[baccounts$week!=0,]


# overall diff abund counts
fungicounts$all <- fungicounts$cs+fungicounts$csswp
baccounts$all <- baccounts$cs+baccounts$csswp

#fungi
hist(fungicounts$all[fungicounts$type=="bulk"]-fungicounts$all[fungicounts$type=="sentinel"])
t.test(fungicounts$all[fungicounts$type=="bulk"], 
       fungicounts$all[fungicounts$type=="sentinel"], paired = TRUE)

#bacteria
hist(baccounts$all[baccounts$type=="bulk"]-baccounts$all[baccounts$type=="sentinel"])
t.test(baccounts$all[baccounts$type=="bulk"], baccounts$all[baccounts$type=="sentinel"], paired = TRUE)


hist(fungicounts$cs[fungicounts$type == "bulk"]-fungicounts$cs[fungicounts$type=="sentinel"])
hist(fungicounts$csswp[fungicounts$type == "bulk"]-fungicounts$csswp[fungicounts$type=="sentinel"])
#CS - no difference
ttestfungics <- t.test(fungicounts$cs[fungicounts$type=="bulk"], 
       y = fungicounts$cs[fungicounts$type=="sentinel"],  
       paired = TRUE)
knitr::kable(c(ttestfungics$statistic, ttestfungics$p.value, 
               ttestfungics$estimate, ttestfungics$conf.int), format = "simple")


#CSSwP - different at p=0.01
#
ttestfungicsswp <- t.test(fungicounts$csswp[fungicounts$type=="bulk"], 
       y = fungicounts$csswp[fungicounts$type=="sentinel"],  
       paired = TRUE)
mean(fungicounts$csswp[fungicounts$type=="bulk"]) #12
mean(fungicounts$csswp[fungicounts$type=="sentinel"]) #8.5
knitr::kable(c(ttestfungicsswp$statistic, ttestfungicsswp$p.value, 
               ttestfungicsswp$estimate, ttestfungicsswp$conf.int), format = "simple")


#CS
ttestbaccs <- t.test(baccounts$cs[baccounts$type=="bulk"], 
            y = baccounts$cs[baccounts$type=="sentinel"],  
            paired = TRUE)
knitr::kable(c(ttestbaccs$statistic, ttestbaccs$p.value, 
               ttestbaccs$estimate, ttestbaccs$conf.int), format = "simple")


#CSSwP
ttestbaccsswp <- t.test(baccounts$csswp[baccounts$type=="bulk"], 
            y = baccounts$csswp[baccounts$type=="sentinel"],  
            paired = TRUE)
knitr::kable(c(ttestbaccsswp$statistic, ttestbaccsswp$p.value, 
               ttestbaccsswp$estimate, ttestbaccsswp$conf.int), format = "simple")

mean(baccounts$csswp[baccounts$type=="bulk"]) # 6
mean(baccounts$csswp[baccounts$type=="sentinel"]) # 15

