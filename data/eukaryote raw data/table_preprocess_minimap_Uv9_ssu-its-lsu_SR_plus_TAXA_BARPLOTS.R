library(tidyverse)
setwd(dir = "/Users/sonya.erlandson/OneDrive - USDA/sterile sentinels/sr analysis/")
## Full ITS ##
tab <- read.delim("minimap_ssu-its-lsu_sr/euk_read_count_table_2020.tsv", header = TRUE)
taxonomy <- read.delim("minimap_ssu-its-lsu_sr/euk_taxonomy_table_2020.tsv", header = TRUE)
samples <- read.csv("../minimap2 statistical analyses/data tables/precursor data tables/sterile_sentinels_sample_data_2020.csv")

# samples and table order
samples$fungi_barcode <- gsub("BC", "barcode", samples$Euk.BC.Well)
samples <- samples[order(samples$fungi_barcode), ]

summary(colnames(tab)==taxonomy$taxon_name_sh)

# get control samples
bag_bulk <- samples$fungi_barcode[c(which(samples$type == "bag"), which(samples$type == "bulk"))]
controls <- samples[!samples$fungi_barcode %in% bag_bulk, ]

control.comm <- tab[row.names(tab) %in% controls$fungi_barcode, ]
control.comm <- control.comm[ ,colSums(control.comm)>0]
control.comm <-as.data.frame(t(control.comm))
control.comm$taxa <- taxonomy$taxon_name_sh[taxonomy$taxon_name_sh %in% row.names(control.comm)]
control.comm$size <- colSums(tab)[colnames(tab) %in% row.names(control.comm)]

#write.csv(control.comm, "mm_Uv9_ssu-its-lsu_SR_control_mock_communities_2020.csv", quote = FALSE)

# remove controls and empty barcodes from table
tab.sub <- tab[row.names(tab) %in% bag_bulk, ]

sort(rowSums(tab.sub))
hist(rowSums(tab.sub))
summary(colSums(tab.sub))
hist(colSums(tab.sub))

# ways to think about abundance filtering
relabund <- colSums(tab.sub)/sum(colSums(tab.sub))

summary(taxonomy$taxon_name_sh == colnames(tab.sub))
taxonomy$size <- colSums(tab.sub)
taxonomy$relabund <- relabund

patab <- tab.sub
patab[patab>0] <- 1
hist(colSums(patab))
summary(colSums(patab)>2)

# mean abundance without zero
meanwithoutzero <- sapply(tab.sub, function(x) mean(x[x != 0], na.rm = TRUE))
taxonomy$means_no0 <- meanwithoutzero


# remove taxa with less than 3 observations (reads)
tab.sub <- tab.sub[, colSums(tab.sub)>=3]


# check controls
control.comm$taxa
testcont <- tab.sub[ ,colnames(tab.sub) %in% row.names(control.comm)] # taxa in samples and controls
control.overlap <- control.comm[match(colnames(testcont), row.names(control.comm)), ] #fix column names?
colSums(testcont) # 41 taxa are still in the table after removing the control samples
rowSums(control.overlap[ ,1:5]); control.overlap$size # taxa - just remove the cryptococcus and saccharamyces?


# remove barcode 37 and barcode67
filt.tab <- tab.sub[!row.names(tab.sub) %in% c("barcode37", "barcode67"), ]
summary(colSums(filt.tab)>0)

#for 2020
filt.tab <- tab.sub

filt.tab <- filt.tab[ ,colSums(filt.tab)>=3]

#subset taxonomy
tax.sub <- taxonomy[taxonomy$taxon_name_sh %in% colnames(filt.tab), ]
#remove unite formatting
tax.sub$taxonomy <- str_remove_all(tax.sub$taxonomy, "[a-z]{1}__")

# separate taxa
tax.sub <- separate(tax.sub, col = taxonomy, sep = ";", into = c("kingdom", "phylum", "class", "order", "family", "genus","species"))
unique(tax.sub$kingdom)


write.csv(filt.tab, "mm_Uv9_ssu-its-lsu_SR_table_2020.csv", quote = F)
write.csv(tax.sub, "mm_Uv9_ssu-its-lsu_SR_taxonomy_2020.csv", quote = F)

write.csv(filt.tab, "mm_Uv9_ssu-its-lsu_SR_table_2021.csv", quote = F)
write.csv(tax.sub, "mm_Uv9_ssu-its-lsu_SR_taxonomy_2021.csv", quote = F)

# load from eukaryotes
filt.tab <- read.csv("other data tables/mm_Uv9_ssu-its-lsu_SR_table_2020.csv", row.names = 1)
tax.sub <- read.csv("other data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_2020.csv", row.names = 1)

summary(colnames(filt.tab)==tax.sub$taxon_name_sh)


filt.tab <- read.csv("other data tables/mm_Uv9_ssu-its-lsu_SR_table_2021.csv", row.names = 1)
tax.sub <- read.csv("other data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_2021.csv", row.names = 1)

# glomeromycota
glomotus <- tax.sub[grep("Glomeromycota", tax.sub$phylum), ]
glom <- filt.tab[ ,colnames(filt.tab) %in% glomotus$taxon_name_sh]

write.csv(glomotus, "Glomeromycota_mm_Uv9_ssu-its-lsu.csv")

# fungi only
# only keep fungi that are in 3 samples or more.
fungi.tab <- filt.tab[ , grep("Fungi", tax.sub$kingdom)]
fungi.tax <- tax.sub[grep("Fungi", tax.sub$kingdom), ]

sort(rowSums(fungi.tab))
summary(colSums(fungi.tab))

fungipa <- fungi.tab
fungipa[fungipa>0] <- 1

#fungi.tab <- fungi.tab[ ,colSums(fungipa)>=3]
#fungi.tax <- fungi.tax[fungi.tax$taxon_name_sh %in% colnames(fungi.tab), ]

# add traits
traits <- readxl::read_xlsx("fungal traits.xlsx")
traits <- as.data.frame(traits)

tax.sub$primary_lifestyle <- traits[match(tax.sub$genus, traits$GENUS), 8]
tax.sub$secondary_lifestyle <- traits[match(tax.sub$genus, traits$GENUS), 9]
tax.sub$decay_type <- traits[match(tax.sub$genus, traits$GENUS), 13]
tax.sub$growth_form <- traits[match(tax.sub$genus, traits$GENUS), 18]

tax <- tax.sub
#animal parasite
reassign <- which(tax$primary_lifestyle == "animal_parasite")
newnames <- tax$secondary_lifestyle[reassign]
tax$primary_lifestyle[reassign] <- newnames

# foliar_endophyte
reassign <- which(tax$primary_lifestyle == "foliar_endophyte")
newnames <- tax$secondary_lifestyle[reassign]
tax$primary_lifestyle[reassign] <- newnames


# epiphyte
reassign <- which(tax$primary_lifestyle == "epiphyte")
newnames <- tax$secondary_lifestyle[reassign]
tax$primary_lifestyle[reassign] <- newnames


fungi.tax <- tax[tax$kingdom=="Fungi", ]


write.csv(fungi.tab, "data tables/mm_Uv9_ssu-its-lsu_SR_table_FUNGI_2020.csv", quote = F)
write.csv(fungi.tax, "data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_FUNGI_2020.csv", quote = F)


write.csv(fungi.tab, "data tables/mm_Uv9_ssu-its-lsu_SR_table_FUNGI_2021.csv", quote = F)
write.csv(fungi.tax, "data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_FUNGI_2021.csv", quote = F)


##### nematodes
nematode_Genera <- unique(tax[tax$phylum == "Nematoda", ]$genus)

write.csv(nematode_Genera, "nematode_genera2.csv", quote = F, row.names = F)
.##### make taxa barplots ##########################################
library(tidyverse)
comm <- fungi.tab
tax <- fungi.tax
colnames(comm) <- tax$taxon_name # rename columns to remove sh

tax$taxonomy <- str_remove_all(tax$taxonomy, "[a-z]{1}__") # remove unite formatting
tax <- separate(tax, col = taxonomy, sep = ";", into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
# separate taxonomy column 

comm$barcode <- row.names(comm)
plotdf <- as.data.frame(pivot_longer(comm, cols = !barcode, names_to = "species", values_to = "counts"))
plotdf <- plotdf[plotdf$counts > 0, ] # remove zeros


# make columns for different taxonomic levels
summary(tax$species[match(plotdf$species, tax$species)] == plotdf$species) # check that they match
plotdf$phylum <- as.factor(tax$phylum[match(plotdf$species, tax$species)])
plotdf$barcode <- as.factor(plotdf$barcode)
unique(plotdf$phylum)

plotdf <- plotdf[plotdf$counts > 10, ]

# remove taxa that don't have a phyla assignment
#plotdf <- plotdf[!plotdf$phylum == "", ]
plotdf$phylum <- droplevels(plotdf$phylum)

# fix counts and check df structure
plotdf$counts <- as.numeric(plotdf$counts)
plotdf <- as.data.frame(plotdf)
str(plotdf)

# for barplot
plotdf$week <- as.factor(samples$week[match(plotdf$barcode, samples$fungi_barcode)])
plotdf$time <- as.factor(samples$time[match(plotdf$barcode, samples$fungi_barcode)])
plotdf$type <- as.factor(samples$type[match(plotdf$barcode, samples$fungi_barcode)])
plotdf$rot <- as.factor(samples$rotation[match(plotdf$barcode, samples$fungi_barcode)])
plotdf$week_rot <- paste0(plotdf$week, "_", plotdf$rot)

# bags/bulk
bagbarcodes <- samples$fungi_barcode[samples$type == "bag"]
bulkbarcodes <- samples$fungi_barcode[samples$type == "bulk"]
plotdf_bags <- plotdf[plotdf$barcode %in% bagbarcodes, ]
plotdf_bulk <- plotdf[plotdf$barcode %in% bulkbarcodes, ]

unique(plotdf_bags$phylum)
plotdf_subset <- plotdf_bags[ plotdf_bags$week == 1 | plotdf_bags$week == 12, ]
plotdf_subset <- plotdf_bags[ plotdf_bags$time == 1 | plotdf_bags$time == 5, ]
plotdf_subset$phylum <- droplevels(plotdf_subset$phylum)
unique(plotdf_subset$phylum)

library(viridis)

# plot bags and bulk
phycols <- viridis(10, option = "turbo", direction = -1)
ggplot(plotdf_subset, aes(x=as.factor(week_rot), y = counts))+
  geom_bar(aes(fill = phylum), stat = "identity", position = "fill")+
  labs(x = "Week Sampled After Placement", 
       y = "Fungal Phylum Relative Abundances", 
       title = "Minimap2 Fungal Phyla Sentinels")+
  scale_fill_manual(values = phycols)+
  theme_minimal()


###############################################################
# GLOMERO
##### make taxa barplots
library(tidyverse)
comm <- glom
tax <- glomotus
colnames(comm) <- tax$taxon_name # rename columns to remove sh

tax$taxonomy <- str_remove_all(tax$taxonomy, "[a-z]{1}__") # remove unite formatting
tax <- separate(tax, col = taxonomy, sep = ";", into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
# separate taxonomy column 

comm$barcode <- row.names(comm)
plotdf <- as.data.frame(pivot_longer(comm, cols = !barcode, names_to = "species", values_to = "counts"))
plotdf <- plotdf[plotdf$counts > 0, ] # remove zeros


# make columns for different taxonomic levels
summary(tax$species[match(plotdf$species, tax$species)] == plotdf$species) # check that they match

summary(tax$taxon_name_sh[match(plotdf$species, tax$taxon_name_sh)] == plotdf$species) # check that they match
plotdf$kingdom<- as.factor(tax$kingdom[match(plotdf$species, tax$taxon_name_sh)])
plotdf$phylum <- as.factor(tax$phylum[match(plotdf$species, tax$taxon_name_sh)])
plotdf$class <- as.factor(tax$class[match(plotdf$species, tax$taxon_name_sh)])
plotdf$order <- as.factor(tax$order[match(plotdf$species, tax$taxon_name_sh)])
plotdf$family <- as.factor(tax$family[match(plotdf$species, tax$taxon_name_sh)])
plotdf$genus <- as.factor(tax$genus[match(plotdf$species, tax$taxon_name_sh)])
plotdf$barcode <- as.factor(plotdf$barcode)
unique(plotdf$family)
unique(plotdf$genus)
summary(plotdf$counts)
plotdf <- plotdf[plotdf$counts > 10, ]

# remove taxa that don't have a phyla assignment
#plotdf <- plotdf[!plotdf$phylum == "", ]
plotdf$family <- droplevels(plotdf$family)
plotdf$kingdom <- droplevels(plotdf$kingdom)

# fix counts and check df structure
plotdf$counts <- as.numeric(plotdf$counts)
plotdf <- as.data.frame(plotdf)
str(plotdf)

# for barplot
plotdf$week <- as.factor(samples$week[match(plotdf$barcode, samples$fungi_barcode)])
plotdf$time <- as.factor(samples$time[match(plotdf$barcode, samples$fungi_barcode)])
plotdf$type <- as.factor(samples$type[match(plotdf$barcode, samples$fungi_barcode)])
plotdf$rot <- as.factor(samples$rotation[match(plotdf$barcode, samples$fungi_barcode)])
plotdf$week_rot <- as.factor(paste0(plotdf$week, "_", plotdf$rot))
plotdf$week_rot <- factor(plotdf$week_rot, levels = c("0_CS", "0_CSSwP", 
                                                      "1_CS", "1_CSSwP", 
                                                      "2_CS", "2_CSSwP", 
                                                      "4_CS", "4_CSSwP", 
                                                      "8_CS", "8_CSSwP",
                                                      "12_CS", "12_CSSwP"))

plotdf$week_rot <- factor(plotdf$week_rot, levels = c("0_CS", "0_CSSwP", 
                                                      "1_CS", "1_CSSwP", 
                                                      "2_CS", "2_CSSwP", 
                                                      "3_CS", "3_CSSwP", 
                                                      "4_CS", "4_CSSwP",
                                                      "5_CS", "5_CSSwP"))


# bags/bulk
bagbarcodes <- samples$fungi_barcode[samples$type == "bag"]
bulkbarcodes <- samples$fungi_barcode[samples$type == "bulk"]

plotdf_bags <- plotdf[plotdf$barcode %in% bagbarcodes, ]
plotdf_bulk <- plotdf[plotdf$barcode %in% bulkbarcodes, ]

unique(plotdf_bulk$kingdom)
plotdf_bulk$family <- droplevels(plotdf_bulk$family)

plotdf_bulk$kingdom <- droplevels(plotdf_bulk$kingdom)

unique(plotdf_bags$genus)
plotdf_bags$genus <- droplevels(plotdf_bags$genus)
plotdf_bags$kingdom <- droplevels(plotdf_bags$kingdom)

plotdf_subset <- plotdf_bags[ plotdf_bags$week == 1 | plotdf_bags$week == 12, ]
plotdf_subset$genus <- droplevels(plotdf_subset$genus)
unique(plotdf_subset$genus)

library(viridis)

nofungi <- plotdf[plotdf$kingdom!="Fungi", ]
nema <- plotdf[plotdf$phylum=="Nematoda", ]
alvea <- plotdf[plotdf$kingdom=="Alveolata", ]
strame <- plotdf[plotdf$kingdom=="Stramenopila", ]

metazoa <- plotdf[plotdf$kingdom=="Metazoa", ]
amoeba <- plotdf[plotdf$kingdom=="Amoebozoa", ]
viri <- plotdf[plotdf$kingdom=="Viridiplantae", ]
unique(plotdf$kingdom)

# plot bags and bulk
unique(viri$class)
unique(plotdf$kingdom)
phycols <- viridis(12, option = "turbo", direction = -1)

ggplot(strame, aes(x=as.factor(week_rot), y = counts))+
  geom_bar(aes(fill = family), stat = "identity", position = "fill")+
  labs(x = "Week Sampled After Placement", 
       y = "Relative Abundances", 
       title = "Stramenopila Classes 2020")+
  scale_fill_manual(values = phycols)+
  facet_grid(vars(type))+
  theme_minimal()

