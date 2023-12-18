# code for making heatmaps of most prevalent taxa in each dataset. 
# the plot annotated with guild or phyla is made from two separate plots glued together
# with patchwork, a plot of taxa names+guild or phyla assignment and a plot of relative abundance. 

# be careful to change names or commands for filtering 2020 or 2021 data
# meant to make this a function, but not a functional function yet.
#+ set up
library(tidyverse,  quietly = T, warn.conflicts = F)
library(vegan, quietly = T, warn.conflicts = F)
library(easyCODA)
library(colorspace)
library(patchwork)
library(ggtext)

comm <- read.csv("data tables/mm_Uv9_ssu-its-lsu_SR_counts_table_FUNGI_2020.csv", header = TRUE, row.names = 1, check.names = F)
tax <- read.csv("data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_FUNGI_2020.csv", header = TRUE)
samples <- read.csv("data tables/its_sample_data_2020.csv")

comm <- read.csv("data tables/mm_Uv9_ssu-its-lsu_SR_counts_table_FUNGI_2021.csv", header = TRUE, row.names = 1, check.names = F)
tax <- read.csv("data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_FUNGI_2021.csv", header = TRUE)
samples <- read.csv("data tables/its_sample_data_2021.csv")

fungiHeatMap <- function(community_data, sample_data, taxonomy_data, amplicon_year, Plot_title)
  
  comm <- read.csv(community_data,  header = TRUE, row.names = 1, check.names = F)
  samples <- read.csv(sample_data)
  tax <- read.csv(taxonomy_data, header = TRUE)
  
  # FUNGI 
  samples$ttr <- paste(samples$type, samples$time, samples$rotation, str_sub(samples$plot,1,1), sep = "-") 
  samples$tt <- paste(samples$type, samples$time, samples$rotation, sep = "-") 
  
  # change name from bag to sentinel
  samples$ttr <- gsub("bag", "sentinel", samples$ttr)
  samples$tt <- gsub("bag", "sentinel", samples$tt)
  
  # order samples from 0,bulk ... 12 sentinel
  samples <- samples[order(samples$ttr, decreasing = F), ]
  
  # fungi match the community table sample order with samples table
  commordered <- as.matrix(comm[match(samples$fungi_barcode, row.names(comm)), ])
  
  # re order the species in bulk time 0 from largest to smallest
  commordered <- commordered[ ,order(commordered[1, ], decreasing = T)]
  # average up counts from each rep
  aggdf <- aggregate(as.data.frame(commordered), by = list(samples$tt) , FUN = mean) 
  # get relative abundance per category
  aggrelabund <- aggdf[,-1]/rowSums(aggdf[,-1])
  # separate the group variables into columns.
  aggrelabund <- cbind(aggrelabund, Group.1 = aggdf$Group.1)
  plotdf <- separate(aggrelabund, "Group.1", into = c("type", "week", "rotation"), sep = "-")
  
  # make long format dataframe
  plotdf <- reshape2::melt(aggrelabund)
  plotdf <- separate(plotdf, "Group.1", into = c("type", "week", "rotation"), sep = "-")
  colnames(plotdf) <- c("type", "week", "rotation", "species", "rel.abund")
  
  plotdf$species <- as.character(plotdf$species)
  
  rotationdata <- read.csv("../../sr analysis/differential abundance data tables/all_taxa_Rotation_diff_abund_no_filter.csv", stringsAsFactors = T)
  #filter to effect size greater or less than 1
  rotationdata <- rotationdata[abs(rotationdata$effect)>1 ,]
  #rotationdata$year <- as.factor(rotationdata$year)
  fungi <- rotationdata[str_detect(rotationdata$amplicon, "fungi_2021"), ] # get df of fungi from 2020
  
  sampletypedata <- read.csv("../../sr analysis/differential abundance data tables/all taxa diff abund TYPE.csv", stringsAsFactors = T)
  sampletypedata <- sampletypedata[abs(sampletypedata$effect)>1, ]
  fungitype <- sampletypedata[str_detect(sampletypedata$amplicon, "fungi_2021"), ]
  fungitype <- fungitype[fungitype$tax_level == "species", ]
  
  taxaoptions <- c(fungi$species, fungitype$species)
  taxaoptions <- droplevels(unique(taxaoptions))
  
  diff_rot <- droplevels(unique(fungi$species))
  diff_type <- droplevels(unique(fungitype$species))
  
  prevtaxa <- colSums(comm)/sum(colSums(comm))
  top40 <- names(sort(prevtaxa, decreasing = T)[1:40])
  
  keeptaxa <- taxaoptions[taxaoptions %in% top40]
  keeptaxa <- droplevels(keeptaxa)
  ################################################
  # subset to the most prevalent taxa
  plotdf_filt <- plotdf[plotdf$species %in% keeptaxa, ]
  plotdf_filt$type_rot <- paste0(plotdf_filt$type, " ", plotdf_filt$rotation)
  
  # force the order of species
  bulkcs0 <- plotdf_filt[plotdf_filt$type=="bulk" & plotdf_filt$week == 0 & plotdf_filt$rotation == "CS", ]
  
  bulkcs0$guild <- tax$primary_lifestyle[match(bulkcs0$species,tax$taxon_name_sh)]
  bulkcs0$guild2 <- tax$secondary_lifestyle[match(bulkcs0$species,tax$taxon_name_sh)]
  
  bulkcs0$guild[str_detect(bulkcs0$guild, "saprotroph")] <- "saprotroph"
  bulkcs0$guild[is.na(bulkcs0$guild)] <- "unassigned"
  
  bulkcs0$guild2[str_detect(bulkcs0$guild2, "saprotroph")] <- "saprotroph"
  bulkcs0$guild2[is.na(bulkcs0$guild2)] <- "unassigned"
  bulkcs0$guild2[str_detect(bulkcs0$guild2, "algal_parasite")] <- "unassigned"
  
  
  bulkcs0 <- bulkcs0[order(bulkcs0$rel.abund), ]
  taxaorderedbyguild <- c(which(bulkcs0$guild == "unassigned"),
                          which(bulkcs0$guild == "arbuscular_mycorrhizal"),
                          which(bulkcs0$guild == "mycoparasite"),
                          which(bulkcs0$guild == "plant_pathogen"),
                          which(bulkcs0$guild == "saprotroph"))
  
  bulkcs0 <- bulkcs0[taxaorderedbyguild, ]
  taxaorder <- as.character(bulkcs0$species)

  # get a dataframe of guilds
  guilds <- data.frame(x = rep(0, length(taxaorder)), 
                       species = taxaorder,
                       guild = tax$primary_lifestyle[match(taxaorder, tax$taxon_name_sh)],
                       guild2 = tax$secondary_lifestyle[match(taxaorder, tax$taxon_name_sh)])
  guilds$guild[str_detect(guilds$guild, "saprotroph")] <- "saprotroph"
  guilds$guild[is.na(guilds$guild)] <- "unassigned"
  
  guilds$guild2[str_detect(guilds$guild2, "saprotroph")] <- "saprotroph"
  guilds$guild2[is.na(guilds$guild2)] <- "unassigned"
  guilds$guild2[guilds$guild2 == "algal_parasite" ] <- "unassigned"
  
  guilds2 <- guilds[, c(1,2,4)]
  guilds2$x <- 1
  colnames(guilds2) <- c("x", "species", "guild")
  guilds <- guilds[,c(1,2,3)]

    plotguilds <- rbind(guilds, guilds2)
  
    plotguilds$x <- as.factor(plotguilds$x)
    plotguilds$species <- factor(plotguilds$species, levels = taxaorder)
    plotguilds$guild <- gsub("_", " ", plotguilds$guild)
  
    # reorder taxa
  plotdf_filt$species <- factor(plotdf_filt$species, levels = taxaorder)
  plotdf_filt$week_plot <- factor(plotdf_filt$week, labels = c(0,1,2,4,8,12))
  
  # make new spp labels
  # get groups
  diffrotation <- droplevels(diff_rot[diff_rot %in% taxaorder])
  difftype <- droplevels(diff_type[diff_type %in% taxaorder])
  bothdiff <- droplevels(diffrotation[diffrotation %in% difftype])
  #only in one
  diffrotonly <- diffrotation[!diffrotation %in% bothdiff]
  difftypeonly <- difftype[!difftype %in% bothdiff]
  
  #make new vector
  diffabundloc <- rep("a", length(taxaorder))
  diffabundloc[match(bothdiff, taxaorder)] <- "b"
  diffabundloc[match(diffrotonly, taxaorder)] <- "r"
  diffabundloc[match(difftypeonly, taxaorder)] <- "t"
  
  #newsppnames <- str_sub(taxaorder, 1, (str_locate(taxaorder,"_SH")[,1]-1))
  newsppnames <- str_sub(taxaorder,1, str_locate(taxaorder, "_")[,1]-1)
 
  newnames <- paste0("*",newsppnames,"*"," ","sp", "<sup>",diffabundloc, "</sup>")
  newnames
  # 
  
  fungiheatmap <- ggplot(plotdf_filt)+
    geom_tile(aes(week_plot, species, fill = sqrt(rel.abund)))+
    scale_fill_viridis_c( name = "sqrt of mean\nrelative\nabundance")+
    #ggtitle(Plot_title)+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.x = element_text(size=9, angle = 90),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 11),
          plot.title = element_text(size = 11, hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.position = "right",
          plot.margin = margin(0,0,0,0))+
    xlab("sampling week")+
    facet_grid(cols=vars(type_rot), scales = "free")
  fungiheatmap
  #"#DB9D85" "#ABB065" "#5CBD92" "#4CB9CC" "#ACA4E2" "#E093C3"
  #orange     yellow     green     teal      blue      purple
  # foliar fungdec       myco      plant     root      sapro
  # 2021
  guildcolors <- c("#76B2DB", qualitative_hcl(n = 6, palette = "Dynamic"), "#BBBBBB")
  # 2020
  guildcolors <- c(qualitative_hcl(n = 6, palette = "Dynamic"), "#BBBBBB")
  
guildplot <-  ggplot(plotguilds)+
    geom_tile(aes(x = x, y = species, fill = guild))+
    #scale_fill_brewer(palette = "Set3")+
    scale_fill_manual(values = guildcolors)+
    #scale_fill_discrete_qualitative(palette = "Dark 3")+
    scale_y_discrete(labels = newnames, expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 9),
          legend.position = "left", 
          legend.key = element_blank(),
          legend.key.size = unit(1, "line") , 
          plot.margin = margin(0,0,0,0))+
    xlab("sampling week")+ ylab("differentially abundant taxa")+
  theme(axis.text.y = element_markdown())

  
  layout <- "
ABBBB
"
  
  fungiplot <- guildplot+fungiheatmap+plot_layout(design = layout, widths = c(0.25,4), guides = "collect")
  fungiplot

  fungi2020heatmap <- fungiplot
  ggsave("fungi2020_heatmap_supp_fig_12.tiff", plot = fungi2020heatmap, width = 8.5, height = 6, units = "in")
  
  fungi2021heatmap <- fungiplot
  ggsave("fungi2021_heatmap_supp_fig_12.tiff", plot = fungi2020heatmap, width = 8.5, height = 6, units = "in")

# fungi2020heatmap <- fungiHeatMap("data tables/mm_Uv9_ssu-its-lsu_SR_counts_table_FUNGI_2020.csv",
#                                  "data tables/its_sample_data_2020.csv",
#                                  "data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_FUNGI_2020.csv", 
#                                  "fungi_2020", "Fungi 2020")
# 
# fungi2021heatmap <- fungiHeatMap("data tables/mm_Uv9_ssu-its-lsu_SR_counts_table_FUNGI_2021.csv",
#                                  "data tables/its_sample_data_2021.csv",
#                                  "data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_FUNGI_2021.csv", 
#                                  "fungi_2021", "Fungi 2021")
# fungi2020heatmap
# fungi2021heatmap

###############################################################################
#+ set up

# load data from bac 2020 and bac 2021 scripts

bacteriaHeatMap <- function(community_data, sample_data, taxonomy_data, amplicon_year, Plot_title)
  
  comm <- read.csv(community_data,  header = TRUE, row.names = 1, check.names = F)
  samples <- read.csv(sample_data)
  tax <- read.csv(taxonomy_data, header = TRUE)
  
  # bacteria
  samples$ttr <- paste(samples$type, samples$time, samples$rotation, str_sub(samples$plot,1,1), sep = ",") 
  samples$tt <- paste(samples$type, samples$time, samples$rotation, sep = ",") 
  
  # change name from bag to sentinel
  samples$ttr <- gsub("bag", "sentinel", samples$ttr)
  samples$tt <- gsub("bag", "sentinel", samples$tt)
  
  # order samples from 0,bulk ... 12 sentinel
  samples <- samples[order(samples$ttr, decreasing = F), ]
  
  # match the community table sample order with samples table
  commordered <- as.matrix(comm[match(samples$bacteria_barcode, row.names(comm)), ])
  
  # re order the species in bulk time 0 from largest to smallest
  commordered <- commordered[ ,order(commordered[1, ], decreasing = T)]
  # average up counts from each rep
  aggdf <- aggregate(as.data.frame(commordered), by = list(samples$tt) , FUN = mean) 
  # get relative abundance per group
  aggrelabund <- aggdf[,-1]/rowSums(aggdf[,-1])
  # separate the group variables into columns.
  aggrelabund <- cbind(aggrelabund, Group.1 = aggdf$Group.1)
  plotdf <- separate(aggrelabund, "Group.1", into = c("type", "week", "rotation"), sep = ",")
  
  # make long format dataframe
  plotdf <- reshape2::melt(aggrelabund)
  plotdf <- separate(plotdf, "Group.1", into = c("type", "week", "rotation"), sep = ",")
  colnames(plotdf) <- c("type", "week", "rotation", "species", "rel.abund")
  
  plotdf$species <- as.character(plotdf$species)
  
  # get list of differentially abundant taxa from sentinels vs. bulk and CS vs. CSSwP t-tests
  # also get most prevalent taxa. 
  # only keep most prevalent taxa that are differentially abundant.
  rotationdata <- read.csv("../../../sr analysis/differential abundance data tables/all_taxa_Rotation_diff_abund_no_filter.csv", stringsAsFactors = T)
  #filter to effect size greater or less than 1
  rotationdata <- rotationdata[abs(rotationdata$effect)>1 ,]
  bacteria <- rotationdata[str_detect(rotationdata$amplicon, "bac_2021"), ]
  
  sampletypedata <- read.csv("../../../sr analysis/differential abundance data tables/all taxa diff abund TYPE.csv", stringsAsFactors = T)
  sampletypedata <- sampletypedata[abs(sampletypedata$effect)>1, ]
  bactype <- sampletypedata[sampletypedata$amplicon == "bacteria_2021", ]
  bactype <- bactype[bactype$tax_level == "species", ]
  
  taxaoptions <- c(bacteria$species, bactype$species)
  taxaoptions <- droplevels(unique(taxaoptions))
  
  diff_rot <- droplevels(unique(bacteria$species))
  diff_type <- droplevels(unique(bactype$species))
  
  prevtaxa <- colSums(comm)/sum(colSums(comm))
  top40 <- names(sort(prevtaxa, decreasing = T)[1:40])
  
  keeptaxa <- taxaoptions[taxaoptions %in% top40]
  keeptaxa <- droplevels(keeptaxa)
  ################################################
  # subset to the most prevalent taxa
  plotdf_filt <- plotdf[plotdf$species %in% keeptaxa, ]
  
  plotdf_filt$type_rot <- paste0(plotdf_filt$type, " ", plotdf_filt$rotation)
  plotdf_filt$week_plot <- factor(plotdf_filt$week, labels = c(0,1,2,4,8,12))
  
  # force the order of species - order by abundance in bulkcs0, grouped by phylum
  bulkcs0 <- plotdf_filt[plotdf_filt$type=="bulk" & plotdf_filt$week == 0 & plotdf_filt$rotation == "CS", ]
  bulkcs0$phyla <- taxtrait[match(bulkcs0$species, taxtrait$species), c(7)]
  bulkcs0 <- bulkcs0[order(bulkcs0$rel.abund), ]
  unique(bulkcs0$phyla)
  
  #2020
  taxaorderedbyphyla <- c(which(bulkcs0$phyla == "Proteobacteria"),
                          which(bulkcs0$phyla == "Nitrospirae"),
                          which(bulkcs0$phyla == "Acidobacteria"),
                          which(bulkcs0$phyla == "Planctomycetes"),
                          which(bulkcs0$phyla == "Verrucomicrobia"),
                          which(bulkcs0$phyla == "Gemmatimonadetes"),
                          which(bulkcs0$phyla == "Bacteroidetes"),
                          which(bulkcs0$phyla == "Actinobacteria"))
  # 2021
  taxaorderedbyphyla <- c(which(bulkcs0$phyla == "Proteobacteria"),
                          which(bulkcs0$phyla == "Nitrospirae"),
                          which(bulkcs0$phyla == "Acidobacteria"),
                          which(bulkcs0$phyla == "Planctomycetes"),
                          which(bulkcs0$phyla == "Verrucomicrobia"),
                          which(bulkcs0$phyla == "Gemmatimonadetes"),
                          which(bulkcs0$phyla == "Bacteroidetes"),
                          which(bulkcs0$phyla == "Firmicutes"),
                          which(bulkcs0$phyla == "Actinobacteria"))
  
  
  
  
  bulkcs0 <- bulkcs0[taxaorderedbyphyla, ]
  taxaorder <- as.character(bulkcs0$species)
  #order plotting dataframes
 str(plotdf_filt)
  plotdf_filt$ordspecies <- factor(plotdf_filt$species, levels = taxaorder)
  phylaplotdf <- data.frame(x = rep(0, length(taxaorder)), 
                       species = factor(taxaorder, levels = taxaorder),
                       phyla = bulkcs0$phyla)
  str(phylaplotdf)
  # make new spp labels
  # get groups
  diffrotation <- droplevels(diff_rot[diff_rot %in% taxaorder])
  difftype <- droplevels(diff_type[diff_type %in% taxaorder])
  bothdiff <- droplevels(diffrotation[diffrotation %in% difftype])
  #only in one
  diffrotonly <- diffrotation[!diffrotation %in% bothdiff]
  difftypeonly <- difftype[!difftype %in% bothdiff]
  
  #make new vector
  diffabundloc <- rep("n", length(taxaorder))
  diffabundloc[match(bothdiff, taxaorder)] <- "b"
  diffabundloc[match(diffrotonly, taxaorder)] <- "r"
  diffabundloc[match(difftypeonly, taxaorder)] <- "t"
  
  #bacsppnames <- gsub(pattern = " ",replacement = "_", taxaorder)
  taxaorder[19] <- "Betaproteobacteria sp. GR16-43" # fix how long this name is 2020 only
  newnames <- paste0("*",taxaorder,"*", "<sup>",diffabundloc, "</sup>")
  newnames
  #newnames <- paste0(test, "^","{",diffabundloc, "}")

 bacteriaheatmap <-  ggplot(plotdf_filt)+
    geom_tile(aes(week_plot, ordspecies, fill = sqrt(rel.abund)))+
    scale_fill_viridis_c( name = "sqrt of mean\nrelative\nabundance")+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.x = element_text(size=9, angle = 90),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 11),
          plot.title = element_text(size = 11, hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.position = "right",
          plot.margin = margin(0,0,0,0))+
    xlab("sampling week")+
    facet_grid(cols=vars(type_rot), scales = "free")
  bacteriaheatmap
  
  # 2020
  phylacolors <- c(qualitative_hcl(n = 9, palette = "Dark2"))
  phylacolors <- phylacolors[-4]
  #2021
  phylacolors <- c(qualitative_hcl(n = 9, palette = "Dark2"))
  
  phylaplot <-  ggplot(phylaplotdf)+
    geom_tile(aes(x = x, y = species, fill = phyla))+
    scale_fill_manual(values = phylacolors)+
    scale_y_discrete(labels = newnames, expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 9),
          legend.position = "left", 
          legend.key = element_blank(),
          legend.key.size = unit(1, "line") , 
          plot.margin = margin(0,0,0,0))+
    xlab("sampling week")+ ylab("differentially abundant taxa")+
    theme(axis.text.y = element_markdown())
  
  phylaplot
  layout <- "
ABBBB
"
  
  bacplot2020 <- phylaplot+bacteriaheatmap+plot_layout(design = layout, widths = c(0.25,4), guides = "collect")
  
  
bacplot2020
ggsave("bac2020_heatmap_supp_fig_12.tiff", plot = bacplot2020, width = 9, height = 7, units = "in")

bacplot2021 <- phylaplot+bacteriaheatmap+plot_layout(design = layout, widths = c(0.25,4), guides = "collect")
ggsave("bac2021_heatmap_supp_fig_12.tiff", plot = bacplot2021, width = 9, height = 7, units = "in")
# bac2020heatmap <- bacteriaHeatMap("data tables/emu_16S_counts_2020.csv",
#                                   "data tables/16S_sample_data_2020.csv",
#                                   "data tables/emu_16S_taxonomy_2020.csv",
#                                   "bacteria_2020", "Bacteria 2020")
# 
# 
# bac2021heatmap <- bacteriaHeatMap("data tables/emu_16S_counts_2021.csv",
#                                   "data tables/16S_sample_data_2021_subset.csv",
#                                   "data tables/emu_16S_taxonomy_2021.csv",
#                                   "bacteria_2021", "Bacteria 2021")






