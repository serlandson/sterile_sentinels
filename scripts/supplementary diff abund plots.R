library(ggplot2)
library(RColorBrewer)

plotdf2020 <- read.csv("differential abundance data tables/trait differential abundance tables/motility diff abund bacteria 2020 bag v bulk.csv")
plotdf2021 <- read.csv("differential abundance data tables/trait differential abundance tables/motility diff abund bacteria 2021 bag v bulk.csv")


plotdf2020 <- read.csv("differential abundance data tables/trait differential abundance tables/guild diff abund eukaryotes 2020 bag v bulk.csv")
plotdf2021 <- read.csv("differential abundance data tables/trait differential abundance tables/guild diff abund eukaryotes 2021 bag v bulk.csv")


plotdf2020 <- read.csv("filtered differential abundance data tables/trait differential abundance tables/bulk_2020_fungi_diffabund_Guild_Rotation.csv")
plotdf2021 <- read.csv("filtered differential abundance data tables/trait differential abundance tables/bulk_2021_fungi_diffabund_Guild_Rotation.csv")

plotdf2020 <- read.csv("filtered differential abundance data tables/trait differential abundance tables/sentinels_2020_fungi_diffabund_Guild_Rotation.csv")
plotdf2021 <- read.csv("filtered differential abundance data tables/trait differential abundance tables/sentinels_2021_fungi_diffabund_Guild_Rotation.csv")


plotdf2020 <- read.csv("filtered differential abundance data tables/trait differential abundance tables/motility diff abund bacteria 2020 CS v CSSwP.csv")
plotdf2021 <- read.csv("filtered differential abundance data tables/trait differential abundance tables/motility diff abund bacteria 2021 CS v CSSwP.csv")


plotdf2020$year <- rep("2020", length(plotdf2020$rab.all))
plotdf2021$year <- rep("2021", length(plotdf2021$rab.all))

plotdf <- rbind(plotdf2020, plotdf2021)

plotdfeffect <- plotdf[abs(plotdf$effect)>1, ]

# motility
plotdf$week <- gsub(pattern = "time", replacement = "",x = plotdf$type)
plotdf$week <- factor(x = plotdf$week, labels = c(1,2,4,8,12))
plotdf$motility <- plotdf$species
plotdf$motility <- gsub(pattern = "yes", replacement = "motile", x = plotdf$motility)
plotdf$motility <- gsub(pattern = "no", replacement = "non-motile", x = plotdf$motility)

cols <- brewer.pal(9, name = "Blues")
cols <- cols[5:9]

alphaeffect <- abs(plotdf$effect)>1
alphaeffect[alphaeffect==TRUE] <- 1
alphaeffect[alphaeffect==FALSE] <- 1/3

ggplot(plotdf)+
  geom_point(aes(y = motility, 
                 x = effect, 
                 color = week,
                 shape = year),
                 alpha = alphaeffect, size = 3)+
  theme_light()+
  scale_color_manual(values = cols, 
                     name = c("sampling\nweek"),
                     labels = c("1","2","4","8","12"))+
  geom_vline(aes(xintercept = 0))+
  geom_vline(xintercept=-1,linetype = 2, alpha = 0.5 )+
  geom_vline(xintercept= 1,linetype = 2, alpha = 0.5 )+
  guides(alpha = "none")+
  labs(x = "effect size",
       y = "potential motility type")
ggsave("motility diff abund.tiff", device = "tiff", units = "in",
       height = 4, width = 6, plot = last_plot())

### 
# Fungal Guilds

plotdf2020$year <- rep("2020", length(plotdf2020$rab.all))
plotdf2021$year <- rep("2021", length(plotdf2021$rab.all))

plotdf <- rbind(plotdf2020, plotdf2021)

plotdfeffect <- plotdf[abs(plotdf$effect)>1, ]
write.csv(plotdfeffect, "fungi guild diff abund sample type.csv", quote = F, row.names = F)

plotdf$week <- gsub(pattern = "time", replacement = "",x = plotdf$type)
plotdf$week <- factor(x = plotdf$week, labels = c(1,2,4,8,12))

cols <- brewer.pal(9, name = "Blues")
cols <- cols[5:9]

ggplot(plotdf)+
  geom_point(aes(y = species, 
                 x = effect, 
                 color = week,
                 size = rab.all, 
                 shape = year,
                 alpha = abs(effect)>1), size = 4)+
#scale_size(name = "median clr\nabundance")+
theme_light()+
  #scale_color_brewer(palette = "Blues")+
  scale_color_manual(values = cols, 
                     name = c("sampling\nweek"),
                     labels = c("1","2","4","8","12"))+
  
  # geom_vline(aes(xintercept = 0))+
  scale_shape_manual(values = c(17,19,2,1),
                     name = c("year"))+
  geom_vline(xintercept=-1,linetype = 2, alpha = 0.5 )+
  geom_vline(xintercept= 1,linetype = 2, alpha = 0.5 )+
  labs(x = "effect size (sentinel vs. bulk soil)",
       y = "differentially abundant taxa")+
  guides(alpha = "none")
ggsave("fungal guilds diff abund.tiff", device = "tiff", units = "in",
       height = 4, width = 6, plot = last_plot())
