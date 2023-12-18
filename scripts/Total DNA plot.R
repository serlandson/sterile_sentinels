#################################
# Total DNA plot
# load sample data from both 2020 AND 2021
setwd("../OneDrive - USDA/sterile sentinels/sr analysis/")
samples20 <- readxl::read_xlsx("../Sterile Sentinel Datasheet 2020.xlsx")
samples21 <- readxl::read_xlsx("../Sterile Sentinel Datasheet 2021.xlsx")

# remove mislabled sample T1A6 - its acually a bulk sample
samples21 <- samples21[-17, ]

colnames(samples20) <- samples20[1, ]
colnames(samples21) <- samples21[1, ]

#subset
samples20 <- samples20[-1, c(1:4, 9:11, 15,17,19)]
samples21 <- samples21[-1, c(1:3,5, 10:12, 17,19,21)]

#add year column
samples20$year <- rep("2020", length(samples20$Plot))
samples21$year <- rep("2021", length(samples21$Plot))

# rename columns
colnames(samples20) <- c("plot", "type", "time", "sample_id", "tin", "wet", "dry", "r1conc", "r2conc", "r3conc", "year")
colnames(samples21) <- c("plot", "type", "time", "sample_id", "tin", "wet", "dry", "r1conc", "r2conc", "r3conc", "year")

# combine both years
samps <- rbind(samples20, samples21)
samps <- samps[!is.na(samps$type), ]
samps <- samps[samps$type %in% c("bag", "bulk"), ]


samps <- as.data.frame(samps)
samps$type <- factor(samps$type, labels = c("sterile sentinels", "bulk soil"))
str(samps)
samps[, 5:10] <- sapply(samps[,5:10], FUN = as.numeric)
# get proportion water content
samps$SWC <- ((samps$wet-samps$tin) - (samps$dry-samps$tin))/(samps$dry-samps$tin)
#get mean dna conc
samps$meanDNAconc <- apply(samps[,8:10], FUN = mean, MARGIN = 1)
# get total dna on a dry soil basis
# dnaconc ng/uL*100uL = total ng dna per 0.25 g soil
# 1-SWC * 0.25 g soil = dry soil g - proportion of soil that is dry soil, multiply by the 0.25 g soil to get approx. amount of dry soil used.
# convert to ug(total dna*.001) / dry soil g = total ug DNA per g dry soil
samps$totalDNAug.drysoilg <- (samps$meanDNAconc*100)*0.001/((1-samps$SWC)*0.25)
samps$totalDNAug.drysoilg[is.na(samps$totalDNAug.drysoilg)] <- 0

# get weeks
samps$week <- as.numeric(as.character(factor(samps$time, levels = c(0,1,2,3,4,5), labels=c(0,1,2,4,8,12))))
weekmeans <- rbind(aggregate(totalDNAug.drysoilg~week, FUN = mean, data = samps, subset = samps$type == "sterile sentinels"), aggregate(totalDNAug.drysoilg~week, FUN = mean, data = samps, subset = samps$type == "bulk soil"))
weekmeans

# make samples easier to subset for linear and non linear models
samps$type_year <- paste0(samps$type, " ", samps$year)
ss20 <- nls()

bs20 <- lm(totalDNAug.drysoilg ~ week, data = samps, subset = samps$type_year == "bulk soil 2020")
bs20
summary(bs20)


bs21 <- lm(totalDNAug.drysoilg ~ week, data = samps, subset = samps$type_year == "bulk soil 2021")
bs21
summary(bs21)

str(samps)
ss20 <- nls(totalDNAug.drysoilg ~ week, start = c(a=1, b=1),
            data = samps, 
            subset = samps$type_year == "sterile sentinels 2020")

ss21 <- nls(totalDNAug.drysoilg ~ a*week/(b+week), start = c(a=0.1, b=0.1),
            data = samps, 
            subset = samps$type_year == "sterile sentinels 2021")
ss21
summary(ss21)

rss <- sum(residuals(ss21)^2)
tss <- sum((samps[samps$type_year == "sterile sentinels 2021", ]$totalDNAug.drysoilg - mean(samps[samps$type_year == "sterile sentinels 2021", ]$totalDNAug.drysoilg))^2)
1-(rss/tss)

plot(x=samps[samps$type_year == "sterile sentinels 2021", ]$week, y=samps[samps$type_year == "sterile sentinels 2021", ]$totalDNAug.drysoilg)
curve(x*10.6/(x+1.13), from = 0,12, add = T)
#get some estimation of goodness of fit
cor(y,predict(m))

library(tidyverse)
## DNA
ggplot(samps)+
  geom_point(aes(x = week, y = totalDNAug.drysoilg, color = as.factor(year)), size = 3)+
  geom_smooth(aes(x = week, y = totalDNAug.drysoilg, group = year, color = year), span = 1, data = subset(samps, type == "sterile sentinels"), method = "loess", se = F)+
  geom_smooth(aes(x = week, y = totalDNAug.drysoilg, group = year, color = year), data = subset(samps, type == "bulk soil"), method = "lm", se = F)+
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50 ))+
  scale_x_continuous(breaks = c(0,1,2,4,8,12))+
  scale_color_manual(values = c("gray","black"))+
  facet_grid(cols = vars(type), scales = "free")+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 12), 
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.title.y = element_text(size = 12, vjust = 4),
        legend.text = element_text(size = 12), legend.title = element_text(size =12),
        plot.margin = margin(12,12,12,12))+
  labs(y = expression(paste("Total DNA  ", over(paste(" \U003BC","g DNA"), "g dry soil"))), 
       x = "Sampling Week", 
       color = "year")

# guides(color = guide_legend(
#   override.aes=list(shape = 15)))

anova(lm(DNA_ng.ul ~ week*percent_wc*type, data = samps))
summary(lm(DNA_ng.ul ~ week*type, data = samps))

mean(bag$DNA_ng.ul); sd(bag$DNA_ng.ul)
mean(bulk$DNA_ng.ul);sd(bulk$DNA_ng.ul)

############
ggplot(bulk)+
  geom_point(aes(x = week, y = DNA_ng.ul, color = as.factor(week), shape = as.factor(year)), size = 3)+
  theme_minimal()+
  scale_color_brewer(palette = "Spectral")+
  labs(y = "extracted DNA (ng uL)", x = "% Water Content", title = "DNA recovered from sentinels and bulk soil")

anova(lm(DNA_ng.ul ~ week, data = bag))

lm.1 <- lm(DNA_ng.ul ~ week + I(week^2) + I(week^3), data = bag)
lm.s <- step(lm.1)

summary(lm(DNA_ng.ul ~ log(week), data = bag))
summary(lm(DNA_ng.ul ~ week + I(week^2), data = bag))


plot(DNA_ng.ul ~ log(week), data = bag)
plot(DNA_ng.ul ~ week + I(week^2), data = bag)

summary(lm(DNA_ng.ul ~ week, data = bulk))
