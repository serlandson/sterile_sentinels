#+ load data
library(tidyverse)

df1 <- read.csv("fungi_2020_adonis_results.csv")
df2 <- read.csv('fungi_2021_adonis_results.csv')
df3 <- read.csv('bacteria_2020_adonis_results.csv')
df4 <- read.csv("bacteria_2021_adonis_results.csv")

r2 <- rbind(df1,df2,df3,df4)

r2$week_cont <- r2$week
r2$week <- factor(r2$week)
r2$year <- factor(r2$year)

# paired t tests
euks <- r2[r2$amplicon=="fungi", ]
euks
euks <- euks[euks$week!=0, ]
mean(euks[euks$type=="bulk",5])
# check normality
mean(euks$rotation_R2[euks$type=="bag"]-euks$rotation_R2[euks$type=="bulk"])
shapiro.test(euks$rotation_R2[euks$type=="bag"]-euks$rotation_R2[euks$type=="bulk"])
#ttest
eukttest <- t.test(rotation_R2~type, data = euks, paired = TRUE) # how to set up paired data frame
eukttest

bacs <- r2[r2$amplicon=="bacteria", ]
bacs <- bacs[bacs$week!=0, ]
# check normality
mean(bacs$rotation_R2[bacs$type=="bag"]-bacs$rotation_R2[bacs$type=="bulk"])
shapiro.test(bacs$rotation_R2[bacs$type=="bag"]-bacs$rotation_R2[bacs$type=="bulk"])
#ttest
bacttest <- t.test(rotation_R2~type, data=bacs, paired = TRUE)
bacttest

# make tables of results
ttesttable <- data.frame(matrix(nrow=2, ncol = 9))
colnames(ttesttable) <- c("amplicon", "t", "df", "p-value", "lower 95 CI", "upper 95 CI", "mean of differences", "mean bulk R2", "mean sentinel R2")
ttesttable[1,1] <- "Bacteria"
ttesttable[2,1] <- "Fungi"
ttesttable[1,2] <- bacttest$statistic #16 t
ttesttable[2,2] <- eukttest$statistic #its t
ttesttable[1,3] <- bacttest$parameter #df
ttesttable[2,3] <- eukttest$parameter #df
ttesttable[1,4] <- bacttest$p.value #16s pval
ttesttable[2,4] <- eukttest$p.value # its pval
ttesttable[1,5] <- bacttest$conf.int[1] # 16s lower CI
ttesttable[2,5] <- eukttest$conf.int[1] # its lower ci
ttesttable[1,6] <- bacttest$conf.int[2] # 16s upper ci
ttesttable[2,6] <- eukttest$conf.int[2] #its upper ci
ttesttable[1,7] <- bacttest$estimate # 16s mean diff
ttesttable[2,7] <- eukttest$estimate # its mean diff
ttesttable[1,8] <- mean(bacs[bacs$type=="bulk",5])
ttesttable[2,8] <- mean(euks[euks$type=="bulk",5])
ttesttable[1,9] <- mean(bacs[bacs$type=="bag",5])
ttesttable[2,9] <- mean(euks[euks$type=="bag",5])
ttesttable

ttesttable_sub <- ttesttable[ ,c(1,8,9,2,7,4)]

##### Paper Plots
r2$year_type <- paste0(r2$year, "_", r2$sample_type)
r2$sample_type <- factor(r2$type, labels = c("sentinel soil", "bulk soil"))
pointplot <- ggplot(r2[r2$week!=0, ])+
  theme(#legend.position = "none", 
    panel.grid.major.y = element_line(color = "gray85", size = 0.25),
    panel.grid.major.x = element_line(color = "gray85", size = 0.1),
    panel.ontop = FALSE,
    axis.text = element_text(size = 10), 
    axis.title = element_text(size = 12), 
    legend.text = element_text(size = 10), 
    strip.text = element_text(size = 12), 
    panel.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(fill = NA, color = "gray85"))+
  #geom_line(aes(x=week, y=r_squared, color = sample_type, group = year_type), size = 1)+
  geom_point(aes(x=week, y=r_squared, color = sample_type), size = 4)+
  scale_color_manual(values = c("#7DB0DD","#70582b"),name=NULL )+
  facet_grid(rows = vars(amplicon), cols = vars(year))+
  labs(y = bquote('Crop Rotation R'^2), x = "sampling week")+
  ylim(0.12,0.26)
pointplot