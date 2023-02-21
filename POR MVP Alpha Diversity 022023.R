
# Load necessary packages 
library(ggbeeswarm)
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(ggpubr)
library(DESeq2); packageVersion("DESeq2") #version ‘1.28.1’
library(RColorBrewer)
library(vegan)
library(dunn.test)
library(scales)
library(lme4)
library(nlme)

# Set your working directory and reef in necessary files 
setwd("/Users/laurenhowe-kerr/Dropbox/LaurenMVP/POR/Rscripts")
ps <- readRDS("PORdinoRNAVMVP.rds")

# Filtering 

ps2 <- prune_samples(sample_sums(ps)>=1000, ps) #remove samples with <1000 reads 
ps3 <- prune_taxa(taxa_sums(ps2)>0, ps2)
#rename back to "ps" for simplicity
ps <-ps3

#add labels, levels
sample_data(ps)$Reef.Type <- factor(sample_data(ps)$Reef.Type, levels = c("FOR", "BAK", "FRG"))
sample_data(ps)$Time2 <- factor(sample_data(ps)$Time2, levels = c("818", "319", "819", "1020"), labels=c("Aug18", "Mar19", "Aug19", "Oct20"))


# rarefy without replacement
set.seed(1)
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)

ps <- ps.rarefied

#### Alpha Diversity Statistics 
erich <- estimate_richness(ps, measures=c("Observed", "Shannon", "Simpson"))
#make dataframe
SampleID <- row.names(erich) 
# makes a new dataframe with the specified columnes
alpha <- data.frame(SampleID, erich) 
s <- data.frame(sample_data(ps)) 
# Change appropriate column title to SampleID to make sure dataframes match 
colnames(s)[2] <- "SampleID" 
colnames(alpha)[1] <- "SampleID" 
# merges metadata df and distances df 
alphadiv <- merge(s,alpha, by = "SampleID")

alphadiv$CT <- factor(alphadiv$CT)
#Reording Time
alphadiv$Time2 <- factor(alphadiv$Time2 , levels=c("818", "319", "819", "1020"))
alphadiv$Reef.Type <- factor(alphadiv$Reef.Type , levels=c("FRG", "BAK", "FOR"))

#test for normal distribution
shapiro.test(alphadiv$Shannon) #not normal (p<0.05)
shapiro.test(alphadiv$Simpson) #not normal (p<0.05)

########################################################################################################################

#assessing differences in diversity according to reef zone & time

#violin
alpha_ObservedTime <- ggviolin(data = alphadiv, x="Reef.Type", y="Observed", add = "boxplot", add.params = list(fill = "white")) +
  xlab("") + ylab("Richness (# aminotypes)") +
  #facet_grid(~Reef.Type) +
  theme(panel.background = element_blank(),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title=element_blank(),
        axis.text = element_text(size =12),
        #axis.text.x = element_blank(),
        #axis.line = element_blank(),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position="none")
 #+geom_beeswarm(alpha=1, size=1)
alpha_ObservedTime 

frg <- filter(alphadiv, Reef.Type=="FRG")
bak <- filter(alphadiv, Reef.Type=="BAK")
forr <- filter(alphadiv, Reef.Type=="FOR")

aug18 <- filter(alphadiv, Time2=="818")
aug19 <- filter(alphadiv, Time2=="819")
mar19 <- filter(alphadiv, Time2=="319")
oct20 <- filter(alphadiv, Time2=="1020")


#test for normal distribution
shapiro.test(alphadiv$Observed) #not normal (p<0.05)
shapiro.test(frg$Observed)
shapiro.test(forr$Observed)
shapiro.test(bak$Observed)

shapiro.test(aug18$Observed)
shapiro.test(aug19$Observed)
shapiro.test(mar19$Observed)
shapiro.test(oct20$Observed)

#non parametric test
kruskal.test(Observed ~ Reef.Type, data=alphadiv)
#data:  Observed by Reef.Type
#Kruskal-Wallis chi-squared = 16.126, df = 2, p-value = 0.000315

kruskal.test(Observed ~ Time, data=alphadiv)
#data:  Observed by Time2
#Kruskal-Wallis chi-squared = 0.92717, df = 3, p-value = 0.8189

pairwise.wilcox.test(alphadiv$Observed, alphadiv$Reef.Type, p.adj = "bonf")

#but each timepoint is normally distributed, so using a linear model:
diversity_mixed = lm(Observed ~ Reef.Type * Time, data = alphadiv)
summary(diversity_mixed)
anova(diversity_mixed)

shapiro.test(resid(diversity_mixed))
hist(resid(diversity_mixed))
qqnorm(resid(diversity_mixed))

#mean values of richness, by different categories 
tapply(alphadiv$Observed, alphadiv$Reef.Type, mean)
tapply(alphadiv$Observed, alphadiv$Reef.Type, median)

####################################################SUPP FIGURE############################################################
#plotting richness vs zoox density
zoox = lm(Observed ~ ZooxDensity, data = alphadiv)
summary(zoox)
anova(diversity_mixed)
qqnorm(resid(zoox))
shapiro.test(resid(zoox))  #normal

zd <- ggplot(alphadiv, aes(x = ZooxDensity, y = Observed)) + 
  geom_point() +
  xlab("Symbiodiniaceae Density (cells per sq cm)") + ylab("Richness (# aminotypes)") +
  theme_classic() +
  stat_smooth(method = "lm", col="black")

#ggsave("richness_zooxdensity.jpeg", width = 4, height = 4)


