
library(reshape2)
library(tidyverse)
library(tidyr)
library(ggbeeswarm)
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("phyloseq")
library(DESeq2); packageVersion("DESeq2") #version ‘1.28.1’
library(vegan); packageVersion("vegan") #version ‘2.5-7’
library(patchwork)

setwd("/Users/laurenhowe-kerr/Dropbox/LaurenMVP/POR/Rscripts")

#read in filtered and rlog transformed phyloseq object
ps <- readRDS("PORdinoRNAVMVP_transformed.rds")
ps

bc_dist <- phyloseq::distance(ps, method = "bray")
sampledf <- data.frame(sample_data(ps))

#########################################################FIGURE 5########################################################################################################################


#############################################
######Between Group Distance Graphs##########
#############################################

#by reef zone
colnames(sampledf)
trt_labels <- sampledf[,c(2,11)]
distances <- bc_dist %>%
  as.matrix() %>%
  melt() %>%
  merge(trt_labels, by.x = "Var1", by.y = "SampleID.1") %>%
  merge(trt_labels, by.x = "Var2", by.y = "SampleID.1") %>%
  filter(Reef.Type.x != Reef.Type.y)

# remove dublicate comparisons
tempa <- distances[which(distances$Reef.Type.x=="FOR" & 
                           distances$Reef.Type.y=="FRG"),]
tempb <- distances[which(distances$Reef.Type.x=="FOR" & 
                           distances$Reef.Type.y=="BAK"),]
tempc <- distances[which(distances$Reef.Type.x=="FRG" & 
                           distances$Reef.Type.y=="BAK"),]
tempr <- rbind(tempa,tempb, tempc)

# add new column
tempr$comp <- paste(tempr$Reef.Type.x, tempr$Reef.Type.y, sep=" v ")
tempr$comp<- factor(tempr$comp, levels = c("FRG v BAK", "FOR v BAK", "FOR v FRG"), labels=c("Fringe v Back", "Fore v Back", "Fore v Fringe"))


bc_group_reef <- ggplot(data = tempr, aes(x = comp, y = value)) +
  theme_classic(base_size=12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_point(position = position_jitter(width = .42, height = 0), shape = 1, colour = "grey") +
  geom_boxplot(alpha=0, outlier.alpha =0) +
  xlab("Reef type comparison") + ylab("Bray Curtis dissimilarity")
bc_group_reef


#############################################
##Within Group Distance Graphs (Dispersion)##
#############################################

#############################################
#reef

bc_dist <- phyloseq::distance(ps, method = "bray")
sampledf <- data.frame(sample_data(ps))

#Adding distance to centroid to mapping dataset
# Run betadisper on the treatments (or whichever groups you are comparing)
betabc <- betadisper(bc_dist, s$Reef.Type, bias.adjust = TRUE)

# Extract within group distances, gives the distances of each group to its centroid by interaction
withdistbc <- betabc$distances

SampleID <- row.names(sampledf) # Extracts the SampleIDs from the dataframe
raw_distances <- data.frame(SampleID, withdistbc) # makes a new dataframe with the specified columnes
s <- data.frame(sample_data(ps)) %>% unite("Reef.Time", Reef.Type:Time2, sep="-", remove=FALSE)
# Makes metadata into a df to work with
colnames(s)[2] <- "SampleID" # Change first column title to SampleID to match distances dataframe
with_distances <- merge(raw_distances, s, by = "SampleID") # merges metadata df and distances df
colnames(with_distances)[2] <- c("withdistbc")
#withdistmax <- read.csv("betadistances_POR_rlog.csv")

with_distances$Reef.Type <- factor(with_distances$Reef.Type, levels = c("FRG", "BAK", "FOR"), labels=c("Fringe", "Back", "Fore"))


bc_with_reef <- ggplot(data = with_distances, aes(x=Reef.Type, y=withdistbc)) +
  geom_boxplot(alpha=0.5) + 
  xlab("Reef type") + ylab("Bray Curtis within group distance") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #stat_summary(fun.y=mean, geom="point", shape=20, size=9, color="red", fill="red") +
  geom_beeswarm(alpha=1, size=2, aes(color=Reef.Type)) +
  scale_color_manual(values=c( "#4F93B8", "#AD8CAE", "#222B4C")) + ylim(0,0.25) 
bc_with_reef

#############################################
#time

bc_dist <- phyloseq::distance(ps, method = "bray")
sampledf <- data.frame(sample_data(ps))

#Adding distance to centroid to mapping dataset
# Run betadisper on the treatments (or whichever groups you are comparing)
betabc <- betadisper(bc_dist, s$Time2, bias.adjust = TRUE)

# Extract within group distances, gives the distances of each group to its centroid by interaction
withdistbc <- betabc$distances

SampleID <- row.names(sampledf) # Extracts the SampleIDs from the dataframe
raw_distances <- data.frame(SampleID, withdistbc) # makes a new dataframe with the specified columnes
s <- data.frame(sample_data(ps)) %>% unite("Reef.Time", Reef.Type:Time2, sep="-", remove=FALSE)
# Makes metadata into a df to work with
colnames(s)[2] <- "SampleID" # Change first column title to SampleID to match distances dataframe
with_distances <- merge(raw_distances, s, by = "SampleID") # merges metadata df and distances df
colnames(with_distances)[2] <- c("withdistbc")

bc_with_time <- ggplot(data = with_distances, aes(x=Time2, y=withdistbc)) +
  geom_boxplot(alpha=0.5, outlier.alpha = 0) + 
  xlab("Sampling time") + ylab("Bray Curtis within group distance") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_beeswarm(alpha=1, size=2, aes(color=Reef.Type)) +
  scale_color_manual(values=c( "#4F93B8", "#AD8CAE", "#222B4C")) + ylim(0,0.25)
bc_with_time


bc_group_reef + bc_with_time + plot_annotation(tag_levels = 'A')
#ggsave("Figure-5-feb23.jpeg", width = 6, height = 4.5)

bc_with_reef
#ggsave("Supplemental_fig_6.jpeg", width = 5.5, height = 4.5)

#############################################
##################Stats######################
#############################################
#r log transformed ps object
ps <- readRDS("PORdinoRNAVMVP_transformed.rds")
s <- data.frame(sample_data(ps))
ps_bc <- phyloseq::distance(ps, method = "bray")

adonis(ps_bc ~ Reef.Type, data = s) #Reef.Type   2   0.02388 0.0119419  2.0254 0.02647  0.007 **
adonis(ps_bc ~ Final.Health, data = s) #Final.Health   1   0.01312 0.0131204  2.2131 0.01454  0.019 *
adonis(ps_bc ~ Time2, data = s) #Time2       3   0.00958 0.0031919 0.46884 0.00811      1
pairwise.adonis.dm <- function(x,factors,stratum=NULL,p.adjust.m="bonferroni",perm=999){
  
  library(vegan)
  
  if(class(x) != "dist") stop("x must be a dissimilarity matrix (dist object)")
  
  co = as.matrix(combn(unique(factors),2))
  
  pairs = c()
  
  F.Model =c()
  
  R2 = c()
  
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    
    sub_inds <- factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))
    
    resp <- as.matrix(x)[sub_inds,sub_inds]
    
    ad = adonis(as.dist(resp) ~
                  
                  factors[sub_inds], strata=stratum[sub_inds], permutations=perm);
    
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    
    R2 = c(R2,ad$aov.tab[1,5]);
    
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
}


pairwise.adonis.dm(ps_bc, s$Reef.Type, p.adjust.m = "BH")
#pairs  F.Model         R2 p.value p.adjusted
#1 FOR vs BAK 1.689494 0.01473103   0.085      0.085
#2 FOR vs FRG 1.727205 0.01559874   0.076      0.085
#3 BAK vs FRG 2.294851 0.01817058   0.013      0.039

beta <- betadisper(ps_bc, s$Reef.Type) #Groups     2 0.01331 0.0066556 1.6339    999  0.199
permutest(beta)
beta <- betadisper(ps_bc, s$Time2) #Groups     3 0.02955 0.0098508 2.4333    999  0.058 .
permutest(beta)

#########################################################FIGURE 6########################################################################################################################

#############################################
##Within Group Distance Graphs (Dispersion)##
#############################################
ps <- readRDS("PORdinoRNAVMVP_transformed.rds")

ps <- subset_samples(ps, Final.Health != "skip") #152 samples, analyzing samples with a non-ambiguous health trajectory

bc_dist <- phyloseq::distance(ps, method = "bray")
s <- data.frame(sample_data(ps))

#Adding distance to centroid to mapping dataset
# Run betadisper on the treatments (or whichever groups you are comparing)
betabc <- betadisper(bc_dist, s$StressYear.Health, bias.adjust = TRUE)

# Extract within group distances, gives the distances of each group to its centroid by interaction
withdistbc <- betabc$distances

SampleID <- row.names(s) # Extracts the SampleIDs from the dataframe
raw_distances <- data.frame(SampleID, withdistbc) # makes a new dataframe with the specified columnes
#s <- data.frame(sample_data(ps)) %>% unite("Reef.Time", Reef.Type:Time2, sep="-", remove=FALSE)
# Makes metadata into a df to work with
colnames(s)[2] <- "SampleID" # Change first column title to SampleID to match distances dataframe
with_distances <- merge(raw_distances, s, by = "SampleID") # merges metadata df and distances df
colnames(with_distances)[2] <- c("withdistbc")
#withdistmax <- read.csv("betadistances_POR_rlog.csv")

with_distances$StressYear.Health <- factor(with_distances$StressYear.Health , levels=c("normal-no", "thermal.stress-no", "normal-yes", "thermal.stress-yes"))

bc_with_HT <- ggplot(data = with_distances, aes(x=StressYear.Health, y=withdistbc, color=Stress.Year)) +
  geom_boxplot(alpha=0.5) +
  xlab("") + ylab("Bray Curtis within group distance") +
  theme_classic(base_size=12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position="none") + 
  scale_color_manual(values=c("black", "#D45769"))+
  geom_beeswarm(alpha=1, size=2)
bc_with_HT 

bc_with_HT
#ggsave("Figure-6-feb23.jpeg", width = 4, height = 4.5)

#stats
ps <- readRDS("PORdinoRNAVMVP_transformed.rds")
s <- data.frame(sample_data(ps)) %>% unite("Reef.Health", Reef.Type, Final.Health, sep="-", remove=FALSE) %>% unite("Time.Health", Time2, Final.Health, sep="-", remove=FALSE) %>% unite("StressYear.Health", Stress.Year, Final.Health, sep="-", remove=FALSE) %>% unite("Reef.Time", Reef.Type, Time2, sep="-", remove=FALSE)
sample_data(ps) <- s

ps <- subset_samples(ps, Final.Health != "skip") #152 samples

bc_dist <- phyloseq::distance(ps, method = "bray")
s <- data.frame(sample_data(ps))

#####

beta <- betadisper(bc_dist, s$Final.Health) #Groups      1 0.00051 0.0005056 0.1461    999  0.707
permutest(beta)
beta <- betadisper(bc_dist, s$Stress.Year) #Groups      1 0.01733 0.017333 5.1355    999  0.033 *
permutest(beta)

adonis(bc_dist ~ Stress.Year, data = s) #1   0.00404 0.0040366 0.71733 0.00476  0.809
adonis(bc_dist ~ Final.Health, data = s) #Final.Health   1   0.01107 0.0110659   1.983 0.01305   0.03 *

adonis(bc_dist ~ Final.Health*Stress.Year, data = s) #Final.Health:Stress.Year   1   0.00258 0.0025770 0.45936 0.00304  0.987

####
ps <- readRDS("PORdinoRNAVMVP_transformed.rds")
s <- data.frame(sample_data(ps)) %>% unite("Reef.Health", Reef.Type, Final.Health, sep="-", remove=FALSE) %>% unite("Time.Health", Time2, Final.Health, sep="-", remove=FALSE) %>% unite("StressYear.Health", Stress.Year, Final.Health, sep="-", remove=FALSE) %>% unite("Reef.Time", Reef.Type, Time2, sep="-", remove=FALSE)
sample_data(ps) <- s

#test both
ps <- subset_samples(ps, Final.Health == "yes") #84 samples
ps <- subset_samples(ps, Final.Health == "no") #68 samples

bc_dist <- phyloseq::distance(ps, method = "bray")
s <- data.frame(sample_data(ps))

adonis(bc_dist ~ Stress.Year, data = s) #1   0.00404 0.0040366 0.71733 0.00476  0.809
beta <- betadisper(bc_dist, s$Stress.Year) 
permutest(beta)

#for yes decline in health
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.008671 0.0086714 3.4786    999  0.057 .
#Residuals 82 0.204411 0.0024928 

#for no decline in health
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.008728 0.0087279 1.9775    999  0.169
#Residuals 66 0.291295 0.0044136 

