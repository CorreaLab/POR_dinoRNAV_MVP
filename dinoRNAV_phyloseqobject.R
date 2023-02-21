library(tidyr)
library(ggbeeswarm)
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("phyloseq")
library(ggpubr)
library(DESeq2); packageVersion("DESeq2") #version ‘1.28.1’
library(RColorBrewer)
library(vegan); packageVersion("vegan")
library(dunn.test)
library(scales)
library(lme4)
library(nlme)
packageVersion("stats")
library(patchwork)
library(plyr)
library(yarrr)

# Set working directory and reef in necessary files 
setwd("/Users/laurenhowe-kerr/Dropbox/LaurenMVP/POR/Rscripts")
taxa <- as.matrix(read.csv("PorMVP_AminoTypes1.csv", row.names = 1))
samdf <- read.csv("POR_dinoRNAV_sampledata1.csv", row.names = 1)
counts <- read.csv("PorMVP_protcounts1.csv", row.names=1, header = TRUE, check.names=FALSE)

# Creat Phyloseq object 
ps <- phyloseq(otu_table(counts, taxa_are_rows=FALSE),sample_data(samdf), tax_table(taxa))
ps

# Save this object to your folder 
#saveRDS(ps, file ="PORdinoRNAVMVP.rds")
ps <- readRDS("PORdinoRNAVMVP.rds")


# Filtering 

ps2 <- prune_samples(sample_sums(ps)>=1000, ps) #remove samples with <1000 reads 
ps3 <- prune_taxa(taxa_sums(ps2)>0, ps2)
#rename back to "ps" for simplicity
ps <-ps3

#add labels, levels
sample_data(ps)$Reef.Type <- factor(sample_data(ps)$Reef.Type, levels = c("FOR", "BAK", "FRG"))
sample_data(ps)$Time2 <- factor(sample_data(ps)$Time2, levels = c("818", "319", "819", "1020"), labels=c("Aug18", "Mar19", "Aug19", "Oct20"))


#rlog transformation, to account for differences in read depth for beta diversity analyses
##first transform counts data using DESeq2
ps.NT <- ps ##this makes a copy of the phyloseq object so that you can still easily reference the NOT TRANSFORMED (ps.NT) version of your ps object
diagdds = phyloseq_to_deseq2(ps, ~ CT) ##chose the main variable of distinction but this is a not super important 

##rlog transform (make sure all dependencies and packages have loaded properly)
rlogCOUNTS<-rlog(diagdds,blind=TRUE)
head(assay(rlogCOUNTS))

dat=as.data.frame(assay(rlogCOUNTS)) 

dat[dat < 0.0] <- 0.0 ## remove negative values before calculating bray curtis distances
otu_table(ps) <- otu_table(dat, taxa_are_rows = TRUE) ## put this transformed data back into phyloseq object as ASV table
XsaveRDS(ps, file ="PORdinoRNAVMVP_transformed.rds")

#this is the ps file used for beta div analyses