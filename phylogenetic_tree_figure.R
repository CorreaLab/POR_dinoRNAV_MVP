### Phylogenetic tree figure
### POR dinoRNAVs

library("treeio")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(ggtree); packageVersion("ggtree")
library(ggtreeExtra) ; packageVersion("ggtreeExtra")
library(ggnewscale)
library(ggstance)
library(tidytree)
library(ape)
library(tidyverse)
library(ggstar)
library(RColorBrewer)
library(nord) #for colors
library(yarrr)

setwd("/Users/laurenhowe-kerr/Dropbox/LaurenMVP/POR/phylo_tree")

#this new file has only POR MVP aminotypes present at at least 10% rel abundance in at least 1 sample, also contains PVID seqs, HcRNAV, Levin, etc 
tree <- read.iqtree("/Users/laurenhowe-kerr/Dropbox/LaurenMVP/POR/phylo_tree/alignmentfiles_treefigure/PorMVP_finalseqset_.treefile")

#make into a phyloseq format
phylo <- as.phylo(tree)

#collapsing branches poorly supported
Badnodes <- which(as.numeric(phylo$node.label) < 50) + length(phylo$tip.label)
Badnodes_indexes <- c()
for(node in Badnodes){
  Badnodes_indexes <- c(Badnodes_indexes, which(phylo$edge[,2] == node))
}

phylo$edge.length[Badnodes_indexes] <- 0 # Make the branch length of the bad nodes = 0
tree_multi <- di2multi(phylo) #Use di2multi function to convert the branches with lenght = 0, to multichotomies


# Read in other components to build phyloseq object, set your working directory 
taxa <- as.matrix(read.csv("PorMVP_AminoTypes_taxa_treeextras.csv", row.names = 1)) #added the "extras" (non POR samples) to this file
samdf <- read.csv("/Users/laurenhowe-kerr/Dropbox/LaurenMVP/POR/Rscripts/POR_dinoRNAV_sampledata1.csv", row.names = 1)
counts <- read.csv("countstable_fortree_0921.csv", row.names=1, header = TRUE, check.names=FALSE) #also added the extras here, with 0 abundance in all samples (to trick phyloseq into keeping these other taxa that I want to plot in the tree)

########################################################
#to make file "countstable_fortree_0921.csv", do the following:
ps2 <- prune_samples(sample_sums(ps)>=1000, ps) #if there are obvious low outliers, remove (we're removing any sequence with less than 1000 reads)
ps3 <- prune_taxa(taxa_sums(ps2)>0, ps2)
ps4 <- subset_samples(ps3, CT != "water" & CT != "neg")
ps5 <- prune_taxa(taxa_sums(ps4)>0, ps4)
ps6 <- subset_samples(ps5, SampleID.1 != "POR25_1020_Plate2" & SampleID.1 != "POR8_1020_Plate2" & SampleID.1 != "POR146_1020_Plate2")
ps7 <- prune_taxa(taxa_sums(ps6)>0, ps6) #the 3 duplicate  samples sequenced on both plates had no unique amino types on plate 2
#water sample had no unique ASVS
#rename back to "ps" for simplicity
ps <- ps7 #should have 778 taxa and 178 samples after all of this filtering
ps  = transform_sample_counts(ps, function(x) x / sum(x) )
#find list of aminotypes that are present at at least 5% rel abundance in at least 1 sample:
psmelted <- psmelt(ps)
psmelted
#find list of aminotypes to keep
psmelted$AminoType[psmelted$Abundance < 0.1] <- "< 10% abund."
psmelted$AminoType <- factor(psmelted$AminoType)
keep <- as.character(unique(psmelted$AminoType))
write.csv(keep, file="aminotypestokeep.csv")
#prune un-melted ps to just include your "keep" aminotypes
ps = prune_taxa(keep, ps)
tab <- otu_table(ps)
write.csv(tab, file="countstable_fortree_0921.csv")
########################################################

ps <- phyloseq(otu_table(counts, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa),
               phy_tree(tree_multi))
ps

#saveRDS(ps, "ps-treefig-MVPPOR.rds")
ps<- readRDS("ps-treefig-MVPPOR.rds")

#order
sample_data(ps)$Time2 <- factor(sample_data(ps)$Time2 , levels=c("818", "319", "819", "1020"))
sample_data(ps)$Reef.Type <- factor(sample_data(ps)$Reef.Type , levels=c("FOR", "BAK", "FRG"))

#count # of colonies per reef zone with the aminotype in abundance > 0.1
melt <- psmelt(ps) %>% 
  filter(Abundance >0.1) %>% 
  group_by(Reef.Type,CT,OTU) %>%
  count() %>% 
  ungroup() %>% 
  group_by(Reef.Type, OTU) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(reef=Reef.Type, ID=OTU) %>% 
  select(ID, reef, n)

melt$reef <- factor(melt$reef, levels=c("FRG", "BAK", "FOR"))


#count # of colonies of each health category with the aminotype in abundance > 0.1
health.count <- psmelt(ps) %>% 
  filter(Final.Health != "skip") %>% 
  filter(Abundance >0.1) %>% 
  group_by(Final.Health, CT_2,OTU) %>%
  count() %>% 
  ungroup() %>% 
  group_by(Final.Health, OTU) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(tissue_loss=Final.Health, ID=OTU) %>% 
  select(ID, tissue_loss, n) 

#count # of colonies of each health category IN EACH REEF ZONE with the aminotype in abundance > 0.1
reef.health.count <- psmelt(ps) %>% 
  filter(Final.Health != "skip") %>% 
  filter(Abundance >0.1) %>% 
  unite("reef.health", Reef.Type, Final.Health, sep="-", remove=TRUE) %>% 
  group_by(reef.health, CT_2,OTU) %>%
  count() %>% 
  ungroup() %>% 
  group_by(reef.health, OTU) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(reefhealth=reef.health, ID=OTU) %>% 
  select(ID, reefhealth, n) 

#need to figure out how to count number of unique CTs per reef health category and divide n/# to get % of colonies that have the aminotype
reef.health.count
reef.health.count$reefhealth <- factor(reef.health.count$reefhealth, levels=c("FRG-no", "FRG-yes", "BAK-no", "BAK-yes", "FOR-no", "FOR-yes"))

#get sequence source info
taxa2 <- read.csv("PorMVP_AminoTypes_taxa_treeextras.csv")
origin <- taxa2 %>% 
  select(X, Study, Source) %>% 
  mutate(source2=Source, study2=Study) %>% 
  select(X, source2, study2)


piratepal(palette = "all")
piratepal(palette = "eternal")

#make basic tree, label nodes with >70 bootstrap support
p <- ggtree(ps, ladderize = TRUE, size=1, color="darkgray") +
  geom_nodepoint(aes(label=label, subset = as.numeric(label) >70), size=2) + geom_treescale(fontsize=4, linesize=1, offset=1)

#add sequence source
r <- p + 
  geom_fruit(data=origin, geom=geom_tile, mapping=aes(y=X, fill=study2, width=0.2), pwidth=0.2, offset = 0.03) + 
  scale_fill_manual(values=c("#473B75FF", "#52194CFF","#AA4371", "#4D709CFF", "#6F766BFF", "#92ADC4FF", "lightgray"), guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  )  

#better colors Jan 2022
r <- p + 
  geom_fruit(data=origin, geom=geom_tile, mapping=aes(y=X, fill=study2, width=0.2), pwidth=0.2, offset = 0.03) + 
  scale_fill_manual(values=c("#473B75FF", "#AA4371", "#4D709CFF", "black", "#6F766BFF", "#92ADC4FF", "lavender"), guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  )  


#add heat map of abundance by reef zone
t <- r + new_scale_fill() +
  geom_fruit(data=melt, geom=geom_tile,
             mapping=aes(y=ID, x=reef, fill=n),
             offset = 0.05)+
  scale_fill_gradientn(colours = c("white", "black", "black"), breaks = c(1,2,5,15), limits=c(0,15), guide = "legend")

#add heat map of abundance by health state
s <- t + new_scale_fill() +
  geom_fruit(data=health.count, geom=geom_tile,
             mapping=aes(y=ID, x=tissue_loss, fill=n),
             offset = 0.08, pwidth=0.1, colnames=TRUE)+
  scale_fill_gradientn(colours = c("white", "black", "black"), breaks = c(1,2,5,15), limits=c(0,15), guide = "legend")


ggsave("fig-3-0722-V4-reefzone.jpeg", width = 12, height = 10)

#add heat map of abundance by health state
u <- r + new_scale_fill() +
  geom_fruit(data=reef.health.count, geom=geom_tile,
             mapping=aes(y=ID, x=reefhealth, fill=n),
             offset = 0.03, pwidth=0.2)+
  scale_fill_gradientn(colours = c("white", "black", "black"), breaks = c(1,2,5,15), limits=c(0,15), guide = "legend")

ggsave("fig-3-0722-V2.jpeg", width = 12, height = 10)