#LULU MVP Porites LSU amplicon sequencing data
#processing ASVs post DADA2
#Sept 2021

library(phyloseq); packageVersion("phyloseq") #version ‘1.36.0’
#install_github("tobiasgf/lulu") 
library(lulu); packageVersion("lulu") #version 0.1.0
library(dplyr); packageVersion("dplyr") #version 1.0.7

########Prepping the two inputs needed for LULU#######
  
#ASV fasta file input name: asvsforLULU_MVP_POR.fasta
seqtab <- as.matrix(read.csv("ASVtable_MVPPOR_LSU_0721.csv", sep=",", row.names=1)) #This used to be seqtab.nochim
path <- "POR_MVP_LSU.fasta"
uniquesToFasta(seqtab, path, ids = NULL, mode = "w", width = 20000)

#	Note: change header names: needs to be just sequence ID # (“>sq1”) and not (“>sq1; size=1321321;”)
#in terminal: cut -d ';' -f 1 POR_MVP_LSU.fasta > asvsforLULU_MVP_POR.fasta

#DADA2 output file name: ASVtableforLULU_MVPPOR_LSU_0721.csv

#######IN TERMINAL#########

#Produce a match list using BLASTn
#First produce a blastdatabase with the ASVs
#>makeblastdb -in asvsforLULU_MVP_POR.fasta -parse_seqids -dbtype nucl

#For two sequences to be merged, they must co-occur >90% of the time and have 84% sequence #similarity. These are the standard parameters
#Then blast the OTUs against the database to produce the match list 
#>blastn -db asvsforLULU_MVP_POR.fasta -outfmt '6 qseqid sseqid pident' -out match_list_MVP_POR.txt -qcov_hsp_perc 90 -perc_identity 84 -query asvsforLULU_MVP_POR.fasta

#######BACK TO R#########

alldat <-read.csv("ASVtableforLULU_MVPPOR_LSU_0721.csv") #add "SampleID" to header of column with sample IDs
head(alldat)

match_list<-read.table("match_list_MVP_POR.txt")
head(match_list)

#Reformat ASV table to desired LULU format
alldat$SampleID <- as.character(alldat$SampleID)
rownames(alldat)<-alldat$SampleID

ASVs<-data.frame(t(alldat[,2:336]))
head(ASVs)

#Now, run the LULU curation
curated_result <- lulu(ASVs, match_list, minimum_relative_cooccurence=0.9) #retains 4 ASVs at 84% similarity and 90% co-occurance
summary(curated_result)

curated_result$minimum_match
curated_result$minimum_relative_cooccurence

#write.csv(curated_result$otu_map,file="LULU_mergedoutput_summary_MVP_POR.csv",row.names=TRUE,quote=FALSE)

#########Now format and make into a phyloseq object#########
#make the results into a dataframe
alldat<-cbind(data.frame(curated_result$curated_table))
alldat.t <-data.frame(t(alldat))
head(alldat.t)

#make into a CSV, add a column name ("Sq.id"") to the "sq" columns for merging
write.csv(alldat.t,file="ASVtablepostLULU_MVP_POR.csv",row.names=TRUE,quote=FALSE)
otu <- as.matrix(read.csv("ASVtablepostLULU_MVP_POR.csv", sep=",", row.names=1)) 

#make taxonomy file with excel- vlookup to carry over original taxonomic IDS from blasting pre-LULU ASVs
taxa<-as.matrix(read.csv("taxa_LULU.csv", row.names = 1))
head(taxa)

#read in meta data file
meta<-read.csv("LSU_POR_MVP_meta.csv", sep=",", row.names=1)

# Construct phyloseq object
psLSU <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                  sample_data(meta), 
                  tax_table(taxa))
psLSU

saveRDS(psLSU, file="LULU_PORpsMVP_LSU_0921.rds")


