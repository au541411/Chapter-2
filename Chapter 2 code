# Chapter-2
# Scripts and code for Chapter 2
# To be executed in R
library(ggforce)
library(tidyverse)
library(limma)
library(plyr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(stackedbarplot)
library(reshape2)
library(ampvis2)
library(wesanderson)
# data from amplicon sequencing in home directory
#dataAD 
#dataKS
# trim primer sequences
cutadapt_path_PE("dataAD","GCSTACSYSATSTACACSTCSGG","SASGTCVCCSGTSCGGTA") cutadapt_path_PE("dataKS","GCNATGGAYCCNCARCARMGNVT","GTNCCNGTNCCRTGNSCYTCNAC")
# variable 'PATH' to directory 
pathAD <- "~/dataAD" 
pathKS <- "~/dataKS"
# trimmed fastq files in list of forward and reverse read files 
# make a list of sample names
fnFsAD <- sort(list.files(pathAD, pattern="_R1_001.trimmed.fastq.gz", full.names = TRUE)) fnRsAD <- sort(list.files(pathAD, pattern="_R2_001.trimmed.fastq.gz", full.names = TRUE)) sample.namesAD <- sapply(strsplit(basename(fnFsAD), "_"), `[`, 1) fnFsKS <- sort(list.files(pathKS, pattern="_R1_001.trimmed.fastq.gz", full.names = TRUE)) fnRsKS <- sort(list.files(pathKS, pattern="_R2_001.trimmed.fastq.gz", full.names = TRUE)) sample.namesKS <- sapply(strsplit(basename(fnFsKS), "_"), `[`, 1)
# plot read quality for the forward reads 
plotQualityProfile(fnFsAD[1:15]) 
plotQualityProfile(fnFsKS[1:15])
# filter and trim reads 
# set file name and location for the filtered/trimmed files
filtFsAD <- file.path(pathAD, "filtered", paste0(sample.namesAD, "_F_filt.fastq.gz")) 
filtRsAD <- file.path(pathAD, "filtered", paste0(sample.namesAD, "_R_filt.fastq.gz"))
filtFsKS <- file.path(pathKS, "filtered", paste0(sample.namesKS, "_F_filt.fastq.gz")) 
filtRsKS <- file.path(pathKS, "filtered", paste0(sample.namesKS, "_R_filt.fastq.gz"))
# trim and filter the reads
outAD <- filterAndTrim(fnFsAD, filtFsAD, fnRsAD, filtRsAD, truncLen=c(240,150),maxN=0, maxEE=c(2,2), truncQ=0, rm.phix=TRUE, compress=TRUE, multithread=4) 
outKS <- filterAndTrim(fnFsKS, filtFsKS, fnRsKS, filtRsKS, truncLen=c(240,150), maxN=0, maxEE=c(2,2), truncQ=0, rm.phix=TRUE, compress=TRUE, multithread=4)
# build error models from your data
errFAD <- learnErrors(filtFsAD, multithread=3) 
errRAD <- learnErrors(filtRsAD, multithread=3) 
errFKS <- learnErrors(filtFsKS, multithread=3)
errRKS<- learnErrors(filtRsKS, multithread=3)
# view error models 
plotErrors(errFAD, nominalQ=TRUE)
plotErrors(errRAD, nominalQ=TRUE) 
plotErrors(errFKS, nominalQ=TRUE)
plotErrors(errRKS, nominalQ=TRUE)
# dereplicate reads and name derep-class object by sample names
derepFsAD <- derepFastq(filtFsAD, verbose=TRUE)
derepRsAD <- derepFastq(filtRsAD, verbose=TRUE) 
names(derepFsAD) <- sample.namesAD 
names(derepRsAD) <- sample.namesAD 
derepFsKS <- derepFastq(filtFsKS, verbose=TRUE)
derepRsKS <- derepFastq(filtRsKS, verbose=TRUE)
names(derepFsKS) <- sample.namesKS 
names(derepRsKS) <- sample.namesKS
# cluster reads into amplicon sequence variants (ASVs) 
dadaFsAD <- dada(derepFsAD, err=errFAD, multithread=4) 
dadaRsAD <- dada(derepRsAD, err=errRAD, multithread=4) 
dadaFsKS <- dada(derepFsKS, err=errFKS, multithread=4) 
dadaRsKS <- dada(derepRsKS, err=errRKS, multithread=4)
# concatenate forward and reverse reads 
mergersAD <- mergePairs(dadaFsAD, derepFsAD, dadaRsAD, derepRsAD, verbose=TRUE, justConcatenate = TRUE) 
mergersKS <- mergePairs(dadaFsKS, derepFsKS, dadaRsKS, derepRsKS, verbose=TRUE, justConcatenate = TRUE)
# sequence table from merged sequences 
seqtabAD <- makeSequenceTable(mergersAD) 
seqtabKS <- makeSequenceTable(mergersKS)
# distribution of sequence lengths 
table(nchar(getSequences(seqtabAD))) 
table(nchar(getSequences(seqtabKS)))
# remove chimeras from dataset 
seqtab.nochimAD <- removeBimeraDenovo(seqtabAD, method="consensus", multithread=3, verbose=TRUE)
seqtab.nochimKS <- removeBimeraDenovo(seqtabKS, method="consensus", multithread=4, verbose=TRUE)
# overview table of the number of reads in each sample at each stage of the process 
getN <- function(x) sum(getUniques(x))
trackAD <- cbind(outAD, sapply(dadaFsAD, getN), sapply(dadaRsAD, getN), sapply(mergersAD, getN), rowSums(seqtab.nochimAD))
colnames(trackAD) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(trackAD) <- sample.namesAD
getN <- function(x) sum(getUniques(x))
trackKS <- cbind(outKS, sapply(dadaFsKS, getN), sapply(dadaRsKS, getN), sapply(mergersKS, getN), rowSums(seqtab.nochimKS))
colnames(trackKS) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(trackKS) <- sample.namesKS
# If downsampling is needed use the following, after trimming
# down-sampling after trimming 
downsample_PE <- function(forward_reads, reverse_reads, number)
{ 
f <- FastqSampler(forward_reads, n=number) 
r <- FastqSampler(reverse_reads, n=number) 
set.seed(123L); foutput <- yield(f) 
set.seed(123L); routput <- yield(r) 
forward_output_filename <-paste(strsplit(forward_reads,"[.]")[[1]][1],"downsampled.fastq.gz",sep='.') 
reverse_output_filename <- paste(strsplit(reverse_reads,"[.]")[[1]][1],"downsampled.fastq.gz",sep='.')
writeFastq(foutput, forward_output_filename, mode="w", compress=TRUE) 
writeFastq(routput, reverse_output_filename, mode="w", compress=TRUE)
}
downsample_path_PE <- function(path, n)
{
forward_reads_files <- sort(list.files(path, pattern="_F_filt.fastq.gz", full.names = TRUE))
reverse_reads_files <- sort(list.files(path, pattern="_R_filt.fastq.gz", full.names = TRUE))
lapply(1:length(forward_reads_files), function(x)
downsample_PE(forward_reads_files[x],reverse_reads_files[x],n))
# save as down-sampling.R
# in MD file
source("down-sampling.R")
downsample_path_PE("dataAD/filtered",26661)
downsample_path_PE("dataKS/filtered",4504)
# set location for the down-sampled files
filtFsADds <- file.path(pathAD, "filtered", paste0(sample.namesAD, "_F_filt.downsampled.fastq.gz"))
filtRsADds <- file.path(pathAD, "filtered", paste0(sample.namesAD, "_R_filt.downsampled.fastq.gz"))
filtFsKSds <- file.path(pathKS, "filtered", paste0(sample.namesKS, "_F_filt.downsampled.fastq.gz"))
filtRsKSds <- file.path(pathKS, "filtered", paste0(sample.namesKS, "_R_filt.downsampled.fastq.gz"))
# repeat protocol from above from step ‘built error models’ 
# files labeled …ds

# This work is done on the server in the terminal
# 6 frame translate ASV nucleotide sequences into amino acid sequences
# filter them with antismash HMM's and then align them with DIAMOND blast against MiBig database.
# 6 frame translation #
# make a folder and place the FASTA file together with Ians script named "translate_sequence_allframes_fasta.py"

python translate_sequence_allframes_fasta.py KS.fasta

# this creates a file called KS.fasta.faa with the aa sequences

# download HMM's from antismash website and put them into the same folder as the aa sequences   https://docs.antismash.secondarymetabolites.org/modules/nrps_pks_domains/

# run hmmpress for all HMMs (example below)

hmmpress PKS_KS.hmm.hmm

hmmpress ketoacyl-synt.hmm

# Use hmmsearch to search your genome for any protein sequences thatmatch the HMM profile

hmmsearch -o hmmer_ketoacyl-synt_KS --cpu 4 --tblout hmmer_ketoacyl-synt_KS.table ketoacyl-synt.hmm KS.fasta.faa

hmmsearch -o hmmer_PKS_KS_KS --cpu 4 --tblout hmmer_PKS_KS_KS.table PKS_KS.hmm KS.fasta.faa

hmmsearch -o hmmer_t2ks_all_KS --cpu 4 --tblout hmmer_t2ks_all_KS.table t2ks_all.hmm KS.fasta.faa

hmmsearch -o hmmer_A-OX_AD --cpu 4 --tblout hmmer_A-OX_AD.table A-OX.hmm AD.fasta.faa

hmmsearch -o hmmer-AMP_binding_AD --cpu 4 --tblout hmmer_AMP-binding_AD.table AMP-binding.hmm AD.fasta.faa

# Remove everything but the sequence IDs

grep -v ^# hmmer_ketoacyl-synt_KS.table | cut -f1 -d'.' > hmmer_ketoacyl-synt_KS.table.IDs

grep -v ^# hmmer_PKS_KS_KS.table | cut -f1 -d'.' > hmmer_PKS_KS_KS.table.IDs

grep -v ^# hmmer_t2ks_all_KS.table | cut -f1 -d'.' > hmmer_t2ks_all_KS.table.IDs

grep -v ^# hmmer_A-OX_AD.table | cut -f1 -d'.' > hmmer_A-OX_AD.table.IDs

grep -v ^# hmmer_AMP-binding_AD.table | cut -f1 -d'.' > hmmer_AMP-binding_AD.table.IDs

# write a new fasta file with the nucleotide sequences where the hmmsearch have filtered the sequences which didn't match
# OBS - the new file will (as an example) be named hmmer_ketoacyl-synt_KS.table.IDs.fasta

python sequence_subset.py hmmer_ketoacyl-synt_KS.table.IDs KS.fasta 

python sequence_subset.py hmmer_PKS_KS_KS.table.IDs KS.fasta 

python sequence_subset.py hmmer_t2ks_all_KS.table.IDs KS.fasta

python sequence_subset.py hmmer_A-OX_AD.table.IDs AD.fasta

python sequence_subset.py hmmer_AMP-binding_AD.table.IDs AD.fasta

# Merge the fasta files and remove redundant sequences !OBS! remember to change input and output files#

cat hmmer_t2ks_all_KS.table.IDs.fasta hmmer_PKS_KS_KS.table.IDs.fasta hmmer_ketoacyl-synt_KS.table.IDs.fasta | awk -v RS=">" '!a[$0]++ { print ">"$0; }' - | grep -Ev '^\s*$|^>\s*$' > merged_KS.fasta

cat hmmer_A-OX_AD.table.IDs.fasta hmmer_AMP-binding_AD.table.IDs.fasta | awk -v RS=">" '!a[$0]++ { print ">"$0; }' - | grep -Ev '^\s*$|^>\s*$' > merged_AD.fasta


# now we can use Vsearch to cluster the sequences at 85% or 0.85 into operational biosynthetic units (OBUs)

# example
vsearch --cluster_size merged_KS.fasta --id 0.85 --uc OBU_merged_KS.txt --otutabout OBU_merged_KS_KS.txt  

vsearch --cluster_size merged_AD.fasta --id 0.85 --uc OBU_merged_AD.txt --otutabout OBU_merged_AD_AD.txt  

# all columns except sequence ID were removed from OTUtable and the new file was saved as OBUKSID.txt#
# use only 1 representative OBU from each cluster by using the script "sequenceExtractor.py)
# OBS remember to change the file names in the python script to filter the correct file

python sequenceExtractor.py OBU_merged_KS_KS.txt KS.fasta OBU_KS_filtered.fasta

python sequenceExtractor.py OBU_merged_AD_AD.txt AD.fasta OBU_AD_filtered.fasta

# download MIBig database in Fasta format and make the database in diamond
wget https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_1.3.fasta
# rename mibig database to MIBif.faa (manualy)
# Make the database 
diamond makedb --in MIBig.faa -d MIBig

# Blast the filtered OBUs against the Mibig database
diamond blastx -d MIBig -q OBU_KS_filtered.fasta --sensitive -f 0 -o OBU_KS_KS.txt --id 85 --threads 6

diamond blastx -d MIBig -q OBU_AD_filtered.fasta --sensitive -f 6 -o OBU_AD_AD.txt --id 85 --threads 6

# Now I want to make the OTUtable into a OBU table where all the abundance data for each cluster is merged
#
# the next part is a bit tricky so I did some of it in excel
# Take the OTUtable that were used for R analysis containing all ASVs and transpose it it excel
# Give the column containing the sequences the name "ASV_ID"
# Take the file called KSdsclass and give the column containing the sequences the name "ASV_ID".

# Upload the transposed OTUtable and the KSdsclass files into R and merge the two files by the sequences
	
merged <- merge( ks_otutable_t, KSdsclass, by = c("ASV_ID"), all = FALSE)
# Write a .CSV file
	
 write.csv(merged, file="merged_ks_otutable.csv")
# open the file in excel
# In the "merged_ks_otutable.csv" delete all rows that should not be in a OTU table and replace the sequences with the names of the ASVs (Now you should have a transposed OTUtable with the names of the ASVs instead of the sequences. Save it as "ks_otutable_IDisASV" as a .csv 

# take the OBU_merged_KS.txt from the vsearch output and delete the bottom rows where the first letter is C (Manually).
# Give the columns names (cluster, cluster_number, ASV_ID and centroid_id and save it as clustertable as a .csv file
# Run the R script "makeOBUtable - this will create the OBU table without IDs.

# Sort the cluster table file and pull the centroid names and transfer them to OBU table. Check if there is any NAs and if they don't have abundance then delete them.

# Transpose the OBU table and save it as final_KS_OBU_table as a .csv

# now we want to assign taxonomy based on the non-redudant protein database and MEGAN
# we need to assign taxonomy to the OBUs by aligning them to the NCBI non redundant protein database.
# Move the fasta file (OBU_AD_filtered.fasta & OBU_KS_filtered.fasta) to the folder where the database is.
# Run the diamond blastx with detached screen
# examples
screen -L  diamond blastx -d nr -q "name.of.file.to.blast".fasta -o matches.m8 --threads 6

screen -L  diamond blastx -d nr -q OBU_KS_filtered.fasta -o OBU_KS_filtered.m8 --threads 6

screen -L  diamond blastx -d nr -q OBU_AD_filtered.fasta -o OBU_AD_filtered.m8 --threads 6

#Download Megan6 and install on server

wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_17_0.sh
chmod + x MEGAN_Community_unix_6_17_0.sh./MEGAN_Community_unix_6_17_0.sh

# Downloading taxonomy for Megan

wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/prot_acc2tax-Jul2019X1.abin.zip
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/nucl_acc2tax-Jul2019.abin.zip
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/acc2interpro-Jul2019X.abin.zip
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/acc2eggnog-Jul2019X.abin.zip 
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/acc2vfdb-Feb2019.map.gz
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/acc2seed-May2015XX.abin.zip

# load the resulting blast files into MEGAN and export the results from MEGAN
# export --> export txt(CSV) --> reads_to_taxon_path
# split the taxonomy in excel use https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi as reference
# data --> import from txt --> choose semicolon as deliminator.
# filter taxa manually...
# Merge the old ASV table with the taxonomy by ID to create the taxonomy file
# example
 merged <- merge( merge2, merge1, by = c('ID'), all = FALSE)
 write.csv(merged, file="OBU_KS_taxon.csv")
 
# now we want to make tables for phyloseq 
# Use this when you need to make phyloseq from two manual made tax table and OBU table

# this is site by species matrix, need row.names = 1 to have phyloseq read it

final_KS_OBU_table <- read.csv("final_KS_OBU_table.csv", header = T, row.names = 1, sep = ';')

# transpose to get a species x site matrix
final_KS_OBU_table_transposed <- t(downsampled_KS_OTUtable)

# NEVER use numbers as a name for Column or Row label in R, R puts an X in front of it
# get rid of Xs that were inserted in front of numbers row.names

rownames(final_KS_OBU_table_transposed) <-gsub("X","",row.names(final_KS_OBU_table_transposed))

# need this to be a matrix
 class(final_KS_OBU_table_transposed)

# make compatible for phyloseq format
final_KS_OBU_table = otu_table(final_KS_OBU_table, taxa_are_rows = TRUE)

# Read taxonomy info in
taxonomy <- read.csv("OBU_KS_taxon.csv.csv", row.names = 1, na.strings=c("","NA"), sep = ";")

# Needs to be a matrix
class(taxonomy)
taxonomy <- as.matrix(taxonomy)

# Make compatible for phyloseq
KS_taxonomy_final = tax_table(taxonomy)

##
OBS IF YOU GET ERROR TAXA/OTU DOES NOT MATCH, RUN CODE:
row.names(KS_taxonomy_final) <- KS_taxonomy_final[,sequence]
head(row.names(KS_taxonomy_final)) 
##

# add sample data (population)
samdfKS <-read.csv("KSsampledata.txt",sep='\t')
rownames(samdfKS) <- samdfKS[, 1]
samdfKS[, 1] <- NULL


#make phyloseq object (with sampledata)
KS_OBU_ps.true <- phyloseq(otu_table(final_KS_OBU_table, taxa_are_rows=TRUE), sample_data(samdfKS), tax_table(KS_taxonomy_final))

#Save phyloseq object
saveRDS(KS_OBU_ps.true, file="KS_OBU_ps.true_final.rds")

# now we want to make a ampvis2 object 
# Combine OBU abundance table and taxonomy table from the phyloseq object "ps":

obutable <- data.frame(OTU = rownames(phyloseq::otu_table(KS_OBU_ps.true)@.Data),
                       phyloseq::otu_table(KS_OBU_ps.true)@.Data,
                       phyloseq::tax_table(KS_OBU_ps.true)@.Data,
                       check.names = FALSE)
# Combine OBU table (Which includes taxonomy information) and load into Ampvis2 
KS_OBU_amp2 <- amp_load(otutable, KSsampledata)

# Data analysis in R
# loading the data

setwd("~/OBU/phyloseq")
KS_OBU_ps <- readRDS("~/OBU/phyloseq/KS/KS_OBU_ps.rds") 
KS_OBU_amp2 <- readRDS("~/OBU/phyloseq/KS/KS_OBU_amp2.rds")
AD_OBU_ps <- readRDS("~/OBU/phyloseq/AD/AD_OBU_ps.rds")
AD_OBU_amp2 <- readRDS("~/OBU/phyloseq/AD/AD_OBU_amp2.rds")
df_ks<-read.csv("~/OBU/phyloseq/KS/KS_OBUtable_phyloseq.csv", row.names=1, sep=";")
names_KS<-read.delim("~/OBU/phyloseq/KS/KSsampledata.txt", row.names=1, comment.char="#")
df_AD<-read.csv("~/OBU/phyloseq/AD/AD_OBUtable_phyloseq.csv", row.names=1, sep=";")
names_AD<-read.delim("~/OBU/phyloseq/AD/ADsampledata.txt", row.names=1, comment.char="#")

# rename Etosha to Otavi in all files

levels(AD_OBU_amp2$metadata$Nest) <- c("Betta", "Otavi", "Karasburg")
levels(AD_OBU_ps@sam_data$Nest) <- c("Betta", "Otavi", "Karasburg")
levels(KS_OBU_amp2$metadata$Nest) <- c("Betta", "Otavi", "Karasburg")
levels(KS_OBU_ps@sam_data$Nest) <- c("Betta", "Otavi", "Karasburg")

# Color palette of choice
color<-wes_palette(name="Darjeeling1", n=18, type="continuous")
color

# Richness diversity measures Ketosynthase
plot_richness(KS_OBU_ps, measures = c ("Observed", "Shannon") ,color = "Nest")+
  scale_x_discrete(limits=c("HK-KS-A","HK-KS-D","HK-KS-H","HK-KS-B","HK-KS-C", "HK-KS-G","HK-KS-I","HK-KS-L","HK-KS-M","HK-KS-E","HK-KS-F","HK-KS-J","HK-KS-K"))+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        theme(axis.title.x = element_text(face = "bold"))+
        theme(axis.title.y = element_text(face = "bold"))+
        geom_point(size=3)

# Making stacked bar plots
# NOTE: Remember to change level and name corresponding to the plot you want to make, as well as colors.
# stacked bar plots kingdom KS
stacked_barplot(KS_OBU_ps, other_threshold = 0.03, level=1) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1], "Fungi"=color[3],"Bacteria"=color[16]), name="Kingdom")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
  
# stacked bar plots phylum KS
   stacked_barplot(KS_OBU_ps, other_threshold = 0.03, level=2) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1], "Unclassified"=color[17],"Ascomycota"=color[2],"Actinobacteria"=color[16],
                                   "Proteobacteria"=color[3],"Bacteroidetes"=color[15]), name="Phylum")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
        
# stacked bar plots class KS
stacked_barplot(KS_OBU_ps, other_threshold = 0.03, level=3) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17], "Actinobacteria"=color[2],"Sordariomycetes"=color[16],"Dothideomycetes"=color[3],"Leotiomycetes"=color[15],"Cytophagia"=color[4]), name="Class")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# stacked bar plots order KS       
stacked_barplot(KS_OBU_ps, other_threshold = 0.03, level=4) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17], "Hypocreales"=color[2],"Pleosporales"=color[16],"Erysiphales"=color[3],"Cytophagales"=color[15],"Dothideales"=color[4]), name="Order")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# stacked bar plots family KS
stacked_barplot(KS_OBU_ps, other_threshold = 0.03, level=5) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17], "Hypocreaceae"=color[2],"Pleosporaceae"=color[16],"Erysiphaceae"=color[3],"Saccotheciaceae"=color[15],"Didymosphaeriaceae"=color[4],"Ophiocordycipitaceae"=color[14]), name="Family")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# stacked bar plots Genus KS
stacked_barplot(KS_OBU_ps, other_threshold = 0.03, level=6) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17], "Trichoderma"=color[2],"Golovinomyces"=color[16],"Paraphaeosphaeria"=color[3],"Alternaria"=color[15],"Purpureocillium"=color[4],"Aureobasidium"=color[14],"Bipolaris"=color[5],"Stemphylium"=color[18]), name="Genus")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# Stacked barplot Species KS
stacked_barplot(KS_OBU_ps, other_threshold = 0.03, level=7) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17], "Golovinomyces magnicellulatus"=color[2],"Paraphaeosphaeria sporulosa"=color[16],"Alternaria oxytropis"=color[3],"Purpureocillium lilacinum"=color[15],"Stemphylium lycopersici"=color[4]), name="Species")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# Fractional abundance of KS OBUs
df_ks<-t(df_ks)
colnames(df_ks)=names_KS$Nest
#df_ks<-t(rowsum(t(df_ks), group = colnames(df_ks), na.rm = T))
sum_ks=data.frame(colSums(df_ks))
df_ks<-data.frame((df_ks))
df_ks<-df_ks %>% mutate(Bettafraq_ks=Betta/sum_ks[1,1]*100) %>% mutate(Otavifraq_ks=Otavi/sum_ks[3,1]*100) %>% mutate(Karasburgfraq_ks=Karasburg/sum_ks[2,1]*100)
colSums(df_ks)
# Rank abundance plot KS OBUs Betta
bettaframe_ks<-df_ks[order(-df_ks$Bettafraq_ks),]
bettaframe_ks<-bettaframe_ks %>% mutate(names_KS=rownames(bettaframe_ks))
bettaframe_ks<-bettaframe_ks[1:20,4:7]
orders<-bettaframe_ks$names_KS

bettaframe_ks<-bettaframe_ks %>% mutate(highlight_flag = names_KS%in%c("KS.ASV.1","KS.ASV.1013","KS.ASV.1017","KS.ASV.1084","KS.ASV.1115","KS.ASV.1134","KS.ASV.1153","KS.ASV.1210","KS.ASV.1214","KS.ASV.1217","KS.ASV.183","KS.ASV.394","KS.ASV.42")) 
plot1<-ggplot()+
  geom_col(data=bettaframe_ks, aes(x=names_KS[1:20], y=Bettafraq_ks[1:20], fill=highlight_flag))+
  scale_x_discrete(limits = orders)+
  scale_y_continuous(name="Relative OBU abundance (%)", limits = c(0,60))+
  scale_fill_manual(values = c("Black", "blue"))+
  ggtitle("Betta")+
  theme(axis.text.x = element_text(angle=90,hjust = 1 ,vjust = 0.5, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color=NA),
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.position = "none",
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x =element_text(color=NA )
        )
plot1

# Rank abundance plot KS OBUs Otavi
Otaviframe_ks<-df_ks[order(-df_ks$Otavifraq_ks),]
Otaviframe_ks<-Otaviframe_ks %>% mutate(names_KS=rownames(Otaviframe_ks))
Otaviframe_ks<-Otaviframe_ks[1:20,4:7]
orders<-Otaviframe_ks$names_KS


Otaviframe_ks<-Otaviframe_ks %>% mutate(highlight_flag = names_KS%in%c("KS.ASV.1","KS.ASV.1013","KS.ASV.1017","KS.ASV.1084","KS.ASV.1115","KS.ASV.1134","KS.ASV.1153","KS.ASV.1210","KS.ASV.1214","KS.ASV.1217","KS.ASV.183","KS.ASV.394","KS.ASV.42"))  
plot2<-ggplot()+
  geom_col(data=Otaviframe_ks, aes(x=names_KS[1:20], y=Otavifraq_ks[1:20], fill=highlight_flag))+
  scale_x_discrete(limits = orders)+
  scale_y_continuous(name="Relative OBU abundance (%)",limits = c(0,60))+
  scale_fill_manual(values = c("Black", "blue"))+
  ggtitle("Otavi")+
  theme(axis.text.x = element_text(angle=90,hjust = 1 ,vjust = 0.5, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color=NA), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.position = "none",
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x =element_text(color=NA )
        )

plot2
# Rank abundance plot KS OBUs Karasburg
karasframe_ks<-df_ks[order(-df_ks$Karasburgfraq_ks),]
karasframe_ks<-karasframe_ks %>% mutate(names_KS=rownames(karasframe_ks))
karasframe_ks<-karasframe_ks[1:20,4:7]
orders<-karasframe_ks$names_KS


karasframe_ks<-karasframe_ks %>% mutate(highlight_flag = names_KS%in%c("KS.ASV.1","KS.ASV.1013","KS.ASV.1017","KS.ASV.1084","KS.ASV.1115","KS.ASV.1134","KS.ASV.1153","KS.ASV.1210","KS.ASV.1214","KS.ASV.1217","KS.ASV.183","KS.ASV.394","KS.ASV.42")) 
plot3<-ggplot()+
  geom_col(data=karasframe_ks, aes(x=names_KS[1:20], y=Karasburgfraq_ks[1:20], fill=highlight_flag))+
  scale_x_discrete(limits = orders)+
  scale_y_continuous(name="Relative OBU abundance (%)",limits = c(0,60))+
  scale_fill_manual(values = c("Black", "blue"))+
  ggtitle("Karasburg")+
  theme(axis.text.x = element_text(angle=90,hjust = 1 ,vjust = 0.5, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color=NA), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.position = "none",
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x =element_text(color=NA )
        )

plot3
# venn diagram KS
amp_venn(KS_OBU_amp2, group_by = "Nest", cut_a = 0, cut_f = 1)+
  ggtitle("KS OBUs")
# Principal coordinate analysis Ketosynthase
amp_ordinate(KS_OBU_amp2, type="PCOA", distmeasure = "bray", transform = "none", 
             filter_species = 0, species_label_taxonomy = "none", print_caption = TRUE, 
             sample_color_by = "Nest", sample_colorframe = "Nest")+
        theme(axis.text.x = element_text(angle=0,hjust = 1 ,vjust = 1, face = "bold", color = "Black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12))+
        theme(axis.title.x = element_text(face = "bold"))+
        theme(axis.title.y = element_text(face = "bold"))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Anosim analysis of ketosynthase domains
anosim_KS<-as.matrix(KS_OBU_ps@otu_table@.Data)
anosim(t(anosim_KS), grouping=KS_OBU_ps@sam_data$Nest, distance = "bray", permutations = 9999)

# Richness plot Adenylation domains
plot_richness(AD_OBU_ps, measures = c ("Observed", "Shannon") ,color = "Nest")+
  scale_x_discrete(limits=c("HK-AD-A","HK-AD-D","HK-AD-H","HK-AD-B","HK-AD-C", "HK-AD-G","HK-AD-I","HK-AD-L","HK-AD-M","HK-AD-E","HK-AD-F","HK-AD-J","HK-AD-K"))+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        theme(axis.title.x = element_text(face = "bold"))+
        theme(axis.title.y = element_text(face = "bold"))+
        geom_point(size=3)
 # Stacked bar plot Kingdom Adenylation domain
 stacked_barplot(AD_OBU_ps, other_threshold = 0.03, level=1) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1], "Unclassified"=color[17], "Fungi"=color[3],"Bacteria"=color[16]), name="Kingdom")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# Stacked bar plot Phylum Adenylation domain
stacked_barplot(AD_OBU_ps, other_threshold = 0.03, level=2) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1], "Unclassified"=color[17],"Ascomycota"=color[2],"Actinobacteria"=color[16],
                                   "Proteobacteria"=color[3],"Bacteroidetes"=color[15],"Firmicutes"=color[4]), name="Phylum")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# Stacked bar plot class Adenylation domain
stacked_barplot(AD_OBU_ps, other_threshold = 0.03, level=3) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17], "Actinobacteria"=color[2],"Sordariomycetes"=color[16],"Dothideomycetes"=color[3],"Eurotiomycetes"=color[15],"Cytophagia"=color[4],"Deltaproteobacteria"=color[18],"Gammaproteobacteria"=color[5],"Alphaproteobacteria"=color[13],"Betaproteobacteria"=color[6]), name="Class")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
 # stacked bar plot order Adenylation domain
 stacked_barplot(AD_OBU_ps, other_threshold = 0.03, level=4) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17],        "Corynebacteriales"=color[2],"Pleosporales"=color[16],"Eurotiales"=color[3],
        "Dothideales"=color[15],"Pseudonocardiales"=color[4],"Myxococcales"=color[18],
        "Cytophagales"=color[5],"Enterobacterales"=color[13],"Streptomycetales"=color[6],
        "Capnodiales"=color[12],"Hypocreales"=color[7],"Burkholderiales"=color[11],"Micrococcales"=color[8],
        "Pseudomonadales"=color[10],"Rhodospirillales"=color[9]), name="Order")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# Stacked bar plot family AD
stacked_barplot(AD_OBU_ps, other_threshold = 0.03, level=5) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17],         "Nocardiaceae"=color[2],"Pleosporaceae"=color[16],"Aspergillaceae"=color[3],
        "Saccotheciaceae"=color[15],"Pseudonocardiaceae"=color[4],"Yersiniaceae"=color[18],
        "Streptomycetaceae"=color[5],"Hymenobacteraceae"=color[13],"Archangiaceae"=color[6],
        "Oxalobacteraceae"=color[12],"Pseudomonadaceae"=color[7],"Rhodospirillaceae"=color[11],"Micrococcaceae"=color[8],
        "Burkholderiaceae"=color[10],"Cordycipitaceae"=color[9],"Ophiocordycipitaceae"=color[14]),name="Family")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
 # Stacked bar plot family AD
 stacked_barplot(AD_OBU_ps, other_threshold = 0.03, level=6) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17], "Rhodococcus"=color[2],"Aureobasidium"=color[16],"Serratia"=color[3],"Hymenobacter"=color[15],"Streptomyces"=color[4],"Pseudomonas"=color[14],"Skermanella"=color[5],"Lecanicillium"=color[18],"Alternaria"=color[6],"Purpureocillium"=color[14],"Kocuria"=color[7]), name="Genus")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# Stacked bar plot species AD
stacked_barplot(AD_OBU_ps, other_threshold = 0.03, level=7) +
  facet_wrap(KS_OBU_ps@sam_data@.Data[[1]], scales = "free_x")+
  theme(axis.text.x = element_text(angle=45,hjust = 1 ,vjust = 1, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color="black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12)
        )+
        scale_fill_manual(values=c("Others"=color[1],"Unclassified"= color[17], "Rhodococcus kyotonensis"=color[2],"Skermanella aerolata"=color[16],"Purpureocillium lilacinum"=color[3]), name="Species")+
        scale_y_continuous(name="Relative OBU abundance")+
        scale_x_discrete(name = "")+
        theme(axis.title.y = element_text(face = "bold"))
# Rank abundance plot AD fractional
df_AD<-t(df_AD)
colnames(df_AD)=names_AD$Nest
df_AD<-t(rowsum(t(df_AD), group = colnames(df_AD), na.rm = T))

sum_AD=data.frame(colSums(df_AD))
df_AD<-data.frame((df_AD))
df_AD<-df_AD %>% mutate(Bettafraq_AD=Betta/sum_AD[1,1]*100) %>% mutate(Otavifraq_AD=Otavi/sum_AD[3,1]*100) %>% mutate(Karasburgfraq_AD=Karasburg/sum_AD[2,1]*100)
colSums(df_AD)

# rank abundance plot AD OBUs Betta
bettaframe_AD<-df_AD[order(-df_AD$Bettafraq_AD),]
bettaframe_AD<-bettaframe_AD %>% mutate(names_AD=rownames(bettaframe_AD))
bettaframe_AD<-bettaframe_AD[1:20,4:7]
orders<-bettaframe_AD$names_AD


bettaframe_AD<-bettaframe_AD %>% mutate(highlight_flag = names_AD%in%c("AD.ASV.1","AD.ASV.1008","AD.ASV.1011","AD.ASV.1013","AD.ASV.1065","AD.ASV.1088","AD.ASV.1114","AD.ASV.1115","AD.ASV.112","AD.ASV.1123","AD.ASV.1230","AD.ASV.1412","AD.ASV.1442","AD.ASV.146","AD.ASV.1523","AD.ASV.2157","AD.ASV.3049","AD.ASV.3081","AD.ASV.360","AD.ASV.3770","AD.ASV.4706","AD.ASV.916")) 
plot11<-ggplot()+
  geom_col(data=bettaframe_AD, aes(x=names_AD[1:20], y=Bettafraq_AD[1:20], fill=highlight_flag))+
  scale_x_discrete(limits = orders)+
  scale_y_continuous(name="Relative OBU abundance (%)", limits = c(0,60))+
  scale_fill_manual(values = c("Black", "blue"))+
  ggtitle("Betta")+
  theme(axis.text.x = element_text(angle=90,hjust = 1 ,vjust = 0.5, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color=NA),
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.position = "none",
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x =element_text(color=NA )
        )
plot11
# rank abundance plot AD OBUs Otavi
Otaviframe_AD<-df_AD[order(-df_AD$Otavifraq_AD),]
Otaviframe_AD<-Otaviframe_AD %>% mutate(names_AD=rownames(Otaviframe_AD))
Otaviframe_AD<-Otaviframe_AD[1:20,4:7]
orders<-Otaviframe_AD$names_AD
Otaviframe_AD<-Otaviframe_AD %>% mutate(highlight_flag = names_AD%in%c("AD.ASV.1","AD.ASV.1008","AD.ASV.1011","AD.ASV.1013","AD.ASV.1065","AD.ASV.1088","AD.ASV.1114","AD.ASV.1115","AD.ASV.112","AD.ASV.1123","AD.ASV.1230","AD.ASV.1412","AD.ASV.1442","AD.ASV.146","AD.ASV.1523","AD.ASV.2157","AD.ASV.3049","AD.ASV.3081","AD.ASV.360","AD.ASV.3770","AD.ASV.4706","AD.ASV.916"))   
plot22<-ggplot()+
  geom_col(data=Otaviframe_AD, aes(x=names_AD[1:20], y=Otavifraq_AD[1:20], fill=highlight_flag))+
  scale_x_discrete(limits = orders)+
  scale_y_continuous(name="Relative OBU abundance (%)",limits = c(0,60))+
  scale_fill_manual(values = c("Black", "blue"))+
  ggtitle("Otavi")+
  theme(axis.text.x = element_text(angle=90,hjust = 1 ,vjust = 0.5, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color=NA), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.position = "none",
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x =element_text(color=NA )
        )
plot22
# rank abundance plot AD OBUs Karasburg
karasframe_AD<-df_AD[order(-df_AD$Karasburgfraq_AD),]
karasframe_AD<-karasframe_AD %>% mutate(names_AD=rownames(karasframe_AD))
karasframe_AD<-karasframe_AD[1:20,4:7]
orders<-karasframe_AD$names_AD
karasframe_AD<-karasframe_AD %>% mutate(highlight_flag = names_AD%in%c("AD.ASV.1","AD.ASV.1008","AD.ASV.1011","AD.ASV.1013","AD.ASV.1065","AD.ASV.1088","AD.ASV.1114","AD.ASV.1115","AD.ASV.112","AD.ASV.1123","AD.ASV.1230","AD.ASV.1412","AD.ASV.1442","AD.ASV.146","AD.ASV.1523","AD.ASV.2157","AD.ASV.3049","AD.ASV.3081","AD.ASV.360","AD.ASV.3770","AD.ASV.4706","AD.ASV.916"))  
plot33<-ggplot()+
  geom_col(data=karasframe_AD, aes(x=names_AD[1:20], y=Karasburgfraq_AD[1:20], fill=highlight_flag))+
  scale_x_discrete(limits = orders)+
  scale_y_continuous(name="Relative OBU abundance (%)",limits = c(0,60))+
  scale_fill_manual(values = c("Black", "blue"))+
  ggtitle("Karasburg")+
  theme(axis.text.x = element_text(angle=90,hjust = 1 ,vjust = 0.5, face = "bold"), 
        panel.background = element_rect(fill="white", color="black"), 
        panel.grid.major = element_line (color=NA), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.position = "none",
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title.y = element_text(face="bold"),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12),
        axis.title.x =element_text(color=NA )
        )

plot33
# putting all 6 rankabundance plots into 1:
gridExtra::grid.arrange(plot1,plot2, plot3, plot11, plot22, plot33, ncol=3, nrow=2)

# Venn diagram AD
amp_venn(AD_OBU_amp2, group_by = "Nest", cut_a = 0, cut_f = 1)+
  ggtitle("AD OBUs")
# Principal coordinate analysis Adenylation domains
amp_ordinate(AD_OBU_amp2, type="PCOA", distmeasure = "bray", transform = "none", 
             filter_species = 0, species_label_taxonomy = "none", print_caption = TRUE, 
             sample_color_by = "Nest", sample_colorframe = "Nest")+
        theme(axis.text.x = element_text(angle=0,hjust = 1 ,vjust = 1, face = "bold", color = "Black"), 
        axis.text = element_text(color="black", face = "bold"),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(face = "bold", size = 12))+
        theme(axis.title.x = element_text(face = "bold"))+
        theme(axis.title.y = element_text(face = "bold"))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Anosim statisctics AD
anosim_AD<-as.matrix(AD_OBU_ps@otu_table@.Data)
anosim(t(anosim_AD), grouping=AD_OBU_ps@sam_data$Nest, distance = "bray", permutations = 9999)
