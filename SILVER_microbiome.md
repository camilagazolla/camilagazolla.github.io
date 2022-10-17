---
title: SILVER microbiome analysis
---

# Demultiplexing steps

We start our journey with multiplexed FASTQ files. 

#### Create mapping file 
We need to create a .txt with connection of samples to barcodes such as [map_2G0089_16S_LibA.txt](https://github.com/camilagazolla/camilagazolla.github.io/blob/main/map_2G0089_16S_LibA.txt). You need to use a file like [2G0089_PCR Profile.xlsx](https://github.com/camilagazolla/camilagazolla.github.io/blob/main/2G0089_PCR%20Profile.xlsx) that is supposed to be at the Dropbox in a folder such as Personal/Next_Generation_Sequencing_Assays/Run/Library. The mapping file needs Linux line breaks. On Notepad++ you can go to edit --> EOL Conversion --> choose LINUX. You can make sure that $ appear at the end with:

```
cat -e map_2G0089_16S_LibA.txt
```

#### Ask for resources on the cluster
```
srun --partition=large-mem --nodes=1 --ntasks=1 --cpus-per-task=1 --mem-per-cpu 800G -t 00-24:00:00 --pty bash
```

#### Place the ITS and 16S library files in separated folders, DO IT FOR EACH LIBRARY TOO!
```
cd [.fastq.gz] # library folder with ITS **OR** 16S 
gunzip *.fastq.gz
```

#### Use the 1_map.sh script to create files
Use [1_map.sh](https://github.com/camilagazolla/camilagazolla.github.io/blob/main/1_map.sh), but ATTENTION using it with different libraries because the same filenames will be created.

```
sh 1_map.sh -m map_2G0089_16S_LibA.txt
```
The outputs should be: 1_1_1_plate_barcodes.txt, 1_1_3_forward.txt, general_barcode_p1.txt, and general_barcode_p2.txt.

#### Use the 2_demultiplex_HPC_ol_skipper_SL.sh script
To use [2_demultiplex_HPC_ol_skipper_SL.sh](https://github.com/camilagazolla/camilagazolla.github.io/blob/main/2_demultiplex_HPC_ol_skipper_SL.sh) you should have in the same folder the 4 outputs from 1_map.sh, 2_demultiplex_HPC_ol_skipper_SL.sh and map_2G0089_16S_LibA.txt.
I recommend to run it with an open screen session:

```
screen -S 2G0089_16S_LibA

srun --partition=large-mem --nodes=1 --ntasks=1 --cpus-per-task=1 --mem-per-cpu 200G -t 00-24:00:00 --pty bash

cd /gs/gsfs0/users/cgazollavo/SILVER/raw_data/2G0089/16S/LibA/ # folder with the files

# you need to give the full path to the fastq files
sh 2_demultiplex_HPC_ol_skipper_SL.sh \
-f /gs/gsfs0/users/cgazollavo/SILVER/raw_data/2G0089/16S/LibA/p8916SV4A_R1_001.fastq \
-r /gs/gsfs0/users/cgazollavo/SILVER/raw_data/2G0089/16S/LibA/p8916SV4A_R2_001.fastq \
-m map_2G0089_16S_LibA.txt 
```

#### Move and rename the samples with the library and run name
```
# cd to the samples folder
cd /gs/gsfs0/users/cgazollavo/SILVER/raw_data/2G0089/16S/LibA/samples

mkdir samples_fastq
mv */*.fastq samples_fastq/
ls -l . | grep -c ^d # count the number of directory (samples+1)
cd samples_fastq 
ls | wc -l # count the number of files

# CHANGE HERE!!!
for FILENAME in *; do mv $FILENAME 2G0089_16S_LibA_$FILENAME; done # rename files

cd ..
ls . | grep -v "samples_fastq" | xargs rm -r
mv samples_fastq/*.fastq .
rm samples_fastq -r
```

#### Evaluate if the demultiplexing worked
IMPORTANT: Demultiplexing could go wrong! You need to make sure that no sample was left behind. Go to the folder with the samples and ls the files:
```
cd /gs/gsfs0/users/cgazollavo/SILVER/raw_data/2G0089/16S/LibA/samples

# CHANGE HERE!!!
ls *R1.fastq > 2G0089_16S_LibA_samples.txt
```

You can put all the files produced in multiple runs in a folder to create something like [2G0089_16S_files_dem.txt](https://github.com/camilagazolla/camilagazolla.github.io/blob/main/2G0089_16S_files_dem.txt)
```
cat 2G0089_16S_Lib* > 2G0089_16S_files_dem.txt
```

Now use R to check if all samples were processed:

```
setwd("/Users/camila/Desktop/SILVER/microbiome_data/demultiplexing")
library(dplyr)

# read the data
files16s <- read.table("2G0089_16S_files_dem.txt", header = F)
head(files16s)

# separate the colum
files16s$V1 <- gsub("_R..fastq","",files16s$V1)
files16s <- separate(files16s, V1, c("Amp","Run","RunLib","BurkID"), sep="\\.", remove = FALSE)

# remove duplicates
files16s <- files16s[!duplicated(files16s$V1),]

View(files16s)
table(files16s$Run)
table(files16s$RunLib)

# read the map files
libA_89 <- read.delim("map_2G0089_16S_LibA.txt", header = T)
libB_89 <- read.delim("map_2G0089_16S_LibB.txt", header = T)

libA_93 <- read.delim("map_2G0093_16S_LibA.txt", header = T)
libB_93 <- read.delim("map_2G0093_16S_LibB.txt", header = T)

# see the differences
length(libA_89$Lab.ID)
dim(files16s %>% filter(Run=="2G0089" & RunLib=="libA"))[1] # SIZES OK

length(libB_89$Lab.ID)
dim(files16s %>% filter(Run=="2G0089" & RunLib=="libB"))[1] # SIZES OK

length(libA_93$Lab.ID) # DEU RUIM
dim(files16s %>% filter(Run=="2G0093" & RunLib=="libA"))[1] # SIZES OK

length(libB_93$Lab.ID)
dim(files16s %>% filter(Run=="2G0093" & RunLib=="libB"))[1] # SIZES OK
```

# Evaluating quality

We start with the paired-end fastq files that have been demultiplexed. To inspect read quality profiles for each SILVER run we are going to use [FastQC](https://github.com/s-andrews/FastQC) to generate log files which are going to be passed to [MultiQC](https://github.com/ewels/MultiQC). A final report will be generated will containg summarising  statistics for all samples.

On Shell:
```
 # activate the conda enviroment with the tools
conda activate biobase

# cd to the folder with fastq files
cd /gs/gsfs0/users/cgazollavo/SILVER/16S/2G00093/ 

mkdir forward_qc reverse_qc
fastqc *fastq -t 30
mv *_R1_fastqc* forward_qc
mv *_R2_fastqc* reverse_qc
multiqc forward_qc/ -o forward_qc/
multiqc reverse_qc/ -o reverse_qc/
```

Check the multiqc_report.html on each folder.

# ASV picking with DADA2 
Considering that every amplicon dataset has a different set of error rates, each one of the two SILVER sequecing runs (2G0089 and 2G0093) will be analysized separately until the classification step.

On R:
```{r}
# load packages
library("dada2")
library("phyloseq")
library("stringr")
library("DECIPHER")

# point to the folder of the run containing the fastq files
path_reads <- "/gs/gsfs0/users/cgazollavo/SILVER/16S/2G00089/"

# define filtering parameters according to fastqc+multiqc results
trimLeftF = 10 # value remove from the start of F (R1) reads
truncLenF = 130 # position to truncate F (R1) reads
trimLeftR = 10 # value remove from the start of R (R2) read
truncLenR = 130 # position to truncate R (R2) reads

# save trim options
A <- c("trimLeftF", "truncLenF", "trimLeftR", "truncLenR")
B <- c(trimLeftF, truncLenF, trimLeftR, truncLenR)
write(rbind(A,B), "trimoptions.txt")

# filtered forward files go into the filtered/ subdirectory
filtpath <- file.path(path_reads, "filtered_strict") 

# get filenames
fastqFs <- sort(list.files(path_reads, pattern="R1.fastq"))
fastqRs <- sort(list.files(path_reads, pattern="R2.fastq"))

# make sure we are only dealing with mate pairs
both_reads <- intersect(gsub("R1.fastq", "", fastqFs), gsub("R2.fastq", "", fastqRs))
names(fastqFs) <- gsub("R1.fastq", "", fastqFs)
names(fastqRs) <- gsub("R2.fastq", "", fastqRs)
fastqFs <- fastqFs[both_reads]
fastqRs <- fastqRs[both_reads]

# filtering and trimming
filterAndTrim_res <- filterAndTrim(fwd=file.path(path_reads, fastqFs), filt=file.path(filtpath, fastqFs),
              rev=file.path(path_reads, fastqRs), filt.rev=file.path(filtpath, fastqRs),
              truncLen=c(truncLenF,truncLenR), trimLeft=c(trimLeftF, trimLeftR),
              maxEE=c(2,2), matchIDs = TRUE, maxN=0, rm.phix=TRUE, compress=TRUE,
              verbose=TRUE, multithread=TRUE)
              
# file parsing
filtFs <- list.files(filtpath, pattern="R1.fastq", full.names = TRUE)
filtRs <- list.files(filtpath, pattern="R2.fastq", full.names = TRUE)

# sample names
sample_names <-gsub("(.*)_R1\\.fastq","\\1",basename(filtFs))
sample_namesR <- gsub("(.*)_R2\\.fastq","\\1",basename(filtRs))

# save environment  
save.image(file='env_1.RData')
```
🛑 STOP: make sure that the samples names are right by checking a few of them:
```{r}
sample_names[1:10]
sample_namesR [1:10]
```
If yes, continue the pipeline:
```{r}
# test if we have paired end data
if(!identical(sample_names, sample_namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample_names
names(filtRs) <- sample_names

set.seed(12345)
# learn forward error rates
errF <- learnErrors(filtFs, multithread=TRUE)

# learn reverse error rates
errR <- learnErrors(filtRs, multithread=TRUE)
 
# dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
# save read number for latter
getN <- function(x) sum(getUniques(x))
track_i <- cbind(filterAndTrim_res, sapply(dadaFs, getN), sapply(dadaRs, getN))
saveRDS(track_i, "track_i.rds")

# save environment  
save.image(file='env_2.RData')
}
```
🛑 STOP: test if the reads should be concatenated or merged. First, try to merge the reads:
```{r}
# merge denoised reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# check if reads were merged
cbind(filterAndTrim_res, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
```
If most of the reads were not merged (**this is the case of SILVER!**), concatenate them:
```{r}
# concatenate denoised reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE)
```
And continue the pipeline:
```{r}
# make the sequence table as save its raw version
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "seqtab_raw.RDS")

# drop reads present in less than 2 samples
seqtab <- seqtab[, colSums(seqtab) > 2]
write(dim(seqtab), "seqtab_drop_only_one_sample_dim.txt")

#check  
head(table(colSums(seqtab)))

saveRDS(seqtab, "seqtab_filt.RDS")
```

#### Taxonomic classification
In SILVER, more than one sequence table was generated. Here we are going to merge them, as well as the track files.

On a fresh R session:
```{r}
# merge data
track1 <- readRDS("track_i_2G00089.rds")
track2 <- readRDS("track_i_2G00093.rds")
track <- rbind(track1, track2)
saveRDS(track, "track_all.rds")

seqtab1 <- readRDS("seqtab_filt_2G00089.RDS")
seqtab2 <- readRDS("seqtab_filt_2G00093.RDS")
seqtab <- mergeSequenceTables(seqtab1, seqtab2)
saveRDS(seqtab, "seqtab_filt_all.rds")

# assign taxonomy
tax <- assignTaxonomy(seqtab, "/gs/gsfs0/users/burk-lab/DB/DADA2/silva_nr_v138_train_set.fa.gz", multithread=TRUE, verbose = TRUE)
saveRDS(tax, "tax.rds")

# create physeq object and save it
physeq_raw <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax))
saveRDS(tax, "physeq_raw.rds")

save.image(file='env_3.RData')
```
 
#### Filter the phyloseq object 
 
```
# load phyloseq
library(phyloseq)
physeq_raw <- readRDS("physeq_raw.rds")

# mantain only Bacteria Kingdom
table(tax_table(physeq_raw)[,"Kingdom"],useNA="always")
physeq_filt <-  phyloseq::subset_taxa(physeq_raw, Kingdom =="Bacteria")

# remove phylum of NA
table(tax_table(physeq_filt)[,"Phylum"],useNA="always")
physeq_filt <- subset_taxa(physeq_filt, Phylum != "NA")

# agglomerate to genus
physeq_filt <- tax_glom(physeq_filt, "Genus")

# replace NA with the last taxa assigned
library(zoo)
tax_tab <- tax_table(physeq_filt)
tax_tab <- na.locf(t(tax_tab), na.rm = FALSE)
tax_table(physeq_filt) <- t(tax_tab)

# filter phyoseq SILVER samples (not controls!) with less than 10K reads
lib_cut = 10000
sample_filt = sample_sums(physeq_raw)>=lib_cut

# if contains "NC, PC or SNEG" in the sample name add FALSE to sample_filt vector
for (i in 1:length(sample_filt)){
    if (grepl("PC|NC|SNEG", names(sample_filt[i]))){
        sample_filt[i] = TRUE
        } 
}

# filter phyloseq object
physeq_filt <- prune_samples(sample_filt, physeq_filt)

# remove taxa with all zeroes
physeq_filt <- prune_taxa(taxa_sums(physeq_filt) != 0, physeq_filt) 

# save phyloseq object
saveRDS(physeq_filt, file = "physeq_filt.rds")
```
