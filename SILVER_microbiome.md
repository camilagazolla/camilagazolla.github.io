---
title: SILVER microbiome analysis
---

# Evaluating quality

On Shell:
```
cd /gs/gsfs0/users/cgazollavo/SILVER/16S/2G00093/ # folder with fastq files
fastqc *fastq -t 30
mkdir forward_qc
mv *_R1_fastqc* forward_qc
mkdir reverse_qc
mv *_R2_fastqc* reverse_qc
multiqc forward_qc/
mv multiqc_report.html forward_qc/
mv multiqc_data/ forward_qc/
multiqc reverse_qc/
mv multiqc_report.html reverse_qc/
mv multiqc_data/ reverse_qc/
```

ADD IMAGE HERE


# ASV picking and taxonomic classification

We start with the paired-end fastq files that have been demultiplexed and with barcodes/adapters removed. 

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
\
\
🛑 STOP: make sure that the samples names are right by checking a few of them:
```{r}
sample_names[1:10]
sample_namesR [1:10]
```
\
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
\
\
🛑 STOP: test if the reads should be concatenated or merged. First, try to merge the reads:

```{r}
# merge denoised reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# check if reads were merged
cbind(filterAndTrim_res, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
```
\
If most of the reads were not merged (**this is the case of SILVER!**), concatenate them:
```{r}
# concatenate denoised reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate=TRUE)
```
\
And continue the pipeline:
```{r}
# make the sequence table as save its raw version
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "seqtab_raw.RDS")

# drop reads present in just one sample
seqtab <- seqtab[, colSums(seqtab) > 1]

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T, minFoldParentOverAbundance = 4)

# save environment  
save.image(file='env_3.RData')
```








## Requirements 

**NOTE:** This analysis requires at least 10GB of RAM to run.
It uses large files not included in the repository and many steps can take a few minutes to run. 

## Parameters

### Analysis input/output

```{r }
```
