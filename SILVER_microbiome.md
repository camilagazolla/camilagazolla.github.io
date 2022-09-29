---
title: SILVER microbiome analysis
---

# Input files

We start our journey with the paired-end fastq files that have been demultiplexed and with barcodes/adapters removed. 


# Evaluating quality


To inspect read quality profiles for each SILVER run we are going to use [FastQC](https://github.com/s-andrews/FastQC) to generate log files which are going to be passed to [MultiQC](https://github.com/ewels/MultiQC). A final report will be generated will containg summarising  statistics for all samples.


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

On R you can also use the plotQC function bellow to cancatenate the data of multiple sequence runs:
```
# load function to plot QC data from multiQC
plotQC <- function(multiqc_dfs=list()){
  require("ggplot2")
  require("reshape2")
  
  for (df_name in names(multiqc_dfs)){
    df <- multiqc_dfs[[df_name]]
    df <- df[c(-1)]
    meandf <- mean(as.matrix(df))
    sddf <- sd(as.matrix(df))
    cat(df_name, "present mean of:", meandf, "and sd of:", sddf, "\n") 
  }
  
  df_plotFinal<- data.frame()
  for (df_name in names(multiqc_dfs)){
    df <- multiqc_dfs[[df_name]]
    df <- df[c(-1)]
    dfmean <- mean(as.matrix(df))
    dfsd <- sd(as.matrix(df))
    df_plot<- colMeans(df)
    df_plot <- melt(df_plot)
    df_plot$variable <- rownames(df_plot)
    df_plot$variable <- gsub("^\\w", "",df_plot$variable)
    df_plot$variable <- as.numeric(df_plot$variable)
    df_plot$df_name <- df_name
    df_plotFinal <- rbind(df_plotFinal,df_plot)
  }
  ggplot(data = df_plotFinal) + 
    geom_line(mapping = aes(x = variable, y = value, color= df_name), size=1, alpha=0.8) +
    theme_classic(12)+ scale_color_brewer(palette="Dark2", direction = 1) +
    labs(y= "Mean quality score", x="Position in the read") +
    theme(legend.title = element_blank(), legend.position = "bottom")+
    scale_y_continuous(limits = c(0, 41))+
    guides(color=guide_legend(ncol=2))
}

# read files for each sequence run and each read orientation
2G00089_F <- read.table("/gs/gsfs0/users/cgazollavo/SILVER/16S/2G00089/reverse_qc/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt", header = T)
2G00089_R <- read.table("/gs/gsfs0/users/cgazollavo/SILVER/16S/2G00089/reverse_qc/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt", header = T)
2G00093_F <- read.table("/gs/gsfs0/users/cgazollavo/SILVER/16S/2G00089/reverse_qc/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt", header = T)
2G00093_R <- read.table("/gs/gsfs0/users/cgazollavo/SILVER/16S/2G00089/reverse_qc/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt", header = T)

multiqc_dfs= list("2G0089_16S"=QC89,"2G0093_16S"=QC93)

# run function to plot and then save it as a svg image
svg("quality_scores.svg", width=3.54331, height=3.54331)
plotQC(multiqc_dfs)
dev.off()

```


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
\

## Taxonomic classification

In SILVER, more than one sequence table was generated. Here we are going to merge them, as well as the track files.

On a fresh R session:

```{r}
# MERGE SEQUENCE TABLES AND TRACK IF MORE THAN 1 RUN

# track1 <- readRDS("path/to/run1/output/track_i.rds")
# track2 <- readRDS("path/to/run2/output/track_i.rds")
# track_ii <- rbind(track1, track2)

# seqtab.nochim1 <- readRDS("path/to/run1/output/seqtab.nochim1.rds")
# seqtab.nochim2 <- readRDS("path/to/run2/output/seqtab.nochim2.rds")
# seqtab.nochim <- mergeSequenceTables(seqtab.nochim1, seqtab.nochim2)


dna <- DNAStringSet(getSequences(seqtab)) # Create a DNAStringSet from the ASVs
load("/gs/gsfs0/users/cgazollavo/GTDB_r207-mod_April2022.RData")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=TRUE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab)

# save raw physeq
physeq_IDTAXA_GTDB<- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxid))
saveRDS(physeq_IDTAXA_GTDB, paste0("raw_physeq_IDTAXA_GTDB_", Sys.Date(),".rds"))
save.image(file='env_4.RData')

################################
```




## Requirements 

**NOTE:** This analysis requires at least 10GB of RAM to run.
It uses large files not included in the repository and many steps can take a few minutes to run. 

## Parameters

### Analysis input/output

```{r }
```
