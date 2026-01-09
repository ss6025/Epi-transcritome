#!/usr/bin/Rscript
## load libraries
library(tidyverse)
library(reshape2)
options(stringsAsFactors = FALSE);
library(magrittr)
library(dplyr)
library(rstatix)
#install.packages("remotes")
library(remotes)
library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
r_lib_path <- args[2]

#remotes::install_github("okg3/MultiDataAnalysis")
#remotes::install_github("okg3/RNAEditingAnalysisTools")
 .libPaths( c( .libPaths(), r_lib_path) )
library(RNAEditingAnalysisTools)
library(data.table)


pdf_dir <- "PDF"

## Add / to directory path if it doesn't exist
output_dir <- ifelse(endsWith(output_dir, "/"), output_dir, paste0(output_dir, "/"))
pdf_dir <- ifelse(endsWith(pdf_dir, "/"), pdf_dir, paste0(pdf_dir, "/"))

anno_dir <- paste0(output_dir, "filtered_AI/annotated")
anno_dir <- ifelse(endsWith(anno_dir, "/"), anno_dir, paste0(anno_dir, "/"))

if (!dir.exists(pdf_dir)) {dir.create(pdf_dir)}


cat("Annotation directory: ",anno_dir)
cat("Output directory: ",output_dir)


#setwd("/work/greenbaum/users/chaconj1/REDItools_unique_bam_workflow/REAT_ready_files_PDAC/bedtools_intersected_files/")
#setwd(output_dir)
# set dirPath to directory containing files
# previous directory----


# Dbashata is split into separate files by sample
filelist<- list.files(anno_dir, pattern = ".txt")


######################################### Need to add annotations for each filtered REDItools out put file. 

## Create vector specifying columns of interest
# these elements represent unique editing events  
idCols <- c("Region", "Position", "Reference", "Strand")

# Create vector specifying columns used to compute editing index 
IndVars <- c("AllSubs", "Coverage-q25", "Frequency", "MeanQ", "BaseCount[A,C,G,T]")

## Combine these datasets (using the MergeIndividualFiles function)
RNAEdData <- MergeIndividualFiles(
  fileDirectory = anno_dir,
  filePattern = ".txt", 
  indVars = IndVars,
  IDcol = idCols,
  na.strings = "-"
)

## "Coverage-q30" could not be recoganzied by the GetEditedReads function. 
## Need to rename it to Coverage-q25, Either change it here or in the previous step for adding annotation

names(RNAEdData)


RNAEdData <- GetEditedReads(RNAEdData)
head(RNAEdData$EditedReads)

table(RNAEdData$Annotation$Strand)
## Global Editing Values by Sample----
globalEd <- GlobalEditing(RNAEdData)
globalEd
freq <- GlobalEditing(RNAEdData, type = "proportion") %>% data.frame() %>% rownames_to_column(var = "Sample_ID") 

colnames(freq)[2] <- "A-I editing"
  
write_tsv(freq, "GlobalEditing_Freq.csv")

######## Global editing by IR 
GlobalEditing(RNAEdData, type = "proportion", by = "IR")
freq <- GlobalEditing(RNAEdData, type = "proportion", by = "IR") %>% data.frame() %>% rownames_to_column(var = "IR_nonIR") 

write_tsv(freq, "GlobalEditing_Freq_IRvsNonIR.csv")

######## Global editing by rep_id 
GlobalEditing(RNAEdData, type = "proportion", by = "RepMask_gid")
freq <- GlobalEditing(RNAEdData, type = "proportion", by = "RepMask_gid") %>% data.frame() %>% rownames_to_column(var = "Sample_ID") 

write_tsv(freq, "GlobalEditing_Freq_bySubFamilies.csv")

###### substitution types - barplot
SubType <- GetSubType(RNAEdData$AllSubs)
subtypeCounts <-table(SubType)
subtypeCounts

# Donut plot 
subtypeCounts <- subtypeCounts[order(subtypeCounts, decreasing = TRUE)]
subtypeLabels <- names(subtypeCounts)
subtypeLabels[-grep("AG|CT|TC|GA", subtypeLabels)] <- NA

pdf(paste0(pdf_dir,"donut.pdf"))
donut(subtypeCounts, labels = subtypeLabels, outer.radius = 1,
      main = "RNA Editing Substitution Types - No Strand Correction")
dev.off()

## Alu regions
aluRepeats <- GetAluRepeat(RNAEdData$Annotation$RepMask_gid)
aluCounts <- table(aluRepeats)
aluCounts
pdf(paste0(pdf_dir,"donut_alus.pdf"))
donut(aluCounts, main = "RNA Editing Events in Alu and Non-Alu Repeats")
dev.off()

## freqency of editing
freqSummary <- SummarizeData(as.numeric(RNAEdData$Frequency), na.rm = TRUE)
round(freqSummary, 3)

pdf(paste0(pdf_dir,"RNA_editing_rate_histogram.pdf"))
par(fig= c(0,1,0.5,0.95))
boxplot(as.numeric(RNAEdData$Frequency), horizontal = TRUE, 
        axes = FALSE)

par(fig= c(0,1,0,0.8), new = TRUE)
hist(as.numeric(RNAEdData$Frequency),
     xlab = "Editing Frequency",
     ylab = "Count", 
     main = NULL,
     col = "grey")
mtext("Distribution of RNA Editing Frequency Rates", 
      side=3, outer=TRUE, line=-3) 
dev.off()
