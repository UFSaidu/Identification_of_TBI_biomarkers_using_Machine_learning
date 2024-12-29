
title: "Data input and manipulation"
author: "Umar Faruk Saidu"
date: "2024-11-10"

# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "WGCNA", "clusterProfiler", 
                       "enrichplot", "edgeR", "org.Rn.eg.db", 
                       "GEOquery", "DOSE", "sva", 
                       "biomaRT", "AnnotationDbi"))

install.packages("ggdendro")
install.packages("ggridges")
install.packages("pheatmap")
devtools::install_github("kevinblighe/CorLevelPlot")

# Load packages into working environment
library(WGCNA)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(edgeR)
library(GEOquery)
library(org.Rn.eg.db)
library(DOSE)
library(biomaRt)
library(readr)
library(tidyverse)
library(AnnotationDbi)
library(ggdendro)
library(ggridges)
library(RColorBrewer)
library(circlize)
library(pheatmap)
library(CorLevelPlot)

# Load first TBI dataset
tbi_count_data <- read_delim("~/WGCNALASSO/Input/GSE173975_TBI_GENE_MATRIX.txt", 
                               col_names = TRUE)

# Ensure no leading/trailing spaces. This is optional depending on your data
tbi_count_data$Gene <- trimws(tbi_count_data$Gene)

# Remove and identify duplicates
duplicates <- tbi_count_data$Gene[duplicated(tbi_count_data$Gene)]
if (length(duplicates) > 0) {
  print(paste("Found duplicates: ", paste(unique(duplicates), collapse = ", ")))
} else {
  print("No duplicates found")
}

# Remove duplicates if present
tbi_count_data <- tbi_count_data[!duplicated(tbi_count_data$Gene),]

# Change tibble to data frame
tbi_count_data <- as.data.frame(tbi_count_data)

# Rearrange column names to match sample names in metadata
tbi_data1 <- tbi_count_data %>% 
  select(c("Gene", "Sham.1.rep1", "Sham.1.rep2", "Sham.14.rep1", "Sham.14.rep2", 
           "Sham.14.rep3", "Sham.14.rep4", "TBI.1.rep1", "TBI.1.rep2", 
           "TBI.1.rep3", "TBI.1.rep4", "TBI.14.rep1", "TBI.14.rep2", 
           "TBI.14.rep3", "TBI.14.rep4"))

# Use geo_accession numbers as sample names
# Rename columns with geo_accession numbers 
# Ensure the order of columns are maintained as in metadata
colnames(tbi_data1) <- c("Gene", "GSM5283787", "GSM5283788", "GSM5283789", 
                         "GSM5283790", "GSM5283791", "GSM5283792", "GSM5283793", 
                         "GSM5283794", "GSM5283795", "GSM5283796", "GSM5283797", 
                         "GSM5283798", "GSM5283799", "GSM5283800")

# Set Gene column as row names
rownames(tbi_data1) <- tbi_data1$Gene
tbi_data1$Gene <- NULL

# Read metadata for first TBI dataset
meta_datatbi1 <- read.csv("~/WGCNALASSO/Input/metadata_tbi1.csv")

# Create trait data for first TBI dataset
trait_TBI1 <- meta_datatbi1 %>% 
  select(2,4) %>% 
  rename("Sample" = "geo_accession")

# Set sample column as row names
rownames(trait_TBI1) <- trait_TBI1$Sample

# Save to file
write.csv(tbi_data1, "~/WGCNALASSO/Output/tbi_data1.csv", row.names = TRUE)
write.csv(trait_TBI1, "~/WGCNALASSO/Output/trait_TBI1.csv", row.names = FALSE)

# Load second TBI dataset
tbi_count_data2 <- read_delim("~/WGCNALASSO/Input/GSE80174_Hippocampus_DEseq2_raw_counts.txt", 
                              col_names = TRUE)

# Rearrange column names to match sample names in metadata
tbi_data2 <- tbi_count_data2 %>% 
  select(c("Gene" = "Ensembl_Gene_ID", "RatID_1_mRNA", 
           "Rat2_sham_hippo_mrna_FC120_3.bam", 
           "Rat3_sham_hippo_mrna_FC120_5.bam", 
           "Rat4_sham_hippo_mrna_FC120_7.bam", 
           "Rat5_sham_hippo_mrna_FC122_1.bam", 
           "Rat6_trauma_hippo_mrna_FC120_2.bam", 
           "Rat7_trauma_hippo_mrna_FC120_4.bam", 
           "Rat8_trauma_hippo_mrna_FC120_6.bam", 
           "Rat9_trauma_hippo_mrna_FC120_8.bam", 
           "Rat10_trauma_hippo_mrna_FC122_2.bam"))

# Use geo_accession numbers as sample names
# Rename columns with geo_accession numbers
# Ensure the order of columns are maintained as in metadata
colnames(tbi_data2) <- c("Gene", "GSM2114177", "GSM2114178", "GSM2114179", 
                         "GSM2114180", "GSM2114181", "GSM2114182", "GSM2114183", 
                         "GSM2114184", "GSM2114185", "GSM2114186")

# Convert Ensembl IDS to Gene symbols
tbi_data2$Gene <- mapIds(org.Rn.eg.db, 
                         keys = tbi_data2$Gene,
                         keytype = "ENSEMBL",
                         column = "SYMBOL",
                         multiVals = "first")

# Remove NA values
tbi_data2 <- tbi_data2[!is.na(tbi_data2$Gene), ]

# Aggregate Gene column to handle duplicate errors
tbi_data2 <- tbi_data2 %>% 
  group_by(Gene) %>% 
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# Convert to dataframe
tbi_data2 <- as.data.frame(tbi_data2)

# Set Gene column as row names
rownames(tbi_data2) <- tbi_data2$Gene

# Remove Gene column now that it is row names
tbi_data2$Gene <- NULL

# Load metadata for second TBI dataset
metadata_tbi2 <- read_csv("~/WGCNALASSO/Input/metadata_tbi2.csv",
                          col_names = TRUE)

# Create trait data for second TBI dataset
trait_TBI2 <- metadata_tbi2 %>% 
  select(2,4)

trait_TBI2 <- as.data.frame(trait_TBI2)
rownames(trait_TBI2) <- trait_TBI2$Sample

# Save to file
write.csv(tbi_data2, "~/WGCNALASSO/Output/tbi_data2.csv", row.names = TRUE)
write.csv(trait_TBI2, "~/WGCNALASSO/Output/trait_TBI2.csv", row.names = FALSE)

##======================Filter & Merge dataset=========================##

# Include genes with expression values greater than 15 across whole samples
tbi_data1 <- tbi_data1 %>% 
  filter(rowSums(. > 15) == ncol(.))

tbi_data2 <- tbi_data2 %>% 
  filter(rowSums(. > 15) == ncol(.))

# Find common genes and then merge the datasets
common_genes_tbi <- intersect(rownames(tbi_data1), rownames(tbi_data2))

tbi_data1 <- tbi_data1[common_genes_tbi, ]
tbi_data2 <- tbi_data2[common_genes_tbi, ]

# Combine datasets
combined_tbi <- cbind(tbi_data1, tbi_data2)

write.csv(combined_tbi, "~/WGCNALASSO/Output/combined_tbi.csv", 
          row.names = TRUE)
saveRDS(combined_tbi, "~/WGCNALASSO/Output/combined_tbi.rds")


# Load third TBI dataset
tbi_count_data3 <- read.csv()


# Note, we'll use the combined data for differential gene expression analysis
# to identify the DEGs and also for the WGCNA analysis.
# However, the third dataset is for validation of our findings.