  #----------------+++++++++++++++-------------------
  # Single-Cell RNA-seq Analysis of Icam1 and Tyrobp in TBI
  #-----------------+++++++++++++++-------------------
  # Install packages
  install.packages("Seurat")
  install.packages("tidyverse")
  install.packages("patchwork")
  install.packages("ggpubr") # For stat_compare_means()

  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(ggpubr)

  #--------------------------------------------------
  # Step 1: Prepare Data for Seurat object
  #--------------------------------------------------
  dge_data <- read.delim("~/WGCNALASSO/Input/GSE101901_DropSeqTBI.digital_expression.txt.gz", 
    header = TRUE, row.names = 1)
  counts <- as(as.matrix(dge_data), "sparseMatrix")

  # Note that steps for data manipulation depends on your data formatting and structure

  # The data we used is on GEO with accession number GSE101901.
  # The data is a DropSeq digital matrix with genes by sample_barcode. 
  # The columns are labeled as Treatment_barcode. 
  # The metadata contain sample characteristics to correctly map metadata to count matrix.

  # Create the sample mapping data frame
  sample_map <- data.frame(
    geo_accession = c("GSM2718317", "GSM2718318", "GSM2718319", 
                      "GSM2718320", "GSM2718321", "GSM2718322"),
    sample_id = c("Sham1", "Sham2", "Sham3", "TBI1", "TBI2", "TBI3"),
    condition = rep(c("Sham", "TBI"), each = 3),
    row.names = c("Sham1", "Sham2", "Sham3", "TBI1", "TBI2", "TBI3")
  )

  # Extract sample IDs from count matrix column names
  sample_ids <- sapply(strsplit(colnames(counts), "_"), `[`, 1)

  # Create comprehensive cell metadata
  cell_metadata <- data.frame(
    cell_barcode = colnames(counts),
    sample_id = sample_ids,
    geo_accession = sample_map[sample_ids, "geo_accession"],
    condition = sample_map[sample_ids, "condition"],
    row.names = colnames(counts)
  )

  # Verify the mapping
  head(cell_metadata)
  saveRDS(counts, file = "~/WGCNALASSO/Output/counts_sparse_matrix.rds")
  saveRDS(cell_metadata, file = "~/WGCNALASSO/Output/cell_metadata.rds")


  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = cell_metadata
  )


  # Add sample-level metadata to the Seurat object
  seurat_obj@meta.data$sample_id <- sample_ids
  seurat_obj@meta.data$geo_accession <- sample_map[sample_ids, "geo_accession"]
  seurat_obj@meta.data$condition <- sample_map[sample_ids, "condition"]

  # Check that the mapping worked correctly 
  table(seurat_obj$sample_id, seurat_obj$condition)

  # Should show:
  #        Sham TBI
  # Sham1    x   0
  # Sham2    x   0
  # Sham3    x   0
  # TBI1     0   x
  # TBI2     0   x
  # TBI3     0   x
  # (where x is the number of cells per sample)

# Since we are validating a known set of genes, it's good to first check if genes exist in your data. 
# This is important so as not to waste time on a data with no gene of interest.

# Check if genes are present
genes_to_check <- c("Icam1", "Tyrobp")
genes_found <- genes_to_check[genes_to_check %in% rownames(seurat_obj)]
if (length(genes_found) < length(genes_to_check)) {
  warning("Missing genes: ", setdiff(genes_to_check, genes_found))
}

print(genes_to_check)

#--------------------------------------------------
  # Step 2: Preprocessing and Quality Control
#--------------------------------------------------

# Calculate percentage of mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentFeatureSet(seurat_obj, 
  pattern = "^MT-|^mt-|^Mt-|^mito-|^Mito-")

VlnPlot(seurat_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "condition"
        ncol = 3)

# Get an idea of the correlation between the QC metrics before filtering
plot1 <- FeatureScatter(seurat_obj, 
  feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "condition")
plot2 <- FeatureScatter(seurat_obj, 
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "condition")
plot1 + plot2

# Our data showed significantly high mitochondrial genes. Some cells have high mitochondrial
# genes above 40%. We'll filter those cells and only keep cells with %mito below 10%.

# Filter out low cells and mitochondrial genes
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

VlnPlot(seurat_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# Get an idea of the correlation between the QC metrics after filtering
plot1 <- FeatureScatter(seurat_obj, 
  feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "condition")
plot2 <- FeatureScatter(seurat_obj, 
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "condition")
plot1 + plot2