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
sample_map <- data.frame(geo_accession = c("GSM2718317", "GSM2718318", 
                                           "GSM2718319", "GSM2718320", 
                                           "GSM2718321", "GSM2718322"),
                         sample_id = c("Sham1", "Sham2", "Sham3", "TBI1", 
                                       "TBI2", "TBI3"),
                         condition = rep(c("Sham", "TBI"), each = 3),
                         row.names = c("Sham1", "Sham2", "Sham3", 
                                       "TBI1", "TBI2", "TBI3"))

# Extract sample IDs from count matrix column names
sample_ids <- sapply(strsplit(colnames(counts), "_"), `[`, 1)

# Create comprehensive cell metadata
cell_metadata <- data.frame(cell_barcode = colnames(counts),
                            sample_id = sample_ids,
                            geo_accession = sample_map[sample_ids, "geo_accession"],
                            condition = sample_map[sample_ids, "condition"],
                            row.names = colnames(counts))

# Verify the mapping
head(cell_metadata)
saveRDS(counts, file = "~/WGCNALASSO/Output/counts_sparse_matrix.rds")
saveRDS(cell_metadata, file = "~/WGCNALASSO/Output/cell_metadata.rds")


# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = cell_metadata)


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

# Since we are validating a known set of genes, 
# it's good to first check if genes exist in your data. 
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
        group.by = "condition",
        ncol = 3)

# Get an idea of the correlation between the QC metrics before filtering
plot1 <- FeatureScatter(seurat_obj, 
                        feature1 = "nCount_RNA", feature2 = "percent.mt", 
                        group.by = "condition")
plot2 <- FeatureScatter(seurat_obj, 
                        feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        group.by = "condition")
plot1 + plot2

# Our data showed significantly high mitochondrial genes. 
# Some cells have high mitochondrial genes above 40%. 
# We'll filter those cells and only keep cells with %mito below 10%.

# Filter out low cells and mitochondrial genes
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

VlnPlot(seurat_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# Get an idea of the correlation between the QC metrics after filtering
plot1 <- FeatureScatter(seurat_obj,
                        feature1 = "nCount_RNA", feature2 = "percent.mt", 
                        group.by = "condition")
plot2 <- FeatureScatter(seurat_obj,
                        feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        group.by = "condition")
plot1 + plot2

#--------------------------------------------------
# Step 3: Follow Normal Seurat Standard Workflow
#--------------------------------------------------

# Normalize data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                            scale.factor = 10000)

# Find highly variable feature (genes)
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)

# Scale data to remove unwanted sources of variation in the data. 
# So that the cells have mean expression of 0 and variance of 1. 
seurat_obj <- ScaleData(seurat_obj)

# Run PCA to remove technical noise from the data
seurat_obj <- RunPCA(seurat_obj)

# Use elbow plot and DimHeatmap to get an ideal of the number of PCs to choose.
# Using Elbowplot(), we observed an elbow around PC15-20, suggesting that 
# the majority of true signal is captured in the first 20 principal components.
ElbowPlot(seurat_obj)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)

# We used 18 PCs to find Gene clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

# Save seurat obj
saveRDS(seurat_obj, file = "~/WGCNALASSO/Output/seurat_obj.rds")
seurat_obj <- readRDS(file = "~/WGCNALASSO/Output/seurat_obj.rds")

# Marker Gene Identification for all Clusters
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.5)

write.csv(all_markers, file = "~/WGCNALASSO/Output/all_markers.csv")

#Get top 5 markers for each cluster
top_markers <- all_markers %>% 
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(top_markers, file = "~/WGCNALASSO/Output/top_markers.csv")


#--------------------------------------------------
# Step 4: Cell Annotation
#--------------------------------------------------

# Load cell markers from Panglao database
panglaoDB_full_markers <- read.delim("~/panglaoDB/PanglaoDB_full_markers_27_Mar_2020 .tsv.gz", 
                                     sep = "\t")

# Filter using Brain and Immune system
panglao_brain_immune_markers <- panglaoDB_full_markers %>%
  filter(organ == "Brain" | organ == "Immune system")

# Merge cluster markers with PanglaoDB
annotation_matches <- top_markers %>%
  mutate(gene_upper = toupper(gene)) %>%  # Convert mouse genes to uppercase
  inner_join(panglao_brain_immune_markers, 
             by = c("gene_upper" = "official.gene.symbol")) %>%  # Match uppercase
  group_by(cluster, cell.type) %>%
  summarise(n_matches = n(), .groups = "drop") %>%
  arrange(cluster, desc(n_matches)) 

# Top candidate per cluster
best_annotations <- annotation_matches %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = n_matches, with_ties = FALSE)

# Note, best annotation returned no cell type for clusters 8 and 9. We used
# manual annotation using known hallmark genes to annotate these clusters. Also,
# cluster 4 was wrongly annotate as Eosinophils, after cross-checking with known
# hallmark genes, we manually re-annotate cluster 4 as Endothelial cells. Next,
# using well known marker genes, we're able to classify interneurons as either
# Excitatory or Inhibitory neurons. It's important to cross-check annotations
# with well know marker genes for each specific cell type.

new_ids <- c(
  "0" = "C0-Oligodendrocytes",
  "1" = "C1-Astrocytes",
  "2" = "C4-Excitatory Neurons",
  "3" = "C3-Microglia",
  "4" = "C4-Endothelial",
  "5" = "C5-Immature Neurons",
  "6" = "C6-Oligodendrocyte progenitor cells",
  "7" = "C7-Inhibitory Neurons",
  "8" = "C8-ChoroidPlexus Epithelial",
  "9" = "C9-Mural",
  "10" = "C10-Meningeal Fibroblast",
  "11" = "C11-Ependymal",
  "12" = "C12-Excitatory Neurons"
)

# Add cell_type column based on RNA_snn_res.0.4 clusters
seurat_obj@meta.data$cell_type <- new_ids[as.character(seurat_obj@meta.data$RNA_snn_res.0.4)]

# Add labels to Seurat object
seurat_obj <- RenameIdents(seurat_obj, new_ids)

# Save final seurat obj
saveRDS(seurat_obj, file = "~/WGCNALASSO/Output/seurat_obj_annotated.rds")
seurat_obj <- readRDS(file = "~/WGCNALASSO/Output/seurat_obj_annotated.rds")

#----------------------------------------------------------
# Cell type specific visualization of target genes
#----------------------------------------------------------

DimPlot(seurat_obj, label = TRUE, repel = TRUE)

FeaturePlot(seurat_obj, features = c("Icam1", "Tyrobp"))

VlnPlot(seurat_obj, 
        features = c("Icam1", "Tyrobp"),
        layer = "data",    
        pt.size = 0.2) 

DotPlot(seurat_obj, 
        features = c("Icam1", "Tyrobp"),
        cols = c("lightgrey", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()

