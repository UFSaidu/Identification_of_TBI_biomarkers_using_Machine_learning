
title: "Differential Gene Expression Analysis using DeSeq2"
author: "Umar Faruk Saidu"
date: "2024-11-30"

# Load count data and associated metadata
counts_data <- read.csv("~/WGCNALASSO/Output/combined_tbi.csv", row.names = 1)
head(counts_data)

# Read in metadata
meta_data <- read.csv("~/WGCNALASSO/Output/combined_trait_tbi.csv", 
                      row.names = 1)
head(meta_data)

# Ensure column names of counts_data matches with row names of meta_data and
# in the same order
all(colnames(counts_data) %in% rownames(meta_data))
all(colnames(counts_data) == rownames(meta_data))


# When we tried to create DESeqDataSet object, DESeq2 throws an error that 
# some values in assay are not integers. Happens due to some data processing
# Check to confirm
any(counts_data %% 1 != 0)

# The above returns True. Check which indexes are non-integers
which(counts_data %% 1 != 0, arr.ind = TRUE)

# The values are raw counts so we rounded them
counts_data <- round(counts_data)

# Now create DESeqData object
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = meta_data,
                              design = ~ Batch + Group)

# Check to ensure we keep genes with at least 15 reads total
rowSums(counts(dds)) >= 15

# set the factor level. This tells DESeq2 which level is the reference. "HC" is
# our control group, thus we'll be comparing our disease group; "TBI" to "HC"
dds$Group <- relevel(dds$Group, ref = "HC")

# Run DESeq
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# Set padj value less than 0.05
res_o.o5 <- results(dds, alpha = 0.05)
summary(res_o.o5)
res_o.o5 <- res_o.o5[order(res_o.o5$pvalue),]
genes.sig <- subset(res_o.o5, padj < 0.05)
genes.sig <- as.data.frame(genes.sig)

write.csv(genes.sig,file = "~/WGCNALASSO/Output/genes.sig.csv")
saveRDS(genes.sig, "~/WGCNALASSO/Output/genes.sig.rds")

# You can save the list of DEGs as vector if you want
id_degs <- rownames(genes.sig)
write.csv(id_degs,file = "~/WGCNALASSO/Output/id_degs.txt", row.names = FALSE)

##================Create DE genes Visualizations=======================##

# Select the top 30 DEGs ranked by p-adjusted values. i.e most significant
top_30_DEGs <- genes.sig[1:30, ]

# Subset expression data to only include the 30 DEGs
heatmap_data <- combined_tbi[rownames(top_30_DEGs), ]

# Create annotation data
anno_col <- data.frame(combined_trait_tbi[, 1])
colnames(anno_col) <- "Group"
rownames(anno_col) <- colnames(combined_tbi)

pheatmap(heatmap_data, scale = "row", annotation_col = anno_col,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         annotation_colors = list(Group=c(TBI = "black", HC = "orange")),
         cutree_cols = 2,
         cutree_rows = 2,
         main = "Expression and clustering of top 30 DE genes",
         legend_title = "Z-score",
         fontsize = 11, cellwidth = 35, cellheight = 10.25)

# Create a volcano plot

vol_degs <- genes.sig

vol_degs$Regulation <- ifelse(vol_degs$log2FoldChange > 0.5,"Upregulated", 
                      ifelse(vol_degs$log2FoldChange < -0.5,"Downregulated", 
                             "NS"))

# Select the top 5 upregulated and downregulated genes
top5_up <- vol_degs[vol_degs$Regulation == "Upregulated", ]
top5_up <- top5_up[order(-top5_up$log2FoldChange), ][1:5, ]

top5_down <- vol_degs[vol_degs$Regulation == "Downregulated", ]
top5_down <- top5_down[order(top5_down$log2FoldChange), ][1:5, ]

# combine top genes
top_genes <- rbind(top5_up, top5_down)

vol_degs$Label <- ifelse(rownames(vol_degs) %in% 
                         rownames(top_genes),
                         rownames(vol_degs), "")

# Create the volcano plot
volcano_plot <- ggplot(vol_degs, aes(x = log2FoldChange, y = -log10(padj),
                                     color = Regulation)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text(aes(label = Label), vjust = -0.5, size = 3, color ="black") +
  theme_minimal() +
  labs(title = "Volcano plot of DE genes") +
  theme(plot.title = element_text(hjust = 0.5))

print(volcano_plot)
