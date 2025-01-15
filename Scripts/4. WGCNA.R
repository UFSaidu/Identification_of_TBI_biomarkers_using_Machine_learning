
title: "WGCNA of TBI"
author: "Umar Faruk Saidu"
date: "2024-12-20"

##========Weighted Gene Co-expression Network Analysis (WGCNA)=====##

# Use vst transformation to normalize the raw data before using it for WGCNA

# The values are raw counts so we rounded them
counts_data <- round(counts_data)

# Create dds object without specifying a model
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = meta_data,
                              design = ~ 1)

# Perform variance stabilization
dds_norm <- vst(dds)
norm_counts <- assay(dds_norm)

##=======================Remove batch effects====================##

# We combined two different datasets, so is better to remove batch effects

# Remove batch effects using ComBat
mod_tbi <- model.matrix(~1, data = meta_data)
combat_tbi <- ComBat(dat = norm_counts, batch = batch_tbi, 
                          mod = mod_tbi)

# Perform WGCNA using the whole expression data. First, filter the data to 
# filter out non-varying genes based on the 75% median absolute deviation (MAD)

gene_MAD <- apply(combat_tbi, 1, mad)
mad_Threshold <- quantile(gene_MAD, 0.75)
filtered_data <- combat_tbi[gene_MAD > mad_Threshold, ]

# Ensure the data contain good samples and genes. i.e no any outliers
gsg = goodSamplesGenes(filtered_data, verbose = 3)
summary(gsg)
gsg$allOK

# If returns true, use the code bellow to confirm
# table(gsg$goodGenes)
# table(gsg$goodSamples)

# If False, use the code below to subset your data to only include non-outliers
# data <- data[gsg$goodGenes == TRUE, ] # remove genes that are outlier
# data <- data[, gsg$goodSamples == TRUE] # remove samples that are outlier

# Keep only genes and samples that are complete with no outliers
if (!gsg$allOK) {
  filtered_data = filtered_data[gsg$goodSamples, gsg$goodGenes]
}

# You can also check for outliers using hierarchical clustering
htree <- hclust(dist(t(filtered_data)), method = "average")

# Save plot
png("~/WGCNALASSO/Output/tree.png", width = 700)
plot(htree)
dev.off()


##=================Choose Power==========================##

# WGCNA requires a transpose data
filtered_data <- t(filtered_data)

# Choose a set of soft-threshold powers
power = c(c(4:10), seq(from = 12, to = 30, by = 2))

# Prepare soft-threshold for Scale Free Topology Model Fit
sft  <- pickSoftThreshold(filtered_data, 
                          powerVector = power,
                          verbose = 5)

sft_data <- sft$fitIndices

# Scale Independence plot
a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(title = "Scale Independence",
       x = "Soft Threshold (Power)", 
       y = "Scale Free Topology Model Fit, signed R^2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# Mean Connectivity plot
a2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(title = "Mean Connectivity",
       x = "Soft Threshold (Power)", 
       y = "Mean Connectivity") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# Arrange plots in an object
softThreshold_plot <- grid.arrange(a1, a2, nrow = 1)
print(softThreshold_plot)

# Save plot
ggsave("~/WGCNALASSO/Output/softThreshold_plot.png", 
       plot = softThreshold_plot, height = 6, width = 8)

##===========Construct Network using blockwiseModules function===============##

# convert data to matrix
filtered_data[] <- sapply(filtered_data, as.numeric)

# Set power based on the soft-threshold power plot
soft_power <- 14
temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseModules(filtered_data,
                          maxBlockSize = 8000,
                          TOMType = "unsigned",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1992,
                          verbose = 3)

# Plot dendrogram and module colors before and after merging
png("~/WGCNALASSO/Output/module_dendro.png", width = 700)
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, 
                    bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

dev.off()
cor <- temp_cor

# Save module labels and module eigengenes
moduleColors <- labels2colors(bwnet$colors)
moduleLabels <- bwnet$colors
moduleEigengenes <- bwnet$MEs

save(moduleColors, moduleLabels, moduleEigengenes, 
     file = "~/WGCNALASSO/Output/network.RData")

#==========================Correlate modules with traits================##

# Binarize categorical variables in trait data
# This is required for module-trait correlation
trait <- meta_data %>% 
  mutate(TBI = ifelse(grepl("TBI", Group), 1, 0)) %>% 
  mutate(HC = ifelse(grepl("HC", Group), 1, 0)) %>% 
  select(3,4)

# Correlate module eigengenes with traits
moduleTraitCor <- cor(moduleEigengenes, trait, 
                      method = "spearman", use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 
                                      nSamples = ncol(filtered_data))

# Visualize module-trait correlation as a heatmap
heatmap_data <- merge(moduleEigengenes, trait, by = 'row.names')
heatmap_data <- heatmap_data %>% column_to_rownames(var = 'Row.names')

png("~/WGCNALASSO/Output/moduleTraitCor_heatmap.png")
CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[9:10],
             y = names(heatmap_data)[1:8],
             col = c("blue", "skyblue", "white", "pink", "red"),
             main = "Module-trait correlation heatmap",
             cexMain = 1.2)
dev.off()

# Save correlation results
write.csv(moduleTraitCor, "~/WGCNALASSO/Output/moduleTraitCor.csv")
write.csv(moduleTraitPvalue, "~/WGCNALASSO/Output/moduleTraitPvalue.csv")

#================Extract genes from TBI significant modules================##

# Based on correlation, the green and brown modules are most significantly 
# correlated with TBI. Gene are saved in colors.
module_gene_mapping <- as.data.frame(bwnet$colors)

genes_green <- module_gene_mapping %>% 
  filter(bwnet$colors == 'green') %>% 
  rownames()

genes_brown <- module_gene_mapping %>% 
  filter(bwnet$colors == 'brown') %>% 
  rownames()

genes_tbi <- module_gene_mapping %>% 
  filter(bwnet$colors %in% c('green', 'brown')) %>% 
  rownames()

# Save genes to file for further analysis
write.csv(genes_green, "~/WGCNALASSO/Output/genes_green.csv", row.names = FALSE)
write.csv(genes_brown, "~/WGCNALASSO/Output/genes_brown.csv", row.names = FALSE)
write.csv(genes_tbi, "~/WGCNALASSO/Output/genes_tbi.csv", row.names = FALSE)

##==============Intramodular analysis: To identify driver genes=============##

# The module membership also known as intramodular connectivity is calculated 
# as the correlation of eigengene and gene expression data

# Get the top significance genes and associated p-values in TBI module
GS_cor <- cor(filtered_data, trait$TBI, method = "spearman", use = "p")
GS_cor_pvals <- corPvalueStudent(GS_cor, nSamples = ncol(filtered_data))

##===============Module-membership Vs Gene Significance=================##

moduleMembership <- cor(filtered_data, moduleEigengenes, 
                        method = "spearman", use = "p")

# Convert to data frame and assign column names
moduleMembership <- as.data.frame(moduleMembership)
colnames(moduleMembership) <- paste0("MM_", names(moduleEigengenes))

# Save module membership
write.csv(moduleMembership, "~/WGCNALASSO/Output/moduleMembership.csv", 
          row.names = TRUE)

# Extract module membership (MM) and gene significance (GS) for green module
MM_green <- moduleMembership[moduleColors == "green", "MM_MEgreen"]
GS_green <- GS_cor[moduleColors == "green", , drop = FALSE]

# Combine GS and MM for plotting
combine_green <- data.frame(GS = abs(GS_green[, 1]), 
                           MM = abs(MM_green),
                           Module = factor(rep("green", 
                                               length(MM_green))))

plot_green <- ggplot(combine_green, 
                        aes(x = MM, y = GS, 
                            colour = Module)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Module Membership (MM)",
       y = "Gene Significance (GS) for TBI",
       title = "Module Membership vs. Gene Significance\n
       cor = 0.73, p = 4.98e-05") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(plot_green)

# Save the plot
ggsave("~/WGCNALASSO/Output/plot_green.png", 
       plot = plot_green, height = 6, 
       width = 8)

MM_brown <- moduleMembership[moduleColors == "brown", "MM_MEbrown"]
GS_brown <- GS_cor[moduleColors == "brown", , drop = FALSE]

# Combine GS and MM for plotting
combine_brown <- data.frame(GS = abs(GS_brown[, 1]), 
                                MM = abs(MM_brown),
                                Module = factor(rep("brown", 
                                                    length(MM_brown))))

plot_brown <- ggplot(combine_brown, 
                         aes(x = MM, y = GS, 
                             colour = Module)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Module Membership (MM)",
       y = "Gene Significance (GS) for TBI",
       title = "Module Membership vs. Gene Significance\n
       cor = 0.82, p = 6.08e-07") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(plot_brown)

# Save the plot
ggsave("~/WGCNALASSO/Output/plot_brown.png", 
       plot = plot_brown, height = 6, 
       width = 8)

# Save important files into R for further analysis
save(genes_green, genes_brown, genes_tbi,
     file = "~/WGCNALASSO/Output/genes.RData")

# Brown module has the highest correlation and gene significance with TBI,
# Thus, get hub genes of brown module for further analysis.

WGCNA_hub <- genes_brown[which(MM_brown > 0.83 & GS_brown > 0.58)]
write.csv(WGCNA_hub, "~/WGCNALASSO/Output/WGCNA_hub.csv", row.names = FALSE)
