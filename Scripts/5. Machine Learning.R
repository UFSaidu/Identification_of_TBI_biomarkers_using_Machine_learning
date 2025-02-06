
title: "Machine Learning"
author: "Umar Faruk Saidu"
date: "2025-01-05"

#======Expression difference of the hub genes in Training data======#

# Subset the expression data and meta data to include only the 16 hub genes
gene_list <- true_hub$Gene
hub_genes_exp <- filtered_data[gene_list, ]
hub_genes_meta <- meta_data
hub_genes_exp <- as.data.frame(t(hub_genes_exp))
hub_genes_exp$Group <- hub_genes_meta$Group

# Reshape the data for ggplot2
hub_genes_melt <- melt(hub_genes_exp,
                       id.vars = "Group",
                       variable.name = "Gene",
                       value.name = "Expression")

# Create boxplot for each gene
plot_list <- list()

hub_genes_exp$Group <- as.factor(hub_genes_exp$Group)
comparisons <- list(c("HC", "TBI"))

for (gene in unique(hub_genes_melt$Gene)) {
  gene_data <- subset(hub_genes_melt, Gene == gene)
  
  p <- ggplot(gene_data, 
              aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    scale_fill_manual(values = c("TBI" = "red", "HC" = "blue")) +
    stat_compare_means(comparisons = comparisons, method = "t.test", 
                       label = "p.signif") +
    theme_minimal() +
    labs(title = paste(gene),
         x = "Group", y = "Expression") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = margin(t = 20),
          legend.position = "none",
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))
  
  plot_list[[gene]] <- p
    
}

# Use grid extra to display the plots
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 4))

ggsave("~/WGCNALASSO/Output/hub_genes_box_plot.png", 
       combined_plot, width = 26, 
       height = 15, background = "white", dpi = 600)

#==========Machine Learning Identification of core genes of TBI===========#

# We'll use DESeq2 variance-stabilized transformation (vst) data instead of
# the raw datasets. Vst data is normalized, corrected for sequencing depth, 
# and has same variance. This is more suitable for downstream analysis like
# MACHINE LEARNING. Our vst transformed training data is saved as filtered_data.

# Install packages
install.packages("glmnet")
library(glmnet)
install.packages("pROC")
library(pROC)
install.packages("caret")
library(caret)
install.packages("randomForest")
library(randomForest)
install.packages("e1071")
library(e1071)
install.packages("xgboost")
library(xgboost)

# Get vector of the 16 true hub genes
gene_list <- true_hub$Gene

#=================Perform LASSO Regression========================#
# Prepare data
X <- t(filtered_data[gene_list, ])
y <- as.numeric(combined_trait_tbi$Group == 'TBI')

# Fit LASSO
set.seed(123)
cv_fit <- cv.glmnet(X, y, alpha = 1, family = 'binomial')
lasso_fit <- glmnet(X, y, alpha = 1, lambda = cv_fit$lambda.min)

# Extract selected genes
selected_genes_lasso <- rownames(coef(lasso_fit))[which(coef(lasso_fit) != 0)][-1]

# Plot
plot(cv_fit)
plot(cv_fit$glmnet.fit, xvar = "lambda", label = TRUE)

#==================Perform Random Forest Algorithm======================#
# Train Random Forest
set.seed(123)
rf_fit <- randomForest(x = X, y = as.factor(y), 
                       importance = TRUE, 
                       ntree = 1000)

# Extract important genes
importance_scores <- importance(rf_fit)
selected_genes_RF <- rownames(importance_scores)[order(importance_scores[, "MeanDecreaseGini"],
                                                       decreasing = TRUE)][1:5]

# Plots
plot(rf_fit, main = "Random Forest")

top_5_genes <- importance_scores[selected_genes_RF, "MeanDecreaseGini"]

df <- data.frame(Gene = selected_genes_RF,
                 Importance = top_5_genes)

# Rank Genes based on importance
df <- df[order(df$Importance, decreasing = TRUE), ]
df$Gene <- factor(df$Gene, levels = df$Gene)

ggplot(df,
       aes(x = Importance, y = Gene)) +
  geom_segment(aes(xend = 0, yend = Gene), colour = "grey", size = 1.5) +
  geom_point(aes(fill = Importance), shape = 21, size = 5, stroke = 1.2) +
  scale_y_discrete(limits = rev(df$Gene)) +
  labs(x = "Importance",
       y = "Gene",
       fill = "Importance") +
  scale_fill_gradient(low = "skyblue", high = "red") +
  theme_bw() +
  theme(legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

