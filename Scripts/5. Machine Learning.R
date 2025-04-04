
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
install.packages("rms")
library(rms)
install.packages("rmda")
library(rmda)

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

#==================Perform SVM-RFE Algorithm======================#
set.seed(123)
svm_fit <- rfe(X, as.factor(y), sizes = c(1:15),
               rfeControl = rfeControl(functions = caretFuncs, method = "cv"))

# Extract selected genes
selected_genes_svm <- predictors(svm_fit, )

# Plots
results <- svm_fit$results

max_accuracy <- max(results$Accuracy)
best_variables <- results$Variables[which.max(results$Accuracy)]

# Plot cross-validation accuracy
ggplot(results, 
       aes(x = Variables, y = Accuracy)) +
  geom_line(color = "blue") +
  geom_point(size = 2, color = "blue") +
  annotate("text", x = best_variables, y = max_accuracy,
           label = paste(round(max_accuracy, 2)),
           color = "red", hjust = 1.3, vjust = 0.3) +
  labs(x = "Number of Features",
       y = "Accuracy (Cross-Validation)") +
  theme(axix.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  theme_bw()

# Plot cross-validation Error
ggplot(results, 
       aes(x = Variables, y = 1 - Accuracy)) +
  geom_line(color = "blue") +
  geom_point(size = 2, color = "blue") +
  labs(x = "Number of Features",
       y = "Error (Cross-Validation)") +
  theme(axix.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  theme_bw()

#=====================Extract common core genes=====================#
# Get common genes present in all three algorithms
ml_gene_lists <- list(selected_genes_lasso,
                      selected_genes_RF,
                      selected_genes_svm)
# Extract core genes
core_genes <- Reduce(intersect, ml_gene_lists)

saveRDS(core_genes, file = "~/WGCNALASSO/Output/core_genes.rds")

# Create Venn diagram
venn_plot <- venn.diagram(x = list(
  LASSO = selected_genes_lasso,
  RandomForest = selected_genes_RF,
  `SVM-RFE` = selected_genes_svm
),
filename = "~/ML_AI/venn.png",
fill = c("blue", "red", "green"),
alpha = 0.5,
cex = 1.5,
cat.cex = 1.5,
cat.pos = c(-35, 25, 180))

#================Perform Nomogram and Decision Curve Analysis=================#
# Prepare data
X <- t(filtered_data[core_genes, ])
y <- as.numeric(combined_trait_tbi$Group == 'TBI')
data <- data.frame(X, Group = y)

# Fit logistic regression
ddist <- datadist(data)
options(datadist = "ddist")
logist_model <- lrm(Group ~ ., data = data)

# Create nomogram
nom <- nomogram(logist_model, fun = plogis, 
                fun.at = c(0.1, 0.5, 0.9),
                funlabel = "Disease Risk")
plot(nom)

# Calibration plot
cal <- calibrate(logist_model, method = "boot", B = 1000)

plot(cal, xlab = "Predicted Probability", legend = FALSE, 
     subtitles =  FALSE, cex.axis = 0.8,
     cex.lab = 0.8, ce.main = 0.9)
abline(0,1, col = "red", lty = 2)
legend("bottomright", 
       legend = c("Apparent", "Bias-corrected", "Ideal"),
       lty = c(3, 1, 2),
       col = c("black", "black", "red"),
       cex = 0.3,
       y.intersp = 0.8)

# Decision analysis
# Add probability to the data
data$pred_prob <- pred_prob
dca_data <- decision_curve(Group ~ pred_prob, data = data,
                           fitted.risk = TRUE, 
                           thresholds = seq(0, 1, by = 0.01))

# Plot
plot_decision_curve(dca_data, col = "black",
                    xlab = "Threshold Probability", 
                    ylab = "Net Benefit",
                    legend.position = "none",
                    curve.names = "two genes",
                    cost.benefit.axis = TRUE,
                    confidence.intervals = FALSE)
legend("topright",
       legend = c("two genes", "All", "None"),
       lwd = c(2, 1, 1),
       col = c("black", "grey", "black"),
       cex = 0.8,
       x.intersp = 0.2,
       y.intersp = 0.3,
       bty = "n")
