
title: "Gene Set Enrichment Analysis"
author: "Umar Faruk Saidu"
date: "2024-12-10"

genes.sig <- readRDS(file = "~/WGCNALASSO/Output/genes.sig.rds")

genes.sig <- as.data.frame(genes.sig)

ranked_gene_list <- genes.sig$log2FoldChange
names(ranked_gene_list) <- rownames(genes.sig)
ranked_gene_list <- sort(ranked_gene_list, decreasing = TRUE)

# Perform GSEA for diabetes DEGs

gsea <- gseGO(geneList = ranked_gene_list,
                       OrgDb = org.Rn.eg.db,
                       ont = "BP",
                       keyType = "SYMBOL",
                       pvalueCutoff = 0.05,
                       verbose = TRUE,
                       pAdjustMethod = "none")

gsea_table <- as.data.frame(gsea)
write.csv(gsea_table,file = "~/WGCNALASSO/Output/gsea_table.csv", 
          row.names = FALSE)
gsea <- saveRDS(gsea, file = "~/WGCNALASSO/Output/gsea.rds")
gsea <- readRDS(file = "~/WGCNALASSO/Output/gsea.rds")

# plot top 5 gsea results
gsea1 <- gseaplot(gsea, title = gsea$Description[1], geneSetID = 1)
gsea2 <- gseaplot(gsea, title = gsea$Description[2], geneSetID = 2)
gsea3 <- gseaplot(gsea, title = gsea$Description[3], geneSetID = 3)
gsea4 <- gseaplot(gsea, title = gsea$Description[4], geneSetID = 4)
gsea5 <- gseaplot(gsea, title = gsea$Description[5], geneSetID = 5)
gsea6 <- gseaplot(gsea, title = gsea$Description[6], geneSetID = 6)
gsea7 <- gseaplot(gsea, title = gsea$Description[7], geneSetID = 7)
gsea8 <- gseaplot(gsea, title = gsea$Description[8], geneSetID = 8)
gsea9 <- gseaplot(gsea, title = gsea$Description[9], geneSetID = 9)
gsea10 <- gseaplot(gsea, title = gsea$Description[10], geneSetID = 10)

# Save objects
png("~/WGCNALASSO/Output/gsea1.png", res =  300, width = 2000, height = 1600)
print(gsea1)
dev.off()

png("~/WGCNALASSO/Output/gsea2.png", res =  300, width = 2000, height = 1600)
print(gsea2)
dev.off()

png("~/WGCNALASSO/Output/gsea3.png", res =  300, width = 2000, height = 1600)
print(gsea3)
dev.off()

png("~/WGCNALASSO/Output/gsea4.png", res =  300, width = 2000, height = 1600)
print(gsea4)
dev.off()

png("~/WGCNALASSO/Output/gsea5.png", res =  300, width = 2000, height = 1600)
print(gsea5)

png("~/WGCNALASSO/Output/gsea6.png", res =  300, width = 2000, height = 1600)
print(gsea6)
dev.off()

png("~/WGCNALASSO/Output/gsea7.png", res =  300, width = 2000, height = 1600)
print(gsea7)
dev.off()

png("~/WGCNALASSO/Output/gsea8.png", res =  300, width = 2000, height = 1600)
print(gsea8)
dev.off()

png("~/WGCNALASSO/Output/gsea9.png", res =  300, width = 2000, height = 1600)
print(gsea9)
dev.off()

png("~/WGCNALASSO/Output/gsea10.png", res =  300, width = 2000, height = 1600)
print(gsea10)
dev.off()
