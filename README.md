# **Traumatic Brain Injury (TBI) Biomarker Discovery**  

## **Project Progress**  
- Identified **16 key genes** from hippocampus bulk microarray data using **DESeq2** (DEGs), **WGCNA** (co-expression networks), **PPI** (Protein Network Analysis) and **GSEA** functional enrichment analysis. Genes were primarily inflammatory mediators and involved in inflammation and immune response.
- Performed feature selection using **Lasso, Random Forest (RF), and SVM-RFE** on the 16 key genes, identifying **2 core genes** as optimal TBI biomarkers.
- Investigatd **diagnostic potential** of core genes using **GLM, SVM, and XGBoost** with **ROC analysis**.
- Next, we constructed a **nomogram model** for risk prediction incorporating core genes.
- **Calibration & Decision Curve Analysis (DCA)** were performed to further assess clinical utility of core genes.
- **single-cell RNA-seq** validation of the core genes revealed predominant microglia expression, with minimal expression in astrocytes and endothelial cells. The core genes were involved in leukocyte recruitment and microglia activation, driving neuroinflammation and immune cell infiltration post-TBI.
- Performed **trajectory analysis** on microglia subtypes using **Slingshot** package in **R** to identify microglia lineages and map pseudotime injury responses.
- Further validated core genes using  microglia RNA-seq data from cortex, incorporating different time-points (acute, subacute, and chronic TBI stages). We observed that core genes exhibited temporal expression patterns across time-points.
- **Age** drive TBI progression and outcomes. Thus, immune cell infiltration driving by core genes were correlated with age and associated with injury progression and neurodegeneration.

## **Next Steps**  
- **Proteomic analysis** to quantify proteins encoded by the core genes in serum and cerebrospinal fluid.


📌 **This README will be updated as the project progresses. Full code will be uploaded to the repository after the paper has been accepted for publication.**  

Last updated: **September 2025**
