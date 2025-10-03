# **Traumatic Brain Injury (TBI) Biomarker Discovery**  

## **Project Progress**  
- Identified **16 key genes** from hippocampus bulk microarray data using **DESeq2** (DEGs), **WGCNA** (co-expression networks), and **PPI** (Protein Network Analysis). Genes were primarily inflammatory mediators and involved in inflammation and immune response.
- Performed feature selection using **Lasso, Random Forest (RF), and SVM-RFE** on the 16 key genes, identifying **2 core genes** as optimal TBI biomarkers.  
- **single-cell RNA-seq** validation of the core genes revealed predominant microglia expression, with minimal expression in astrocytes and endothelial cells. The core genes are involved in leukocyte recruitment and microglia activation, driving neuroinflammation and immune cell infiltration post-TBI.
- Performed **Trajectory analysis** on microglia subtypes using **Slingshot** package in **R** to identify microglia lineages and map pseudotime injury response.
- Further validation of core genes using  microglia RNA-seq data from cortex, incorporating different time-points (acute, subacute, and chronic TBI stages). Core genes exhibited temporal expression patterns across time-points.
- **Age** drive TBI progression and outcomes. Core genes were correlated with age and associated with injury progression and neurodegeneration.

## **Next Steps**  
- **Proteomic analysis** quantify proteins encoded by the core genes in serum and cerebrospinal fluid.
- **Diagnostic potential** of core genes using **GLM, SVM, and XGBoost** with **ROC analysis**.
- **Nomogram construction** for risk prediction.  
- **Calibration & Decision Curve Analysis (DCA)** to assess clinical utility of core genes.  

📌 **This README will be updated as the project progresses.**  

Last updated: **September 2025**
