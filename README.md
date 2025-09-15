# PVP DMAC Resource Repository  

This repository serves as a **centralized resource hub for the PVP Data Management and Analysis Core (DMAC)**.  
It is intended to collect and share commonly used scripts, pipelines, and documentation to support standardized and reproducible analysis workflows across projects.  

---

## 📂 Current Contents  

### Transcriptomics  

This folder contains scripts and pipelines related to transcriptomics (RNA-seq) data processing and analysis.  

- **Preprocessing**  
  - `rnaseq_pipeline.sh`: Shell script documenting the RNA-seq data processing pipeline currently used for the EPICHIPC project.  
    - Includes major steps for preprocessing, alignment, and downstream quantification.  

- **Analysis**  
  - `Transcriptomics_00_pca_functions.R`: Functions for performing PCA on transcriptomics data.  
  - `Transcriptomics_01_de_analysis_condition.R`: Script for differential expression (DE) analysis based on condition.  

---

## 🚀 Future Plans  

As the repository grows, it will include additional resource folders (e.g., **Metabolomics**, **Proteomics**, **Multi-omics Integration**) with standardized pipelines and functions to support PVP DMAC projects.  
