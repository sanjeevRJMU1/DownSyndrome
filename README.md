# Companion Analysis Code – *Myocardial reprogramming by HMGN1 underlies heart defects in trisomy 21*

This repository provides the analysis scripts used in support of the manuscript:

> **Ranade SS\*, Li F\*, Whalen S, Pelonero A, Ye L, Huang Y, et al.**  
> *Myocardial reprogramming by HMGN1 underlies heart defects in trisomy 21.*  
> *Nature* (2025). https://doi.org/10.1038/s41586-025-09593-9  
> (\*equal contribution)

---

## Purpose

These scripts are shared to promote transparency and reproducibility. They represent generalized workflows used for data analysis in the published study.  

---

## Repository contents

### CUT&RUN (Extended Data Fig. 8)
- `cutNrun_nfcore_postprocessing_workflow.sh` 
- `cutNrun_create_regions_bed.R` 
- `cutNrun_stats.R` 

### scATAC-seq (Extended Data Fig. 8)
- `scATAC_workflow_DS.R`

### scRNA-seq (Figs. 1–4, Extended Data Figs. 1,2,4,6-9)
- `scRNA_generalized_workflow_Opus4.1.R` 
- `stats_lmer_or_wilcox.R` 
- `WGCNA_and_lmerStats.R`

---

## Data availability

All sequencing data are available through GEO under accession: [GSE271447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271447).  
Processed outputs for CUT&RUN, scRNA-seq, and scATAC-seq are also deposited.

---

## Code availability

All scripts here are shared under an open license for academic use.  
Please cite the above publication if you use or adapt any portion of this code.

---

## Contact

For questions, please contact:

- **Sanjeev S. Ranade** – [sranade@sbpdiscovery.org](mailto:sranade@sbpdiscovery.org)  
- **Angelo Pelonero** – [angelo.pelonero@gladstone.ucsf.edu](mailto:angelo.pelonero@gladstone.ucsf.edu)

---

## How to use this repository

1. Clone or download the repository from GitHub.  
   ```bash
   git clone https://github.com/<your-org>/<your-repo>.git
   cd <your-repo>
