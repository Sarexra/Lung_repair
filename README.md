# Lung_repair
Code and data for the manuscript "Tidal  ventilation induces divergent repair responses in lung epithelial cells and fibroblasts mediated by IL-6  signalling"
The folder contains the following files:

- "Lung_repair.R": R script used to load and process all the data files described below. All the analyses and figures shown in the manuscript are generated with this script.
  
- "Wound_healing_BEAS.xlsx": Sheet 1. Raw data used to quantify the wound healing assay in BEAS-2B cells supplemented with BALF from MV or CPAP patients.

- "Samples_BEAS": Sheet 1. Contains ColData (sample metadata) used for BEAS-2B RNA-seq analysis.
  
- "transcripts_BEAS_github.rds": Count matrix from BEAS-2B RNA-seq samples supplemented with BALF from MV or CPAP patients.
  STRING protein-protein interaction network was downloaded from the STRING website: https://version-11-5.string-db.org/cgi/download?sessionId=b4uePWoT3fVG
  
- "Wound_healing_MRC5.xlsx": Sheet 1. Raw data used to quantify the wound healing assay in MRC5 cells supplemented with BALF from MV or CPAP patients.

- "Samples_MRC5": Sheet 1. Contains ColData (sample metadata) used for MRC5 RNA-seq analysis.
  
- "transcripts_MRC5_github.rds": count matrix from MRC5 RNA-seq samples supplemented with BALF from MV or CPAP patients.
  STRING protein-protein interaction network was downloaded from the STRING website: https://version-11-5.string-db.org/cgi/download?sessionId=b4uePWoT3fVG
  
- "Wound_healing_BEAS.xlsx": Sheet 2. Raw data used to quantify the wound healing assay in BEAS-2B cells with BALF from MV or CPAP patients treated with tocilizumab.

- "Samples_BEAS": Sheet 2. Contains ColData (sample metadata) used for BEAS-2B with tocilizumab RNA-seq analysis.
  
- "transcripts_BEAS_toci_github.rds": Count matrix from BEAS-2B RNA-seq samples supplemented with BALF from MV or CPAP patients trated with tocilizumab.
  
- "Wound_healing_MRC5.xlsx": Sheet 2. Raw data used to quantify the wound healing assay in MRC5 cells with BALF from MV or CPAP patients treated with tocilizumab.

- "Samples_MRC5": Sheet 2. Contains ColData (sample metadata) used for MRC5 with tocilizumab RNA-seq analysis.
  
- "transcripts_MRC5_toci_github.rds": Count matrix from MRC5 RNA-seq samples supplemented with BALF from MV or CPAP patients treated with tocilizumab.
  
- "proteins_abundance_github.rds": Protein abundance matrix obtained from BALF proteomic analysis from patients under MV or CPAP.

- "human_lr_pair.rds":Database of human ligand–receptor pairs from CellTalkDB.
  STRING protein-protein interaction network was downloaded from the STRING website: https://version-11-5.string-db.org/cgi/download?sessionId=b4uePWoT3fVG
  
- "qPCR_BEAS_CsA": Raw qPCR data from IL-6 and PPIA relative expression in BEAS-2B cell line supplemented with BALF from MV or CPAP patients treated with cyclosporine-A (CsA).
