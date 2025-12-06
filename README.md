# Hearing Gene Prioritization Pipeline

## Overview
This project implements an end-to-end bioinformatics pipeline to prioritize therapeutic gene targets for hearing loss. It integrates multi-species single-cell RNA-seq data (human and mouse) across development, maturation, and aging stages to identify genes that are:
1.  **Highly specific to hair cells** (Target specificity).
2.  **Conserved between human and mouse** (Translational potential).
3.  **Supported by GWAS evidence** for hearing traits (Clinical relevance).
4.  **Dynamically regulated** during development or aging (Regenerative potential).

## Datasets
The pipeline utilizes a comprehensive set of public datasets:
*   **[GSE213796](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213796)**: Human fetal and adult inner ear snRNA-seq.
*   **[GSE114157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114157)**: Mouse mature cochlear hair cells (Smart-seq2).
*   **[GSE60019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60019)**: Mouse cochlear and utricle development (Bulk RNA-seq).
*   **[GSE274279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274279)**: Mouse cochlea and utricle aging snRNA-seq (Xia et al., 2025).
*   **[GSE137299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137299)**: Mouse developmental scRNA-seq (E14-P7).
*   **[GSE165502](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165502)**: Mouse Spiral Ganglion Neuron (SGN) scRNA-seq.
*   **[GSE210215](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210215)**: Runx1 cKO SGN atlas.
*   **[GSE210216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210216)**: Complementary Runx1 SGN study.
*   **[GSE132925](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132925)**: Multi-age SGN transcriptome.

## Therapeutic Priority Score (TPS)
The TPS is a composite metric derived from four key biological modules:
```
TPS = 0.30 * HC_Regulatory_Module + 
      0.20 * SC_Plasticity_Module + 
      0.25 * SGN_Vulnerability_Module + 
      0.25 * Human_Transferability
```
*   **HC Regulatory Module:** Correlation with key drivers (*ATOH1*, *POU4F3*, *GFI1*).
*   **SC Plasticity Module:** Correlation with regeneration markers (*SOX2*, *HES1*, *LGR5*).
*   **SGN Vulnerability Module:** Specificity to spiral ganglion neurons and noise sensitivity.
*   **Human Transferability:** Conservation in human organoid atlas.

Genes are ranked by TPS to suggest the most promising candidates for gene therapy or drug targeting. Top candidates include **OTOF**, **SLC26A5**, and **MYO7A**.

## Advanced Analysis
Beyond the base TPS, the pipeline includes:
1.  **Cell-Type Specific Scoring:** Identifies targets unique to Hair Cells, Supporting Cells, SGNs, and Stria Vascularis.
2.  **Noise Injury Integration:** Incorporates data from noise-exposed cochleae to prioritize rescue targets.
3.  **Pathway Enrichment:** GO & KEGG analysis to reveal biological mechanisms (e.g., Sensory Perception, Ion Transport).
4.  **Druggability:** Integrates pLI scores and subcellular localization for clinical feasibility.

## Manuscript & Validation
This repository includes a `manuscript/` directory with auto-generated drafts for a scientific paper:
*   `abstract.md`, `methods.md`, `results.md`
*   `sanity_check_report.md`: Validates that known deafness genes (e.g., *OTOF*, *MYO7A*) are highly ranked.

## Installation & Usage

### Prerequisites
*   Python 3.10+
*   pip

### Setup
1.  Clone the repository:
    ```bash
    git clone https://github.com/omrozmn/hearing-gene-prioritization.git
    cd hearing-gene-prioritization
    ```
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### Running the Pipeline
The pipeline scripts are located in the `src` and `scripts` directories. To reproduce the analysis:
```bash
# Ensure data is downloaded first (see scripts/01_download_data.py)
python scripts/02_pipeline.py
python scripts/04_advanced_analysis.py
```
Results will be saved in the `results/` directory, including:
*   `therapeutic_priority_scores.csv`: The ranked list of genes.
*   `TPS_final.csv`: Final prioritized list with druggability and module info.
*   `figures/`: Generated UMAPs, heatmaps, and dotplots.

## Deployment
This project includes a web dashboard to visualize the results.
*   **Live Demo:** [https://omrozmn.github.io/hearing-gene-prioritization/](https://omrozmn.github.io/hearing-gene-prioritization/)

To deploy your own:
1.  Push the repository to GitHub.
2.  Go to Settings > Pages.
3.  Select `main` branch and `/` root folder.

## Citation
If you use this pipeline or results, please cite:
> "Integrative Multi-Omics Prioritization of Hearing Loss Therapeutic Targets." (2025).

---
*Created by Antigravity AI*

## 📊 Dataset Sources (Automatically Generated)

This project integrates multiple validated inner-ear transcriptomic datasets across development, species, neuronal subtypes, inflammation models, and human inner ear organoids.  
A full manifest is stored in **data/dataset_manifest.json**.

To add new datasets, update the manifest file and rebuild the analysis pipeline.

