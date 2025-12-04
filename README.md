# Hearing Gene Prioritization Pipeline

## Overview
This project implements an end-to-end bioinformatics pipeline to prioritize therapeutic gene targets for hearing loss. It integrates multi-species single-cell RNA-seq data (human and mouse) across development, maturation, and aging stages to identify genes that are:
1.  Highly specific to hair cells.
2.  Conserved between human and mouse.
3.  Supported by GWAS evidence for hearing traits.
4.  Dynamically regulated during development or aging.

## Datasets
The pipeline utilizes the following public datasets:
*   **GSE213796**: Human fetal and adult inner ear snRNA-seq.
*   **GSE114157**: Mouse mature cochlear hair cells (Smart-seq2).
*   **GSE60019**: Mouse cochlear and utricle development (RNA-seq).
*   **GSE274279**: Mouse cochlea and utricle aging snRNA-seq.

## Therapeutic Priority Score (TPS)
The TPS is a composite metric calculated as:
```
TPS = 0.25 * HC_Specificity + 
      0.20 * Dev_Peak_Expression + 
      0.20 * Aging_Decline + 
      0.15 * Species_Conservation + 
      0.20 * GWAS_Support
```
Genes are ranked by TPS to suggest the most promising candidates for gene therapy or drug targeting.

## Installation & Usage

### Prerequisites
*   Python 3.10+
*   pip

### Setup
1.  Clone the repository:
    ```bash
    git clone https://github.com/yourusername/hearing-gene-prioritization.git
    cd hearing-gene-prioritization
    ```
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### Running the Pipeline
Execute the main script to download data (if configured) and run the full analysis:
```bash
python src/main.py
```
Results will be saved in the `results/` directory.

## Deployment
To deploy the results dashboard to GitHub Pages:
1.  Commit all changes, including the `results/` folder.
2.  Push to GitHub.
3.  Go to repository Settings > Pages.
4.  Select the `main` branch and root folder (or `docs` if configured).
5.  The site will be live at `https://yourusername.github.io/hearing-gene-prioritization/`.

## Citation
If you use this pipeline or results, please cite:
> [Author Name], et al. "Integrative Multi-Omics Prioritization of Hearing Loss Therapeutic Targets." (2025).

---
*Created by Antigravity AI*
