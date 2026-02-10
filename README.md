# Hearing Gene Prioritization Pipeline (TPS v2.0)

## Overview
This project implements the **Therapeutic Priority Score (TPS)**, an integrative multi-omics framework to prioritize gene therapy targets for sensorineural hearing loss. It synthesizes 15 datasets across human and mouse models (2023-2025) to identify targets that are:
1.  **Strictly Specific** to the inner ear (Safety Axis).
2.  **Developmentally Potent** for regeneration (Efficacy Axis).
3.  **Vulnerable** to aging and noise trauma.

**Key Innovation (v2.0):** Unlike traditional rankings, TPS explicitly **decouples** biological discovery from human genetic evidence (GWAS). GWAS is used only for independent validation, preventing circular reasoning.

## Datasets (Key Examples)
The pipeline integrates 15 validated datasets, including:
*   **Wang et al. (2024)**: Diseased Human Utricle Atlas (Regenerative signal).
*   **Kelley et al. (2023)**: Human Inner Ear Organoid Atlas.
*   **Xia et al. (2025)**: Mouse Cochlea Aging Atlas.
*   **GSE194089**: Spatial Tonotopy Gradients.
*   *...and 11 others spanning development and injury models.*

## Therapeutic Priority Score (TPS)
The score is calculated using a **Weighted Multi-Criteria Decision Analysis (MCDA)** framework structured along two clinical axes.

### 1. Safety Axis (40%)
*   **Human-Integrated Specificity ($S_{hc}$, 40%):** Penalizes off-target expression to minimize toxicity.

### 2. Efficacy Axis (60%)
*   **Developmental Peak ($D_{peak}$, 15%):** Regenerative potential (E14-P7).
*   **Global Vulnerability ($A_{vuln}$, 15%):** Genes lost during aging/trauma.
*   **Human Support ($S_{sup}$, 15%):** Organoid validation.
*   **Spatial Tonotopy ($T_{spatial}$, 15%):** Frequency map relevance.

### Scoring Formula
$$ TPS = 0.40 \cdot S_{hc} + 0.15 \cdot S_{sup} + 0.15 \cdot A_{vuln} + 0.15 \cdot D_{peak} + 0.15 \cdot T_{spatial} $$

### Top Candidates
*   **Rank #1: POU4F3**
*   **Rank #2: SOX2**
*   **Rank #5: ATOH1**
*   *Note: This revised ranking emphasizes combinatorial reprogramming factors over simple structural genes.*

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

