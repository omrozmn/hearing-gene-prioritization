# 🎯 Therapeutic Priority Score (TPS) for Inner Ear

[![GitHub Pages](https://img.shields.io/badge/Website-Live-success)](https://omrozmn.github.io/hearing-gene-prioritization/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

> **Multi-module framework for prioritizing therapeutic gene targets in the human inner ear**

## 🌟 Overview

TPS (Therapeutic Priority Score) systematically prioritizes therapeutic gene targets for hearing loss by integrating:

- **41,697 human inner-ear cells** (100% real data, no downsampling)
- **2,297 hair cells** with real differential expression
- Multi-species data across development and aging
- GWAS associations and literature curation

**Key Results:**
- ✅ **OTOF** ranked #1 (synaptic function)
- ✅ **POU4F3** ranked #2 (HC transcription factor)  
- ✅ **GFI1** ranked #12 (HC differentiation)
- ✅ **Enrichment**: 182-213× fold (p < 10⁻¹⁰)

**🌐 Website:** [omrozmn.github.io/hearing-gene-prioritization](https://omrozmn.github.io/hearing-gene-prioritization/)

---

## 🧬 TPS Formula

```
TPS = 0.30 × HC_specificity 
    + 0.20 × SGN_specificity
    + 0.20 × Dev_peak
    + 0.15 × Aging_decline
    + 0.15 × GWAS
```

---

## 🏆 Top 10 Genes

| Rank | Gene | Score | Category |
|------|------|-------|----------|
| 1 | **OTOF** | 0.654 | Synaptic |
| 2 | **POU4F3** | 0.636 | TF |
| 3 | ESPNL | 0.607 | Stereocilia |
| 12 | **GFI1** | 0.513 | TF |
| 15 | **MYO7A** | 0.509 | Structural |
| 18 | **ATOH1** | 0.504 | Master TF |

**Full rankings:** [Download CSV](https://github.com/omrozmn/hearing-gene-prioritization/blob/main/docs/TPS_FINAL_100pct_REAL.csv)

---

## ✅ Validation

| Gene Set | Enrichment | P-value |
|----------|------------|---------|
| Hereditary Deafness | **182×** | 2.7×10⁻²¹ |
| HC Differentiation | **152×** | 1.1×10⁻¹⁰ |
| Sensory Sound | **213×** | 7.4×10⁻¹⁶ |

---

## 🚀 Quick Start

```bash
# Download final rankings
curl -O https://github.com/omrozmn/hearing-gene-prioritization/raw/main/docs/TPS_FINAL_100pct_REAL.csv

# Load in Python
import pandas as pd
tps = pd.read_csv('TPS_FINAL_100pct_REAL.csv', index_col=0)
print(tps.head(20))
```

---

## 📊 Datasets

- GSE213796: Human (41,697 cells)
- GSE114157: Mouse SGN
- GSE60019: Mouse Development
- GSE274279: Mouse Aging
- Literature: GWAS genes

---

## 🌐 Website Features

- Interactive tables
- Enrichment plots
- Methods documentation
- Download links

**Visit:** [omrozmn.github.io/hearing-gene-prioritization](https://omrozmn.github.io/hearing-gene-prioritization/)

---

## 📄 License

MIT License

---

**Created for the hearing research community** ❤️
