<div align="center">
  
  ![Bioinformatics Portfolio Banner](banner.jpg)
  
  <br/>
[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)]()
[![PLINK](https://img.shields.io/badge/PLINK-1.9%20%2F%202.0-blue?style=for-the-badge)]()
[![GCTA](https://img.shields.io/badge/GCTA-1.93-green?style=for-the-badge)]()
[![ggplot2](https://img.shields.io/badge/ggplot2-3.4.2-orange?style=for-the-badge)]()
[![Phases](https://img.shields.io/badge/Phases-16-purple?style=for-the-badge)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)

# 🧬 Bioinformatics Pipeline Portfolio & Tutorial

*A complete, 16-phase bioinformatics workflow starting from raw genomic data quality control to advanced metabolite-QTL (mQTL) mapping and pathway enrichment.*

</div>

---

> [!NOTE]
> I have structured this repository not just as a showcase of my code, but as a **tutorial** to help you understand the *biology* and *methodology* behind each step of a modern Genome-Wide Association Study (GWAS).

## 📑 Table of Contents
- [🧬 Pipeline Overview](#-pipeline-overview)
- [🖼️ Sample Outputs](#️-sample-outputs)
- [📚 Curriculum & Content](#-curriculum--content)
  - [Task 1: Population Genetics & GWAS](#task-1-population-genetics--gwas-pipeline)
  - [Task 2: Kinship, Heritability, and Metabolite GWAS](#task-2-kinship-heritability-and-metabolite-mixed-model-gwas)
- [🚀 Quick Start](#-quick-start)
- [👨‍💻 Author](#-author)

---

## 🧬 Pipeline Overview

The pipeline is split into two interconnected tasks. Task 1 focuses on population genetics and basic GWAS, while Task 2 incorporates complex trait analysis via mixed models.

```mermaid
graph TD
    classDef default fill:#f9f9f9,stroke:#333,stroke-width:2px;
    classDef highlight fill:#d4e6f1,stroke:#2980b9,stroke-width:2px,color:#2c3e50;

    A[Task 1: Population Genetics & GWAS]:::highlight -->|QC'd Genomic Data| B(Task 2: Metabolite Mixed-Model GWAS):::highlight
    
    subgraph Task 1
        QC[Phase 1: QC] --> PCA[Phase 2/3: PCA & Clustering]
        PCA --> GWAS[Phase 4/5: GWAS]
        GWAS --> Annot[Phase 6/7: Annotation & Enrichment]
    end
    
    subgraph Task 2
        MetabQC[Phase 11: Metabolite QC] --> Assoc[Phase 12: Metabolite-Diabetes Association]
        QC --> Kinship[Phase 8-10: LD, Kinship, Heritability]
        Kinship --> mQTL[Phase 15/16: mQTL Mapping & Discovery]
        Assoc --> mQTL
    end
    
    class QC,PCA,GWAS,Annot,MetabQC,Assoc,Kinship,mQTL default;
```

---

## 🖼️ Sample Outputs

<div align="center">
  <table>
    <tr>
      <td align="center">
        <img src="outputs/Phase2_PCA/Phase2_PCA_Scatter_Plot-1.png" alt="PCA Scatter Plot" width="400"/>
        <br/><em>Population Structure via PCA</em>
      </td>
      <td align="center">
        <img src="outputs/Phase6_Annotation/Phase6_Manhattan_PC1_Annotated-1.png" alt="Annotated Manhattan Plot" width="400"/>
        <br/><em>Annotated Manhattan Plot for PC1</em>
      </td>
    </tr>
    <tr>
      <td align="center">
        <img src="outputs/Phase14_Network/Phase14_Partial_Correlation_Network-1.png" alt="Metabolite Network" width="400"/>
        <br/><em>Metabolite Partial Correlation Network</em>
      </td>
      <td align="center">
        <img src="outputs/Phase15_mQTL/Phase15_Manhattan_Adjusted-1.png" alt="mQTL Manhattan Plot" width="400"/>
        <br/><em>Adjusted mQTL Manhattan Plot</em>
      </td>
    </tr>
  </table>
</div>

---

## 📚 Curriculum & Content

### [Task 1: Population Genetics & GWAS Pipeline](Task1_Population_Genetics.md)

| Phase | Name | Key Output | Status |
| :---: | :--- | :--- | :---: |
| **1** | Quality Control (QC) | Cleaned PLINK binaries | ✅ |
| **2** | Principal Component Analysis | PCA Eigenvectors | ✅ |
| **3** | Sub-population Clustering | Cluster Assignments | ✅ |
| **4** | Structure Association | Linear Model SumStats | ✅ |
| **5** | Sex Association | Logistic Model SumStats | ✅ |
| **6** | Functional Annotation | Gene/Variant Mapping | ✅ |
| **7** | Pathway Enrichment | Over-represented Pathways | ✅ |

> [!TIP]
> Task 1 builds the foundation. Understanding population structure (Phases 2-3) is critical to avoiding confounding in your final GWAS!

### [Task 2: Kinship, Heritability, and Metabolite Mixed-Model GWAS](Task2_Metabolite_GWAS.md)

| Phase | Name | Key Output | Status |
| :---: | :--- | :--- | :---: |
| **8** | Linkage Disequilibrium | LD Heatmaps | ✅ |
| **9** | Kinship Matrix | GRM Matrix | ✅ |
| **10** | Heritability Estimation | $h^2$ Estimates | ✅ |
| **11** | Metabolite Data QC | Normalized Metabolites | ✅ |
| **12** | Phenotype Association | Volcano Plots | ✅ |
| **13** | Feature Classification | Variable Importance | ✅ |
| **14** | Correlation Network | Network Graph | ✅ |
| **15** | mQTL Mapping | GWAS Summary Statistics | ✅ |
| **16** | Functional Discovery | Annotated mQTL Pathways | ✅ |

> [!IMPORTANT]
> In Task 2, we use a mixed-model approach incorporating the Genetic Relationship Matrix (GRM) from Phase 9 to properly account for cryptic relatedness when mapping mQTLs.

---

## 🚀 Quick Start

<details>
<summary><b>Click here to view setup instructions</b></summary>
<br>

To get started with running this pipeline locally, follow these steps:

1. **Clone the repository**
   ```bash
   git clone https://github.com/omar-ahmed/bioinformatics_portfolio.git
   cd bioinformatics_portfolio
   ```

2. **Install Dependencies**
   Ensure you have `PLINK 1.9`, `GCTA`, and `R` installed. Inside R, install required packages:
   ```R
   install.packages(c("tidyverse", "ggplot2", "sommer", "qqman"))
   BiocManager::install(c("clusterProfiler", "biomaRt"))
   ```

3. **Run Pipeline Scripts**
   Navigate to each phase directory to execute the bash and R scripts in order.

</details>

---

## 👨‍💻 Author

**Omar Ahmed**  
*Bioinformatics Developer & Researcher*

*Created as part of an intensive bioinformatics training internship.*
