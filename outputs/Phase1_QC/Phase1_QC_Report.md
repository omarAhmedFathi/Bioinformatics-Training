# Phase 1: Quality Control (QC) of Genomics Data

## Objective
Clean the raw genotype data and characterize allele frequency and missingness before any downstream analysis.

## Summary of Dataset
- **Initial Samples**: 156 (founders only, as loaded from `.fam`)
- **Initial SNPs**: 67,735 (as loaded from `.bim`)
- **Min MAF**: 0.05118
- **Max MAF**: 0.50000
- **Total Genotyping Rate**: 99.88%

*Note: The dataset name (`Qatari156_filtered_pruned`) and the minimum MAF value indicate that the dataset has already been pre-filtered for Minor Allele Frequency >= 5%.*

## Histograms
The following plots have been generated and saved as PDF files in this directory:
- **MAF Distribution**: `MAF_histogram.pdf`
- **Missingness per SNP**: `SNP_missingness_histogram.pdf`
- **Missingness per Sample**: `Sample_missingness_histogram.pdf`

## QC Filtering Thresholds & Results
We applied the standard QC filters along with multiple different threshold combinations. Since the data was pre-filtered, only highly stringent thresholds caused further SNP removal.

| Thresholds | Reason for Threshold | Remaining SNPs | Remaining Samples |
| :--- | :--- | :--- | :--- |
| `--maf 0.05 --geno 0.05 --hwe 1e-6` | **Standard QC**: Removes rare variants, SNPs with >5% missingness, and SNPs deviating from Hardy-Weinberg Equilibrium. | 67,735 | 156 |
| `--maf 0.01 --geno 0.05 --hwe 1e-6` | **Relaxed MAF**: Allows slightly rarer variants. | 67,735 | 156 |
| `--maf 0.10 --geno 0.05 --hwe 1e-6` | **Stricter MAF**: Focuses only on common variants (>10%). | 51,129 | 156 |
| `--maf 0.05 --geno 0.01 --hwe 1e-6` | **Stricter Missingness**: Requires 99% call rate per SNP. | 67,735 | 156 |
| `--maf 0.05 --geno 0.05 --hwe 1e-5` | **Relaxed HWE**: Allows slightly more deviation from HWE. | 67,735 | 156 |

### Consistency Note
The number of SNPs remained exactly the same for most standard filters because the input dataset was already pruned and filtered for MAF > 5% and low missingness. We only observed a reduction when raising the MAF threshold to 10% (`0.10`), which removed ~16,606 SNPs.
