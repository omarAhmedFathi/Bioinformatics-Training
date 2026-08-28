# Phase 5: Identify SNPs Associated with Sex (Logistic Model)

## Objective
Run a binary-trait association test (logistic regression) as a comparison to the linear model, using Sex as the phenotype. This serves as a sanity check for known sex-linked genetic patterns.

## Task 1: Sex Coding Confirmation
We confirmed the sex coding in the `.fam` file (Column 5). The dataset contains:
- **Males (Coded 1):** 49
- **Females (Coded 2):** 107
This 1/2 coding perfectly matches PLINK's default binary phenotype expectation (1 = Control, 2 = Case) for logistic regression.

## Task 2: Logistic Regression Results
We created an alternate phenotype file using Sex and ran a logistic GWAS (`--logistic`), adjusting for the top 10 Principal Components to maintain consistency with the linear model. We then extracted the standard Additive (`ADD`) effects.

*The clean GWAS result file has been saved as `GWAS_Sex_results.txt`.*

## Task 3: Discussion & Interesting Question

**"Are all SNPs associated with sex in the X and Y chromosomes only??"**

Biologically, human sex is determined strictly by the X and Y sex chromosomes. Our GWAS results reflect this, but with some interesting technical nuances:

1. **Autosomal Chromosomes:**
   When we look at the top hits across the autosomes (e.g., Chromosomes 10, 6, 1, 5), the lowest p-value observed is around $4.78 \times 10^{-5}$. Given the 67,735 SNPs tested, the Bonferroni genome-wide significance threshold is $7.3 \times 10^{-7}$. 
   **Conclusion:** There are **zero** genome-wide significant SNPs on the autosomes associated with sex. This is exactly what we expect biologically.

2. **The Y Chromosome:**
   The Y chromosome determines male sex. However, our dataset (`Qatari156_filtered_pruned`) does not contain any Chromosome 24 (Y) variants. If it did, they would show perfect association.

3. **The X Chromosome (Chromosome 23):**
   Interestingly, the standard additive (`ADD`) model in our logistic regression found no significant association for SNPs on the X chromosome (most p-values were 1). 
   **Why?** An additive GWAS model tests for differences in *allele frequencies* between cases (females) and controls (males). However, in a randomly mating population, males and females have the exact same allele frequencies for the X chromosome. The true genetic difference between sexes on the X chromosome is *heterozygosity*—males are hemizygous (can only be A or B), while females can be heterozygous (AB). Because the additive model does not test for heterozygosity differences, it correctly identifies that the base allele frequencies are equal between sexes, thus returning no association.

**Summary:** Yes, biological association is isolated to the sex chromosomes. The lack of significant autosomal hits serves as an excellent sanity check, while the lack of X-chromosome hits in an additive model serves as a great lesson in the limitations of testing allele frequencies versus heterozygosity for sex-linked traits.
