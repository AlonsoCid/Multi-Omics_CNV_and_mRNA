 # Multi-Omics Integration Pipeline: Breast Cancer CNA and Gene Expression
 
This project showcases a computational pipeline designed to integrate and analyze DNA copy number alterations (CNA) and gene expression data. The analysis focuses on breast cancer samples, demonstrating the application of advanced N-dimensional approaches to identify correlated biological patterns across different omics layers.

## 1. Data and structure

The data comes from [Pollack J. R. et al. (2002)](https://pubmed.ncbi.nlm.nih.gov/12297621/). It was generated using a two-color microarray containing 6,691 genes. Each sample was hybridized against the same reference (a normal female leukocyte sample), meaning all expression and CNA values are presented as ratios. The experiment includes 41 samples, 4 breast cancer cell lines, 28 clinical samples from Norway, and 9 clinical samples from Stanford. 

The core data is stored in the pollack.RData file, which contains two primary data frames:

- pollack: Contains information for all samples across all chromosomes.
- pollack.nox: Contains information for all samples across all chromosomes, except for the X chromosome to avoid gender bias.

The datasets contains the following variables:
- CLID: Original microarray ID.
- NAME: Information about the gene and its genomic location, separated by a pipe ("|").
- X1X-X5X: Copy number alterations of five baseline samples.
- BT474-T47D: Copy number alterations of the four breast cancer cell lines.
- NORWAY.7-NORWAY.112: Copy number alterations of the Norway sample.
- STANFORD.2-STANFORD.A: Copy number alterations of the Stanford sample.
- BT474.mRNA-T47D.mRNA: Gene expression profiles of the four breast cancer cell lines.
- NORWAY.7.mRNA-NORWAY.112.mRNA: Gene expression profiles of the Norway sample.
- STANFORD.2.mRNA-STANFORD.A.mRNA: Gene expression profiles of the Stanford sample.
- gene: Gene HGNC Symbol.
- chr: Chromosomal location.
- start: Starting position in the chromosome.
- end: Ending position in the chromosome.

## 2. Analysis and Results

### Data analysis using FactomineR R-package to the Pollack experiment

Before combining the data into a Multy Factor Analysis (MFA), lets explore the structure of each omics layer to make sure there is no outliers or clear barch effect. CNA indicate if cancer cells have change the structure of the DNA by deleting or amplifiying gene copies, meanwhile the expression analysis measure the mRNAs,the level of expression of each gene.

#### CNA PCA
![pca_plot](results/pca.png)

Although the PCA analysis reveals clear differences between the cell lines and the clinical samples, it shows a lot of overlapping between Stanford and Norway samples.

#### Expression PCA
![pca_exp_plot](results/pca_exp.png)

The expression PCA clearly shows differences between the cell lines and the rest of the samples (in dim 1) and also between Stanford and Norway samples (in dim 2).

#### MFA
Now lets perform a MFA. We can't simply merge the RNA and CNV tables and ran a standard PCA, the dataset with higher variance (or more variables) would completely dominate the results. The MFA() function prevents this by mathematically weighting each block (dividing each by its first principal eigenvalue) so that the RNA block and the CNV block have an equal say in the final geometric space.

![mfa_layer](results/mfa_layer.png)

The MFA Groups representation plot reveals how each omics block contributes to the overall diferences. Dimension 1 (13.87% variance) is driven equally by both the transcriptomic (RNA) and genomic (CNA) blocks, representing a shared biological signature. This first dimension also captures the variance introduced by the biological diferences between cell origin (cond). Dimension 2 (7.12% variance) is driven almost entirely by the CNV block, capturing structural genomic variance that is independent of gene expression.

This tables show the top 10 most relevant genes for each dimension of the FMA. For dimension 1 the genes come from the RNA expression block and for dimension 2 from the CNA block. This selection is based on the MFA results, since each block seems to have more effect on its respective dimension.

#### Genes dimension 1
| Gene | Chr | Start | Value |
| :--- | :--- | :--- | :--- |
| *ACTA2* | 10 | 88,935,074 | 0.9139 |
| *ANXA1* | 9 | 73,151,757 | 0.8836 |
| *C1R* | 12 | 7,080,209 | 0.9071 |
| *COL3A1* | 2 | 188,974,320 | 0.8747 |
| *COL4A2* | 13 | 110,305,812 | 0.8625 |
| *FCGR2B* | 1 | 161,663,147 | 0.8730 |
| *IGFBP7* | 4 | 57,030,773 | 0.8891 |
| *IL6* | 7 | 22,725,884 | 0.8865 |
| *LRP1* | 12 | 57,128,493 | 0.8644 |
| *LUM* | 12 | 91,102,629 | 0.8912 |

#### Genes dimension 2
| Gene | Chr | Start | Value |
| :--- | :--- | :--- | :--- |
| *ACTN3* | 11 | 66,546,395 | 0.6886 |
| *BLVRA* | 7 | 43,758,680 | 0.7626 |
| *CHGB* | 20 | 5,911,430 | 0.7274 |
| *ENTPD6* | 20 | 25,195,693 | 0.7943 |
| *HOOK1* | 1 | 59,814,786 | 0.7569 |
| *IGFBP1* | 7 | 45,888,357 | 0.7230 |
| *IGFBP3* | 7 | 45,912,245 | 0.6922 |
| *LFNG* | 7 | 2,512,529 | 0.7447 |
| *MADD* | 11 | 47,269,161 | 0.6878 |
| *SLC2A1* | 1 | 42,925,375 | 0.6726 |

The Value column represent the correlation coefficient between the gene's expression or CNA ratio and the dimension of the MFA, indicating its contribution to the variance.

Check [this circos plot](results/mfa_circos.pdf) to see a visual representation of the chromosome position and effect (the file is in PDF because uses a vector format and allows for greater resolution than a PNG or JPG).

### Apply the DIABLO method using mixOmics R-package
...
