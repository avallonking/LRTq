# LRT-q
LRT-q is used for identifying regulatory effects of rare variants on genes with likelihood ratio test.

To install the R package:
```R
library(devtools)
install_github("avallonking/LRTq")
```

To see the manual:
```R
library(LRTq)
?LRTq
```

**Note that LRT-q performs the rare variant association test for one quantitative phenotype/trait (the expression levels for one gene). So users need to run LRT-q function for multiple times if there are different phenotypes or genes.**

Input: 
- Phenotypes **_E_**: a vector of quantitative traits of _N_ individuals, such as gene expression levels. It should be standardized and normalized.
- Genotypes **_G_**: an _N_ by _k_ matrix of individual genotypes, where _N_ represents the population size and _k_ stands for the number of rare variants. Rows are individuals, and columns are variants.
- Weights **_W_**: a vector of the weights of _k_ rare variants.
- Permutations **_perm_**: the number of permutations to perform to calculate the p-value.

Output:
- The p-value for this rare variant association test. If testing the association between gene expression and rare variants, then a significant p-value, such as a p-value smaller than 0.05, means that the gene expression is regulated by rare variants. Otherwise, there is no regulatory effect of rare variants on this gene.

Example:
```R
library(LRTq)

## use sample data provided by SKAT
require(SKAT)
data("SKAT.example")
attach(SKAT.example)

## use the quantitative trait (y.c), and extract the genotypes of rare variants (Z)
E = y.c # Phenotypes
maf = colMeans(Z) / 2
Z = Z[, maf > 0 & maf < 0.05]
G = Z # Genotypes
## weight all variants equally with 0.30
W = rep(0.30, ncol(G)) # Weights
## and use 1000 permutations
perm = 1000 # Permutations

## run LRT-q with the simulated inputs
LRTq(expr = E, geno = G, causal_ratio = rep(0.30, ncol(G)), perm = 1000)
## the results could be 0.000999001, 0.001998002, or 0.002997003, 
## due to the randomness in the permutation test
```

#### Organization of this repository
Scripts for generating all figures are in `LRTq/blob/master/analysis_scripts/figures/`
