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

To reproduce the figures in the paper, please download the `analysis_scripts` folder and the "[LRTq-Data](https://github.com/avallonking/LRTq-Data)" repository. Decompose the "LRTq-Data" repository under `analysis_scripts/figures`, and rename it as `data/`. **Note that it is important to use the correct folder name (`data/`), otherwise the scripts cannot work.** For example, a user with a Linux or MasOS machine can do
```bash
# download the scripts
git clone https://github.com/avallonking/LRTq
# go to the directory with the scripts to generate figures
cd LRTq/analysis_scripts/figures
# download the required data from Google Drive with the provided link
# rename the folder as data
mv LRTq-Data data
# make the directory for storing the generated plots
mkdir ../materials
# run the scripts to generate figures and tables
```

To re-run the simulation study, users can use the R scripts in ```LRTq/analysis_scripts/simulation/```:
- ```power.simulate.new.R```: power simulation. It runs LRT-q and other methods on the simulated data assuming there are at least one rare variants regulating gene expression. Usage: ```Rscript power.simulate.new.R [simulated haplotypes] [repeats] [causal ratio] [a (constant)] [output file name]```The simulated haplotypes could be ```data/simulation/haplotype/len5k_110var/processed.sim.hap1.100var.tsv``` 
- ```typeIerror.simulation.new.R```: type I error simulation. It runs LRT-q and other methods on the simulated data assuming there are no rare variants affecting gene expression. Usage: ```Rscript typeIerror.simulation.new.R [simulated haplotyes] [permutations] [repeats] [output file name]```The simulated haplotypes could be ```data/simulation/haplotype/len5k_110var/processed.sim.hap1.100var.tsv```

To re-run the analysis of GTEx, users can use the R scripts in ```LRTq/analysis_scripts/gtex/association_tests/```:
- ```acat.R```: run ACAT on the GTEx dataset
- ```acat.regress_out_common_eqtls.R```: run ACAT on the GTEx dataset and regress out the effects of common eQTLs
- ```faster.gtex_power_test.modified.fixed.more_perm.speed_up.other_tissues.R```: run LRT-q, SKAT-O, and VT on the GTEx dataset
- ```faster.gtex_power_test.modified.fixed.more_perm.speed_up.other_tissues.maf01.R```: run LRT-q, SKAT-O, and VT on the GTEx dataset, only considering rare variants with MAF < 0.01
- ```faster.gtex_power_test.modified.fixed.more_perm.speed_up.other_tissues.regress_out_common_eqtls.R```: run LRT-q, SKAT-O, and VT on the GTEx dataset and regress out the effects of common eQTLs


#### Organization of this repository
- Source codes for `LRTq` R package are in `LRTq/src/LRTq.cpp`
- Scripts for generating all figures in the paper are in `LRTq/analysis_scripts/figures/`
- Parameter settings for simulating genotypes are in `LRTq/analysis_scripts/cosi`
- R scripts for running LRT-q and other methods on the GTEx dataset are in `LRTq/analysis_scripts/gtex/association_tests`
- R scripts for analyzing the results of the GTEx dataset are in `LRTq/analysis_scripts/gtex/statistical_analysis`
- R scripts for running the simulation experiments are in `LRTq/analysis_scripts/simulation`
