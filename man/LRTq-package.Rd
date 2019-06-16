\name{LRTq-package}
\alias{LRTq-package}
\alias{LRTq}
\docType{package}
\title{
  LRTq: group-wise rare variant association test for quantitative traits with likelihood ratio test
}
\description{
    LRTq, a powerful statistical method that infers causal statuses of rare variants with prior information and observational data and then aggregate rare variants in a nonlinear manner. With the assumption that individuals with abnormal gene expression levels are likely to enrich for causal rare variants, LRTq performs a likelihood ratio test to identify whether there are at least one causal rare variants in the variant group. It first computes the ratio between the null model that assumes no causal rare variants in a variant set and the alternative model that asserts at least one causal rare variants. LRTq then calculates p-value for the likelihood ratio between these two models using permutation test.

    This package is very easy to use. It only has one function, LRTq, which takes four variables: a vector of quantitative traits, genotypes matrix, a vector of causal probability and the number of permutation. The output would be a p-value for this association test.
}
\details{
    \tabular{ll}{
        Package: \tab LRTq\cr
        Type: \tab Package\cr
        Version: \tab 1.0\cr
        Date: \tab 2019-06-17\cr
        License: \tab MIT License\cr
    }
}
\author{
Jiajin Li (Albert Lee).

Maintainer: Jiajin Li (Albert Lee) <lijj36@ucla.edu>
}
\references{
  This optional section can contain literature or other references for
  background information.
}
\usage{
    LRTq(expr, geno, causal_ratio, perm = 1000)
}
\arguments{
    \item{expr}{A vector of quantitative traits of individuals. It should be standardized and normalized.}
    \item{geno}{A matrix of individual genotypes. Rows are individuals, and columns are variants.}
    \item{causal_ratio}{A vector with length equal to the number of variants. The probability of variants being causal, or simply, the weights of variants.}
    \item{perm}{An integer. The number of permutations.}
}
\examples{
    ## use sample data provided by SKAT
    require(SKAT)
    data("SKAT.example")
    attach(SKAT.example)
    ## use the quantitative trait (y.c), and extract the genotypes of rare variants (Z)
    E = y.c
    maf = colMeans(Z) / 2
    Z = Z[, maf > 0 & maf < 0.05]
    G = Z
    ## weight all variants equally with 0.30 and use 1000 permutations
    LRTq(expr = E, geno = G, causal_ratio = rep(0.30, ncol(G)), perm = 1000)
    ## the results could be 0.000999001, 0.001998002, 0.002997003, 
    ## 0.003996004, 0.004995005 or 0.005994006, due to randomness
}
