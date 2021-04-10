# this script is to see how many egenes are outliers enriching for rare variants
.libPaths("/u/project/eeskin2/k8688933/R/x86_64-pc-linux-gnu-library/3.5")
library(data.table)
library(RNOmni)

options(stringsAsFactors = FALSE)

remove.missing.impute <- function(geno.matrix) {
  # remove missing rate >= 0.05
  if (dim(geno.matrix)[1] == 0) {
      return(NULL)
  }
  good.var.idx <- rowSums(geno.matrix == -1) < 0.05 * (dim(geno.matrix)[2] - 1)
  if (sum(good.var.idx) == 0) {
      return(NULL)
  } else {
      geno.matrix <- geno.matrix[good.var.idx, ]
      geno.matrix[geno.matrix == -1] <- 0
      return(geno.matrix)
  }
}

args <- commandArgs(trailingOnly=T)

tissue <- args[1]

lrt.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/qvals/", tissue, ".lrt.q", sep="")
skat.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/qvals/", tissue, ".skat.q", sep="")
vt.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/qvals/", tissue, ".vt.q", sep="")
acat.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/qvals/", tissue, ".acat.q", sep="")
acat.o.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/qvals/", tissue, ".acat.o.q", sep="")
covar.file <- paste("/u/project/eeskin2/k8688933/rare_var/covariates/v8/", tissue, ".v8.covariates.transposed.txt", sep="")
gene.set.file <- paste("/u/project/eeskin2/k8688933/rare_var/gene_set/tss_20k_v8/chr", c(1:22), ".set", sep="")
expr.file <- "/u/project/eeskin2/k8688933/rare_var/expression_matrices/gene_tpm/v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
geno.file <- paste("/u/project/eeskin2/k8688933/rare_var/genotypes/v8/all_eur_samples_matrix_maf0.05/chr.", c(1:22), ".genotypes.matrix.tsv", sep="")
gene.list.file <- paste("/u/project/eeskin2/k8688933/rare_var/debug/gene_list_v8/", "/", tissue, "/chr", c(1:22), ".gene.list", sep="")
attribute.file <- "/u/project/eeskin2/k8688933/rare_var/expression_matrices/gene_tpm/v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
sample.match.file <- "/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/tissue.name.match.csv"

rare.individual.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/indiv.with.rare.var.idx.rv.egenes.only/indiv.with.rare.var.egenes.only", tissue, sep=".")
rare.indiv <- data.frame()

lrt.rare.z.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/log2.standardized.corrected.tpm.rv.egenes.only/log2.standardized.corrected.lrt.rv.egenes.tpm", tissue, sep=".")
skat.rare.z.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/log2.standardized.corrected.tpm.rv.egenes.only/log2.standardized.corrected.skat.rv.egenes.tpm", tissue, sep=".")
vt.rare.z.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/log2.standardized.corrected.tpm.rv.egenes.only/log2.standardized.corrected.vt.rv.egenes.tpm", tissue, sep=".")
acat.rare.z.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/log2.standardized.corrected.tpm.rv.egenes.only/log2.standardized.corrected.acat.rv.egenes.tpm", tissue, sep=".")
acat.o.rare.z.file <- paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/log2.standardized.corrected.tpm.rv.egenes.only/log2.standardized.corrected.acat.o.rv.egenes.tpm", tissue, sep=".")

lrt.rare.z <- data.frame()
skat.rare.z <- data.frame()
vt.rare.z <- data.frame()
acat.rare.z <- data.frame()
acat.o.rare.z <- data.frame()

lrt.q <- fread(lrt.file)
skat.q <- fread(skat.file)
vt.q <- fread(vt.file)
acat.q <- fread(acat.file)
acat.o.q <- fread(acat.o.file)

skat.q.egene <- skat.q$Gene_ID[apply(skat.q, 1, function(x) {any(as.numeric(x[-1]) < 0.05)})]
lrt.q.egene <- lrt.q$Gene_ID[apply(lrt.q, 1, function(x) {any(as.numeric(x[-1]) < 0.05)})]
vt.q.egene <- vt.q$Gene_ID[apply(vt.q, 1, function(x) {any(as.numeric(x[-1]) < 0.05)})]
acat.q.egene <- acat.q$Gene_ID[apply(acat.q, 1, function(x) {any(as.numeric(x[-1]) < 0.05)})]
acat.o.q.egene <- acat.o.q$Gene_ID[apply(acat.o.q, 1, function(x) {any(as.numeric(x[-1]) < 0.05)})]

gene.list <- lapply(gene.list.file, function(x) {fread(x, header=F, data.table=F)})
all.egenes <- union(lrt.q.egene, union(skat.q.egene, vt.q.egene))
all.egenes <- union(all.egenes, union(acat.q.egene, acat.o.q.egene))
filtered.gene.list <- lapply(gene.list, function(x) {intersect(x$V1, all.egenes)})


tissue.expr <- fread(expr.file, data.table=F)
# tissue.expr$Name <- gsub("\\..*", "", tissue.expr$Name)
sample.attr <- fread(attribute.file)
sample.match <- fread(sample.match.file)
tissue.query = sample.match$tissue.name[which(sample.match$tissue == tissue)]
sample.attr.tissue = sample.attr[which(sample.attr$SMTSD == tissue.query), ]
# sample.attr %>% filter(SMTSD == tissue.query) -> sample.attr.tissue

tissue.expr <- tissue.expr[, c("Name", intersect(sample.attr.tissue$SAMPID, colnames(tissue.expr)))]
new.sample.name <- paste(tstrsplit(x=colnames(tissue.expr)[-1], split="-", fixed=T, keep=1)[[1]], 
                         tstrsplit(x=colnames(tissue.expr)[-1], split="-", fixed=T, keep=2)[[1]],
                         sep="-")
colnames(tissue.expr) <- c("Name", new.sample.name)

# loop over 22 autosomes
for (i in 1:22) {
    gene.set <- fread(gene.set.file[[i]], header=F, data.table=F)
    geno.matrix <- fread(geno.file[[i]], data.table=F)
    covar <- fread(covar.file, data.table=F)

    selected.idv.idx <- which(colnames(tissue.expr) %in% colnames(geno.matrix))
    selected.idv.idx1 <- which(colnames(geno.matrix) %in% colnames(tissue.expr)[selected.idv.idx])
    geno.matrix <- geno.matrix[, c(1, selected.idv.idx1)]
    selected.idv.idx2 <- which(covar$ID %in% colnames(geno.matrix))
    covar <- covar[selected.idv.idx2, -1]
    sample.names <- colnames(geno.matrix)[-1]

    for (selected.gene in filtered.gene.list[[i]]) {
        print(paste("Gene: ", selected.gene))
        print(paste("Chromosome: ", i))

        included.snp <- gene.set$V2[gene.set$V1 == selected.gene]
        G <- geno.matrix[geno.matrix$ID %in% included.snp, ]
        G <- remove.missing.impute(G)

        selected.snp <- G[, 1]
        G <- G[, -1]
        G <- t(G)
        maf.list <- colMeans(G, na.rm = T) / 2
        G[, maf.list >= 0.50] = 2 - G[, maf.list >= 0.50]
        maf.list <- colMeans(G, na.rm = T) / 2

        rare.var.idx <- as.logical(maf.list > 0 & maf.list < 0.05)
        G <- as.matrix(G[, rare.var.idx])

        expr <- as.numeric(tissue.expr[tissue.expr$Name==selected.gene, selected.idv.idx])
        expr <- log2(expr + 1)
        expr <- (expr - mean(expr)) / sd(expr)
        E <- lm(expr ~ as.matrix(covar))$residuals

        outlier.idx <- abs(E) > 2
        individual.with.rare.var.idx <- rowSums(G) > 0
        rare.outlier.idx <- outlier.idx & individual.with.rare.var.idx

        rare.indiv <- rbind(rare.indiv, c(selected.gene, individual.with.rare.var.idx))

        if (selected.gene %in% skat.q.egene) {
            skat.rare.z <- rbind(skat.rare.z, c(selected.gene, E))
        }
        if (selected.gene %in% vt.q.egene) {
            vt.rare.z <- rbind(vt.rare.z, c(selected.gene, E))
        }
        if (selected.gene %in% lrt.q.egene) {
            lrt.rare.z <- rbind(lrt.rare.z, c(selected.gene, E))
        }
        if (selected.gene %in% acat.q.egene) {
            acat.rare.z <- rbind(acat.rare.z, c(selected.gene, E))
        }
        if (selected.gene %in% acat.o.q.egene) {
            acat.o.rare.z <- rbind(acat.o.rare.z, c(selected.gene, E))
        }

        print(warnings())
    }
}

if (nrow(skat.rare.z) > 0) {
    colnames(skat.rare.z) <- c("gene", sample.names)
}
if (nrow(vt.rare.z) > 0) {
    colnames(vt.rare.z) <- c("gene", sample.names)
}
if (nrow(lrt.rare.z) > 0) {
    colnames(lrt.rare.z) <- c("gene", sample.names)
}
if (nrow(acat.rare.z) > 0) {
    colnames(acat.rare.z) <- c("gene", sample.names)
}
if (nrow(acat.o.rare.z) > 0) {
    colnames(acat.o.rare.z) <- c("gene", sample.names)
}

write.csv(skat.rare.z, skat.rare.z.file, quote=F, row.names=F)
write.csv(vt.rare.z, vt.rare.z.file, quote=F, row.names=F)
write.csv(lrt.rare.z, lrt.rare.z.file, quote=F, row.names=F)
write.csv(acat.rare.z, acat.rare.z.file, quote=F, row.names=F)
write.csv(acat.o.rare.z, acat.o.rare.z.file, quote=F, row.names=F)

if (nrow(rare.indiv) > 0) {
    colnames(rare.indiv) <- c("gene", sample.names)
}
write.csv(rare.indiv, rare.individual.file, quote=F, row.names=F)

print(warnings())
