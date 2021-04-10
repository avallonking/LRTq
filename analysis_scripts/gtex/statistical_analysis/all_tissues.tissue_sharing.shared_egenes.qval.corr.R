library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=T)
folder = args[1]
threshold = args[2]

method = 'lrt'
qvals = lapply(dir(folder, pattern=paste(method, ".q", sep=""), full.names=T), fread)
tissue.info = fread("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/tissue.list", header=F)
tissues = tissue.info$V1
names(qvals) = tissues
egenes = lapply(qvals, function(x) x$Gene_ID[apply(x, 1, function(y) {any(as.numeric(y[-1]) < as.numeric(threshold))})])
corr.matrix = as.data.frame(matrix(nrow=length(tissues), ncol=length(tissues)))
colnames(corr.matrix) = tissues
rownames(corr.matrix) = tissues
for (i in 1:length(tissues)) {
    for (j in 1:length(tissues)) {
        shared.egenes = intersect(egenes[[i]], egenes[[j]])
        first.qval = qvals[[i]] %>% filter(Gene_ID %in% shared.egenes) %>% dplyr::select(-1) %>% as.list() %>% unlist() %>% as.numeric()
        second.qval = qvals[[j]] %>% filter(Gene_ID %in% shared.egenes) %>% dplyr::select(-1) %>% as.list() %>% unlist() %>% as.numeric()
        corr.matrix[i, j] = cor(first.qval, second.qval, method="pearson")
    }
}
saveRDS(corr.matrix, paste("tissue_sharing.shared_egenes.qval_corr", method, threshold, "rds", sep="."))
