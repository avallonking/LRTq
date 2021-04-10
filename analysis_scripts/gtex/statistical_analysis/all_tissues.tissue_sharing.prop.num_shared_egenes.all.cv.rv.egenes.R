library(data.table)

args = commandArgs(trailingOnly=T)
folder = args[1]
q.thres = args[2]
num.egene.thres = as.numeric(args[3])

if (is.na(num.egene.thres)) {
    num.egene.thres = 100
}

# RV eGenes
method = 'lrt'
qvals = lapply(dir(folder, pattern=paste(method, ".q", sep=""), full.names=T), fread)
tissue.info = fread("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/result_summary/tissue.list", header=F)
tissues = tissue.info$V1
names(qvals) = tissues
egenes = lapply(qvals, function(x) x$Gene_ID[apply(x, 1, function(y) {any(as.numeric(y[-1]) < as.numeric(q.thres))})])
egenes = egenes[unlist(lapply(egenes, function(x) length(x) >= num.egene.thres))]

# count the number of RV eGenes in each tissue
rv.egene.tissue.count = lapply(egenes, length)

shared.matrix = as.data.frame(matrix(nrow=length(egenes), ncol=2))
colnames(shared.matrix) = c("num.tissues", "prop.shared")
egenes.shared.count = as.data.frame(table(unlist(egenes)))
for (i in 1:length(egenes)) {
    shared.matrix[i, 1] = i
#    if (i == 1) {
#        shared.matrix[i, 2] = 0
#    } else {
        shared.matrix[i, 2] = sum(egenes.shared.count$Freq == i) / nrow(egenes.shared.count)
#    }
}
saveRDS(shared.matrix, paste("tissue_sharing.prop_shared_egenes.by.num_tissues.lrt", "q.thres", q.thres, "num.egene.thres", num.egene.thres, "rv.egenes", "rds", sep="."))

# CV eGenes
load.data <- function(tissue) {
    filenames <- paste("/u/project/eeskin2/k8688933/meta-tissue/gtex_v8_eqtls/GTEx_Analysis_v8_eQTL/", tissue, ".v8.egenes.txt.gz", sep="")
    file <- fread(cmd = paste("zcat", filenames))
    return(file)
}

get.egenes <- function(qvals, num.egene) {
    qvals = qvals[order(qvals$qval), ]
    egenes = qvals$gene_id[1:num.egene]
    return(egenes)
}

tissues = names(egenes)
q.data = lapply(tissues, load.data)
names(q.data) = tissues
egenes = vector("list", length(tissues))
names(egenes) = tissues
for (i in 1:length(tissues)) {
    num.egene = rv.egene.tissue.count[[i]]
    egenes[[i]] = get.egenes(q.data[[i]], num.egene)
}

shared.matrix = as.data.frame(matrix(nrow=length(egenes), ncol=2))
colnames(shared.matrix) = c("num.tissues", "prop.shared")
egenes.shared.count = as.data.frame(table(unlist(egenes)))

for (i in 1:length(egenes)) {
    shared.matrix[i, 1] = i
#    if (i == 1) {
#        shared.matrix[i, 2] = 0
#    } else {
        shared.matrix[i, 2] = sum(egenes.shared.count$Freq == i) / nrow(egenes.shared.count)
#    }
}
saveRDS(shared.matrix, paste("tissue_sharing.prop_shared_egenes.by.num_tissues.lrt", "q.thres", q.thres, "num.egene.thres", num.egene.thres, "cv.egenes", "rds", sep="."))
