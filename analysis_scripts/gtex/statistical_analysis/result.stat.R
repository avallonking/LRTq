.libPaths("/u/project/eeskin2/k8688933/R/x86_64-pc-linux-gnu-library/3.5/")

library(data.table)
library(fdrtool)
library(ACAT)

args = commandArgs(trailingOnly = T)
tissue.list.file = args[1]
out.folder = args[2]

convert.q <- function(df) {
    for (j in 2:ncol(df)) {
        df[, j] = fdrtool(df[, j], statistic="pvalue", plot=F)$qval
    }
    return(df)
}

count.egenes <- function(df, gtex.egenes) {
    egenes <- df$Gene_ID[apply(df, 1, function(x) { any(as.numeric(x[-1]) < 0.05) })]
    return(list(num.egenes=length(egenes), num.distinct.egenes=sum(!egenes %in% gtex.egenes)))
}

stat.one.tissue <- function(folder, tissue) {
    vt = rbindlist(lapply(dir(folder, pattern="vt", full.name=T), fread))
    vt = na.omit(as.data.frame(vt))

    lrt = rbindlist(lapply(dir(folder, pattern="lrt", full.name=T), fread))
    lrt = na.omit(as.data.frame(lrt))

    skat = rbindlist(lapply(dir(folder, pattern="skat", full.name=T), fread))
    skat = na.omit(as.data.frame(skat))

    acat = rbindlist(lapply(dir(folder, pattern="acat", full.name=T), fread))
    acat = na.omit(as.data.frame(acat))

    idx = match(skat$Gene_ID, acat$Gene_ID)
    acat = acat[idx, ]

    acat.o = acat
    for (i in 1:nrow(acat)) {
        for (j in 2:ncol(acat)) {
            if (is.na(acat$Gene_ID[i])) {
                next
            }
            acat.o[i, j] = ACAT(c(vt[i, j], skat[i, j], acat[i, j]))
        }
    }

    egenes = fread(paste("zcat ", "/u/project/eeskin2/k8688933/meta-tissue/gtex_v8_eqtls/GTEx_Analysis_v8_eQTL/", tissue, ".v8.egenes.txt.gz", sep=""))
    egenes.name = egenes$gene_id[egenes$qval < 0.05]
    # egenes.name = gsub("\\..*", "", egenes.name)

    vt = convert.q(vt)
    lrt = convert.q(lrt)
    skat = convert.q(skat)
    acat = convert.q(na.omit(acat))
    acat.o = convert.q(na.omit(acat.o))

    fwrite(vt, paste(out.folder, paste(tissue, "vt.q", sep="."), sep="/"))
    fwrite(lrt, paste(out.folder, paste(tissue, "lrt.q", sep="."), sep="/"))
    fwrite(skat, paste(out.folder, paste(tissue, "skat.q", sep="."), sep="/"))
    fwrite(acat, paste(out.folder, paste(tissue, "acat.q", sep="."), sep="/"))
    fwrite(acat.o, paste(out.folder, paste(tissue, "acat.o.q", sep="."), sep="/"))

    vt.egenes.count = count.egenes(vt, egenes.name)
    lrt.egenes.count = count.egenes(lrt, egenes.name)
    skat.egenes.count = count.egenes(skat, egenes.name)
    acat.egenes.count = count.egenes(acat, egenes.name)
    acat.o.egenes.count = count.egenes(acat.o, egenes.name)

    stat.df = data.frame(tissue=tissue, gtex.egenes=length(egenes.name), 
                         vt.egenes=vt.egenes.count$num.egenes, vt.distinct.egenes=vt.egenes.count$num.distinct.egenes,
                         lrt.egenes=lrt.egenes.count$num.egenes, lrt.distinct.egenes=lrt.egenes.count$num.distinct.egenes,
                         skat.egenes=skat.egenes.count$num.egenes, skat.distinct.egenes=skat.egenes.count$num.distinct.egenes,
                         acat.egenes=acat.egenes.count$num.egenes, acat.distinct.egenes=acat.egenes.count$num.distinct.egenes,
                         acat.o.egenes=acat.o.egenes.count$num.egenes, acat.o.distinct.egenes=acat.o.egenes.count$num.distinct.egenes
                         )
    return(stat.df)
}

result.df = data.frame()
tissue.list = fread(tissue.list.file, header=F)
for (tissue in tissue.list$V1) {
    tissue.folder = paste("/u/project/eeskin2/k8688933/rare_var/results/tss_20k_v8/sungoohw", tissue, "modified_fixed/maf0.05/perm10k/all_weights", sep="/")
    result.df = rbind(result.df, stat.one.tissue(tissue.folder, tissue))
}
fwrite(result.df, paste(out.folder, "egene.counts.csv", sep="/"))
