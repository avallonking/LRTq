{
    require(data.table)
    require(ggpubr)
    require(tidyr)
    require(dplyr)
    require(Cairo)
    setwd("./data")
    # load data - RV eGenes
    {
        indiv.with.rare.var <- lapply(dir("indiv.with.rare.var.idx/", full.names = T), function(x) {try(fread(x), T)})
        vt.rare.z <- lapply(dir("log2.standardized.corrected.tpm/", pattern = "vt", full.names = T), 
                            function(x) {try(fread(x), T)})
        skat.rare.z <- lapply(dir("log2.standardized.corrected.tpm/", pattern = "skat", full.names = T), 
                              function(x) {try(fread(x), T)})
        lrt.rare.z <- lapply(dir("log2.standardized.corrected.tpm/", pattern = "lrt", full.names = T), 
                             function(x) {try(fread(x), T)})
        # names(vt.rare.z) <- c("AVO", "AT", "EM", "Lung", "MS", "WB")
        # names(skat.rare.z) <- c("AVO", "AT", "EM", "Lung", "MS", "WB")
        # names(lrt.rare.z) <- c("AVO", "AT", "EM", "Lung", "MS", "WB")
        # names(indiv.with.rare.var) <- c("AVO", "AT", "EM", "Lung", "MS", "WB")
        names(indiv.with.rare.var) <- tstrsplit(dir("indiv.with.rare.var.idx/"), "\\.")[[7]]
        names(lrt.rare.z) <- tstrsplit(dir("log2.standardized.corrected.tpm/", pattern = "lrt"), "\\.")[[8]]
        names(vt.rare.z) <- tstrsplit(dir("log2.standardized.corrected.tpm/", pattern = "vt"), "\\.")[[8]]
        names(skat.rare.z) <- tstrsplit(dir("log2.standardized.corrected.tpm/", pattern = "skat"), "\\.")[[8]]
    }

    # number of outliers for each tissue
    {
        get.outliers.per.tissue <- function(tpm, rare.idx, threshold = 2) {
            result = data.frame()
            for (tissue in names(rare.idx)) {
                sample.tpm = tpm[[tissue]]
                if (!any(class(sample.tpm) == "data.frame")) {
                    result = rbind(result, rep(NA, ncol(result)))
                    next()
                }
                sample.rare = rare.idx[[tissue]]
                selected.gene.idx = match(sample.tpm$gene, sample.rare$gene)
                selected.gene = sample.tpm$gene
                sample.tpm = sample.tpm[, -1]
                sample.rare = sample.rare[selected.gene.idx, -1]
                sample.tpm = t(scale(t(sample.tpm)))

                outliers.idx = abs(sample.tpm) > threshold
                rare.outliers.idx = outliers.idx & sample.rare
                common.outliers.idx = outliers.idx & (!sample.rare)
                rare.controls.idx = (!outliers.idx) & sample.rare
                common.controls.idx = (!outliers.idx) & (!sample.rare)

                num.total.outliers = rowSums(outliers.idx)
                num.total.controls = rowSums(!outliers.idx)
                num.common.outliers = rowSums(common.outliers.idx)
                num.common.controls = rowSums(common.controls.idx)
                num.rare.outliers = rowSums(rare.outliers.idx)
                num.rare.controls = rowSums(rare.controls.idx)
                prop.rare.outliers = num.rare.outliers / num.total.outliers
                result = rbind(result, 
                               c(
                                 mean(num.total.outliers), 
                                 sd(num.total.outliers) / sqrt(length(num.total.outliers)),
                                 mean(num.total.controls),
                                 sd(num.total.controls) / sqrt(length(num.total.controls)),
                                 mean(num.rare.outliers), 
                                 sd(num.rare.outliers) / sqrt(length(num.rare.outliers)),
                                 mean(num.rare.controls),
                                 sd(num.rare.controls) / sqrt(length(num.rare.controls)),
                                 mean(num.common.outliers), 
                                 sd(num.common.outliers) / sqrt(length(num.common.outliers)),
                                 mean(num.common.controls),
                                 sd(num.common.controls) / sqrt(length(num.common.controls)),
                                 mean(prop.rare.outliers),
                                 sd(prop.rare.outliers) / sqrt(length(prop.rare.outliers)),
                                 sum(num.total.outliers, na.rm = T),
                                 sum(num.total.controls, na.rm = T),
                                 sum(num.common.outliers, na.rm = T),
                                 sum(num.common.controls, na.rm = T),
                                 sum(num.rare.outliers, na.rm = T),
                                 sum(num.rare.controls, na.rm = T),
                                 num.genes = length(num.total.outliers[!is.na(num.total.outliers)])
                                 )
                               )

            }
            colnames(result) = c(
                                 "mean.total.outliers", "se.total.outliers", 
                                 "mean.total.controls", "se.total.controls", 
                                 "mean.rare.outliers", "se.rare.outliers",
                                 "mean.rare.controls", "se.rare.controls",
                                 "mean.common.outliers", "se.common.outliers",
                                 "mean.common.controls", "se.common.controls",
                                 "mean.prop.rare.outliers", "se.prop.rare.outliers",
                                 "num.total.outliers", "num.total.controls",
                                 "num.common.outliers", "num.common.controls",
                                 "num.rare.outliers", "num.rare.controls", "num.genes"
                                 )
            result$tissue = names(rare.idx)
            return(result)
        }
        outliers.per.tissue <- vector("list", 3)
        names(outliers.per.tissue) <- c("vt", "skat", "lrt")
        outliers.per.tissue[["vt"]] <- get.outliers.per.tissue(vt.rare.z, indiv.with.rare.var)
        outliers.per.tissue[["skat"]] <- get.outliers.per.tissue(skat.rare.z, indiv.with.rare.var)
        outliers.per.tissue[["lrt"]] <- get.outliers.per.tissue(lrt.rare.z, indiv.with.rare.var)
        outliers.per.tissue[["vt"]]$method = "VT"
        outliers.per.tissue[["skat"]]$method = "SKAT-O"
        outliers.per.tissue[["lrt"]]$method = "LRT-q"
        result = rbindlist(outliers.per.tissue)

        result %>% filter(method == "LRT-q") -> lrt.q.outliers
        sort(lrt.q.outliers$mean.total.outliers)
        sort(lrt.q.outliers$mean.rare.outliers)
        sort(lrt.q.outliers$mean.prop.rare.outliers)
    }

    {
    CairoPDF("../../materials/outliers.prop.by.total.outliers.per.tissue.lrt.only.pdf", width = 12, height = 8)
    p <- ggplot(na.omit(result) %>% filter(method == "LRT-q") %>% arrange(mean.prop.rare.outliers) %>% 
                mutate(tissue=factor(tissue, levels=tissue)), 
            aes(x=tissue, y=mean.prop.rare.outliers, colour=tissue)) + 
            geom_pointrange(aes(x = tissue, ymin = mean.prop.rare.outliers + 1.96 * se.prop.rare.outliers, 
                    ymax = mean.prop.rare.outliers - 1.96 * se.prop.rare.outliers, 
                    colour = tissue), 
                position = position_dodge(width = 0.6)) +
            ylab("Proportion of outliers with rare variants") + theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size = 10),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 15), 
              legend.position = "none") + xlab("Tissue") + coord_flip()
            # ggtitle("Outliers carrying rare variants in RV eGenes")
            print(p)
            dev.off()
    }
}
