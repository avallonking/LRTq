{
  require(data.table)
  require(ggpubr)
  require(tidyr)
  require(dplyr)
  require(Cairo)
  setwd("/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/")
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
    
    compute.outliers.enrichment <- function(sample.tpm, sample.rare) {
      if (!any(class(sample.tpm) == "data.frame")) {
        return()
      }
      # not compute for each gene, but for the whole dataset
      result = data.frame()
      selected.gene.idx = match(sample.tpm$gene, sample.rare$gene)
      sample.tpm = sample.tpm[, -1]
      sample.rare = sample.rare[selected.gene.idx, -1]
      sample.tpm = t(scale(t(sample.tpm)))
      for (threshold in seq(1, 10)) {
        outliers.idx = abs(sample.tpm) > threshold
        rare.outliers.idx = outliers.idx & sample.rare
        common.outliers.idx = outliers.idx & (!sample.rare)
        rare.controls.idx = (!outliers.idx) & sample.rare
        common.controls.idx = (!outliers.idx) & (!sample.rare)
        num.total.outliers = sum(outliers.idx)
        num.total.controls = sum(!outliers.idx)
        num.common.outliers = sum(common.outliers.idx)
        num.common.controls = sum(common.controls.idx)
        num.rare.outliers = sum(rare.outliers.idx)
        num.rare.controls = sum(rare.controls.idx)
        relative.risk <- (num.rare.outliers / num.total.outliers) / (num.rare.controls / num.total.controls)
        log.se <- sqrt(1/num.rare.outliers - 1/num.total.outliers + 1/num.rare.controls - 1/num.total.controls)
        max.ci <- relative.risk * exp(1.96 * log.se)
        min.ci <- relative.risk * exp(-1.96 * log.se)
        pval <- pnorm(log(relative.risk) / log.se, lower.tail = F)
        result = rbind(result, c(num.total.outliers, num.total.controls, num.rare.outliers, num.rare.controls, 
                                 num.common.outliers, num.common.controls, relative.risk, log.se, max.ci, min.ci, 
                                 pval, threshold))
        
      }
      colnames(result) = c("total.outliers", "total.controls", "rare.outliers", "rare.controls", "common.outliers", 
                           "common.controls", "relative.risk", "log.se", "max.ci", "min.ci", "pval", "threshold")
      return(result)
    }
    
    # standard.z <- function(x) {
    #   return((x - mean(x)) / sd(x))
    # }
    outliers.enrichment = vector("list", length(indiv.with.rare.var))
    names(outliers.enrichment) <- names(indiv.with.rare.var)
    for (tissue in names(indiv.with.rare.var)) {
      outliers.enrichment[[tissue]] = compute.outliers.enrichment(lrt.rare.z[[tissue]], 
                                                                  indiv.with.rare.var[[tissue]])
    }
  }

    {
    summary.outliers.enrichment = Reduce("+", lapply(outliers.enrichment, function(x) {x[, 1:6]}))
    summary.outliers.enrichment$threshold = seq(1, 10)
    summary.outliers.enrichment$log.se = sqrt(1 / summary.outliers.enrichment$rare.outliers -
                                                1 / summary.outliers.enrichment$total.outliers +
                                                1 / summary.outliers.enrichment$rare.controls -
                                                1 / summary.outliers.enrichment$total.controls)
    summary.outliers.enrichment$relative.risk = 
      (summary.outliers.enrichment$rare.outliers / summary.outliers.enrichment$total.outliers) /
      (summary.outliers.enrichment$rare.controls / summary.outliers.enrichment$total.controls)
    summary.outliers.enrichment$max.ci = 
      summary.outliers.enrichment$relative.risk * exp(1.96 * summary.outliers.enrichment$log.se)
    summary.outliers.enrichment$min.ci = 
      summary.outliers.enrichment$relative.risk * exp(-1.96 * summary.outliers.enrichment$log.se)
    summary.outliers.enrichment$pval = 
      pnorm(log(summary.outliers.enrichment$relative.risk) / summary.outliers.enrichment$log.se, lower.tail = F)
    
    CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/log_tpm.egenes.enrichment.for.rare.var.pdf")
    temp = summary.outliers.enrichment
    temp = temp[1:8, ]
    temp$threshold = as.character(temp$threshold)
    temp$threshold = factor(temp$threshold, levels = temp$threshold)
    p = ggplot(temp, aes(x = threshold, y = relative.risk)) +
      theme_bw() + ylab('Enrichment of rare variants') + xlab("Z-score threshold") +
      geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
      geom_pointrange(aes(x = threshold, ymin = min.ci, ymax = max.ci, colour = threshold)) +
      # geom_text(aes(label = sprintf("%s (%s)", total.outliers, formatC(pval, digits = 2, format = "e"))), 
      #           hjust = 0.05, vjust = -1) +
      geom_text(aes(label = sprintf("%s", total.outliers)), 
                hjust = 0.05, vjust = -1) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), legend.position = "none",
            axis.text.x = element_text(size = 13), text = element_text(size = 15),
            axis.text.y = element_text(size = 13),
            axis.line = element_line(colour = "black")) + coord_flip() + ylim(0.9, 1.6)
    
    print(p)
    dev.off()
  }
}
