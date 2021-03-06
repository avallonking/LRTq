# egenes vs sample size
{
  library(data.table)
  library(ggpubr)
  library(Cairo)
  library(dplyr)
  
  egene.count = fread("./data/egene.counts.csv")
  sample.size = fread("./data/sample.size.csv")
  egene.count$sample.size = sample.size$total.samples
  egene.count$total.genes = sample.size$genes
  egene.count$lrt.egenes.per.gene = egene.count$lrt.egenes / egene.count$total.genes
  egene.count$lrt.distinct.egenes.per.gene = 
    egene.count$lrt.distinct.egenes / egene.count$total.genes
  egene.count$cv.egenes.per.gene = egene.count$gtex.egenes / egene.count$total.genes
  
  # Spplemental Fig S6A
  p = ggscatter(egene.count, x="sample.size", y="lrt.egenes.per.gene", color="tissue")
  CairoPDF("../materials/egenes.per.gene.vs.sample.size.pdf")
  ggpar(p, legend="none", ylab="Number of eGenes detected by LRT-q / total expressed genes",
        xlab="Number of samples", font.tickslab = 13, font.legend = 13, 
        font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  
  # Supplemental Fig S6B
  CairoPDF("../materials/novel.rv.egenes.each.tissue.pdf", 
           width = 12)
  t = melt(egene.count %>% filter(lrt.distinct.egenes > 0) %>% as.data.table(), 
           id.vars = c("tissue", "sample.size"), 
           measure.vars = c("vt.distinct.egenes", "lrt.distinct.egenes", 
                            "skat.distinct.egenes", "acat.distinct.egenes", 
                            "acat.o.distinct.egenes"))
  t$tissue.i = rep(c(1:length(unique(t$tissue))), 5)
  t$variable = c(rep("VT", length(unique(t$tissue))), 
                 rep("LRT-q", length(unique(t$tissue))), 
                 rep("SKAT-O", length(unique(t$tissue))), 
                 rep("ACAT-V", length(unique(t$tissue))), 
                 rep("ACAT-O", length(unique(t$tissue))))
  p = ggbarplot(t %>% filter(variable != "ACAT-V"), 
                x="tissue", y="value", fill="variable", position = position_dodge()) + rotate_x_text(45)
  ggpar(p, ylab="Number of novel RV eGenes", legend.title = "Method", legend = "bottom",
        xlab="Tissue", font.ytickslab = 16, font.xtickslab = 8, font.legend = 16, 
        font.main = 16, font.x = 16, font.y = 16)
  dev.off()
}

