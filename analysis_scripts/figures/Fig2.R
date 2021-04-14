{
  library(data.table)
  library(ggpubr)
  library(Cairo)
  library(dplyr)
  
  egene.count = fread("projects/rare_variant/gtex/all_tissues/gtex8/egene.counts.csv")
  sample.size = fread("projects/rare_variant/gtex/all_tissues/gtex8/sample.size.csv")
  egene.count$sample.size = sample.size$total.samples
  egene.count$total.genes = sample.size$genes
  egene.count$lrt.egenes.per.gene = egene.count$lrt.egenes / egene.count$total.genes
  egene.count$lrt.distinct.egenes.per.gene = 
    egene.count$lrt.distinct.egenes / egene.count$total.genes
  egene.count$cv.egenes.per.gene = egene.count$gtex.egenes / egene.count$total.genes
  
  # Fig 2A
  p = ggscatter(egene.count, x="sample.size", y="lrt.egenes", color="tissue")
  CairoPDF("projects/rare_variant/materials/egenes.vs.sample.size.pdf")
  ggpar(p, legend="none", ylab="Number of eGenes detected by LRT-q",
        xlab="Number of samples", font.tickslab = 13, font.legend = 13, 
        font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  
  # Fig 2B
  CairoPDF("projects/rare_variant/materials/rv.egenes.each.tissue.pdf", width = 12)
  t = melt(egene.count %>% filter(lrt.egenes > 0) 
           %>% as.data.table(), 
           id.vars = c("tissue", "sample.size"), 
           measure.vars = c("vt.egenes", "lrt.egenes", "skat.egenes", 
                                "acat.egenes", "acat.o.egenes"))
  t$tissue.i = rep(c(1:length(unique(t$tissue))), 5)
  t$variable = c(rep("VT", length(unique(t$tissue))), 
                 rep("LRT-q", length(unique(t$tissue))), 
                 rep("SKAT-O", length(unique(t$tissue))), 
                 rep("ACAT-V", length(unique(t$tissue))), 
                 rep("ACAT-O", length(unique(t$tissue))))
  p = ggbarplot(t %>% filter(variable != "ACAT-V"), 
                x="tissue", y="value", fill="variable", position = position_dodge()) + rotate_x_text(45)
  ggpar(p, ylab="Number of RV eGenes", legend.title = "Method", legend = "bottom",
        xlab="Tissue", font.xtickslab = 8, font.legend = 16, font.ytickslab = 16,
        font.main = 16, font.x = 16, font.y = 16)
  dev.off()
}
