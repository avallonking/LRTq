# number of rare variants within 1mb and 20kb from TSS
{
  library(data.table)
  library(Cairo)
  library(ggpubr)
  one_mb = rbindlist(lapply(dir("./data/rare_var_per_gene/1mb/", full.names = T), fread))
  twenty_kb = rbindlist(lapply(dir("./data/rare_var_per_gene/20kb/", full.names = T), fread))
  one_mb$dis = "1 Mb"
  twenty_kb$dis = "20 Kb"
  CairoPDF("../materials/rare_var_per_gene.pdf")
  p = ggboxplot(rbind(one_mb, twenty_kb), x = "dis", y = "V2", fill = "dis",
            xlab = "Distance from TSS", ylab = "Number of rare variants per gene")
  ggpar(p, legend.title = "Distance from TSS")
  dev.off()
}
