# decision boundary
{
  library(data.table)
  library(ggpubr)
  library(Cairo)
  library(dplyr)
  file.dir = "./data/simulation/power/decision_boundary"
  
  c05 = rbindlist(lapply(dir(file.dir, pattern = "c05", full.names = T), fread))
  c1 = rbindlist(lapply(dir(file.dir, pattern = "c1", full.names = T), fread))
  
  c05 = na.omit(as.data.frame(c05))
  c05[, 7:14] = if_else(c05[, 7:14] < 0.05, "Significant", "Not significant")
  colnames(c05)[11:14] = c("SKATO", "ACATV", "ACATO", "LRTq")
  c05$truth = if_else(c05$causal1 == 0 & c05$causal2 == 0, 
                      "CE = 0", "CE > 0")
  
  # use transparency to avoid overplotting
  ggplot(c05, aes(x=reg.t.stat1, y=reg.t.stat2, color=LRTq, shape=truth)) +
    geom_point(alpha=0.3)
  
  # use sampling to avoid overplotting

  # Supplemental Fig S4A
  {
    folder = "../materials/"
    pic.name = paste(folder, "decision_boundary.LRTq.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", 
                   color = "LRTq", shape = "truth", ellipse.alpha = 0.3)
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = quote("LRT-q ("~c[i]~"= 0.5)"))
    dev.off()
  }
  # Supplemental Fig S4C
  {
    folder = "../materials/"
    pic.name = paste(folder, "decision_boundary.VT.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", 
                   color = "VT", shape = "truth")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = quote("VT ("~c[i]~"= 0.5)"))
    dev.off()
  }
  # Supplemental Fig S4B
  {
    folder = "../materials/"
    pic.name = paste(folder, "decision_boundary.SKATO.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", 
                   color = "SKATO", shape = "truth")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("SKAT-O ("~c[i]~"= 0.5)"))
    dev.off()
  }
  # Supplemental Fig S4D
  {
    folder = "../materials/"
    pic.name = paste(folder, "decision_boundary.CMC.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", 
                   color = "CMC", shape = "truth")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("CMC ("~c[i]~"= 0.5)"))
    dev.off()
  }
}

