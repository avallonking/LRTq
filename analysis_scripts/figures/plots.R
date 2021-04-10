# power analysis - including ACAT
{
  require(data.table)
  require(dplyr)
  require(tidyr)
  require(ggpubr)
  require(Cairo)
  require(ACAT)
  setwd("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/")
  
  # type I error
  {
    typeIerror.check <- function(p.value) {
      fdr.table <- data.frame("alpha.level" = c(0.05, 0.01, 1e-3, 1e-4, 1e-5, 2.5e-6, 5e-6, 1e-6),
                              "fdr" = 0
      )
      for(i in 1:8) {
        fdr.table$fdr[i] <- sum(p.value < fdr.table$alpha.level[i]) / length(p.value)
      }
      return(fdr.table)
    }
    res = rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/type_I_error/acat/", full.names = T), fread))
    res$`ACAT-O` = apply(res, 1, function(x) ACAT(c(x[4], x[5], x[6])))
    fdr = apply(res, 2, typeIerror.check)
    fdr.table = cbind(fdr[[1]]$alpha.level, fdr[[1]]$fdr, fdr[[2]]$fdr, fdr[[3]]$fdr, 
                      fdr[[4]]$fdr, fdr[[5]]$fdr, fdr[[6]]$fdr, fdr[[7]]$fdr, 
                      fdr[[8]]$fdr)
    colnames(fdr.table) = c("alpha", names(fdr))
    fdr.table
  }
  
  # power
  {
    power = data.frame()
    for (a in c(0.3, 0.6, 0.9, 1.2, 1.5)) {
      for (r in c(0.03, 0.05, 0.07, "0.10")) {
        for (p in c("0.50")) {
          causal = paste("causal", r, sep = "")
          eff = paste("a", a, sep = "")
          p.rate = paste("protective_rate", p, sep = "")
          file.pattern = paste(causal, eff, "repeat1000", p.rate, "\\d+\\.csv", sep = ".")
          raw.data = rbindlist(lapply(dir("power/acat/", 
                                          pattern = file.pattern, full.names = T), 
                                      fread))
          raw.data$`ACAT-O` = apply(raw.data, 1, function(x) ACAT(c(x[4], x[5], x[6])))
          dfrow = c(colSums(raw.data < 0.05) / nrow(raw.data), a, r, p)
          names(dfrow) = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O", "ACAT-V", 
                           "LRT-q", "ACAT-O", "a", "causal", "protective.rate")
          power = bind_rows(power, dfrow)
        }
      }
    }
    power.transform = data.frame()
    for (method in colnames(power)[1:(ncol(power)-3)]) {
      sub.df = power[, c(method, "a", "causal", "protective.rate")]
      sub.df$method = method
      colnames(sub.df)[1] = "power"
      power.transform = rbind(power.transform, sub.df)
    }
    power.transform$power = as.numeric(power.transform$power)
    power.transform$method = factor(power.transform$method, 
                                    levels = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O",
                                               "ACAT-V", "LRT-q", "ACAT-O"))
    
    # protective rate = 50%
    CairoPDF("../materials/acat.effect.size.4.95.protective50.pdf")
    power.transform %>% filter(a == "1.5" & protective.rate == "0.50") %>% 
      ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
             ylab = "Power", legend.title = "Method", shape = "method", 
             title = expression(paste(alpha, "\ = 0.05, effect size"<=4.95))) -> p
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13)
    dev.off()
    
    CairoPDF("../materials/acat.effect.size.0.99.protective50.pdf")
    power.transform %>% filter(a == "0.3" & protective.rate == "0.50") %>% 
      ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
             ylab = "Power", legend.title = "Method", shape = "method", 
             title = expression(paste(alpha, "\ = 0.05, effect size"<=0.99))) -> p
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13)
    dev.off()
    
    CairoPDF("../materials/acat.effect.size.2.97.protective50.pdf")
    power.transform %>% filter(a == "0.9" & protective.rate == "0.50") %>% 
      ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
             ylab = "Power", legend.title = "Method", shape = "method", 
             title = expression(paste(alpha, "\ = 0.05, effect size"<=2.97))) -> p
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13)
    dev.off()
    
    CairoPDF("../materials/acat.causal0.10.protective50.pdf")
    power.transform %>% filter(causal == "0.10" & protective.rate == "0.50") %>%
      mutate(effect = as.numeric(a) * 3.3) %>%
      ggline(x = "effect", y = "power", color = "method",  shape = "method", 
             xlab = "Effect size upper bound", ylab = "Power", legend.title = "Method",
             title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%"))) -> p
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13)
    dev.off()
  }
}

# QQ-plot
{
  ggd.qqplot = function(pvector, main=NULL, ...) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
         xlab=expression(Expected~~-log[10](italic(p))),
         ylab=expression(Observed~~-log[10](italic(p))),
         xlim=c(0,max(e)), ylim=c(0,max(o)))
    lines(e,e,col="red")
  }
  
  lambda <- function(p.values) {
    return(median(qchisq(1 - na.omit(p.values), 1)) / qchisq(0.5, 1))
  }
  
  ggd.qqplot(lrt$Unweighted, main = paste("LRT (unweighted) lambda=", lambda(lrt$Unweighted), sep=""))
  ggd.qqplot(lrt$LINSIGHT, main = paste("LRT (LINSIGHT) lambda=", lambda(lrt$LINSIGHT), sep=""))
  ggd.qqplot(lrt$MAF, main = paste("LRT (MAF) lambda=", lambda(lrt$MAF), sep=""))
  ggd.qqplot(lrt$CADD, main = paste("LRT (CADD) lambda=", lambda(lrt$CADD), sep=""))
  
  ggd.qqplot(skat$Unweighted, main = paste("SKAT-O (unweighted) lambda=", lambda(skat$Unweighted), sep=""))
  ggd.qqplot(skat$LINSIGHT, main = paste("SKAT-O (LINSIGHT) lambda=", lambda(skat$LINSIGHT), sep=""))
  ggd.qqplot(skat$MAF, main = paste("SKAT-O (MAF) lambda=", lambda(skat$MAF), sep=""))
  ggd.qqplot(skat$CADD, main = paste("SKAT-O (CADD) lambda=", lambda(skat$CADD), sep=""))
}

# how number of permutations affect FDR
{
  require(fdrtool)
  setwd("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/Whole_Blood/modified_fixed/more_perm/")
  lrt.p <- fread("lrt.p", data.table = F)
  # only 100000 permutations
  lrt.p[lrt.p <= 1/100001] = 1/100001
  print(sum(fdrtool(lrt.p$Unweighted, statistic = "pvalue", plot = F)$qval < 0.05))
  print(sum(fdrtool(lrt.p$LINSIGHT, statistic = "pvalue", plot = F)$qval < 0.05))
  print(sum(fdrtool(lrt.p$MAF, statistic = "pvalue", plot = F)$qval < 0.05))
  print(sum(fdrtool(lrt.p$CADD, statistic = "pvalue", plot = F)$qval < 0.05))
}

{
  stat = 43
  m = 2000
  p = 2000
  pf(stat * (m-p+1) / (p*m), df1 = p, df2 = m)
}

# power analysis
{
  require(data.table)
  require(ggpubr)
  c05a3 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/mean_var_fixed/", pattern="causal0.05.a0.3", full.names = T), fread))
  c05a6 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/mean_var_fixed/", pattern="causal0.05.a0.6", full.names = T), fread))
  c10a3 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/mean_var_fixed/", pattern="causal0.10.a0.3", full.names = T), fread))
  c10a6 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/mean_var_fixed/", pattern="causal0.10.a0.6", full.names = T), fread))
  df <- data.frame(colSums(c05a3 < 0.05) / dim(c05a3)[1], 
                   colSums(c05a6 < 0.05) / dim(c05a6)[1], 
                   colSums(c10a3 < 0.05) / dim(c10a3)[1],
                   colSums(c10a6 < 0.05) / dim(c10a6)[1])
  colnames(df) <- c("causal ratio=0.05\na=0.3", "causal ratio=0.05\na=0.6", "causal ratio=0.10\na=0.3", "causal ratio=0.10\na=0.6")
  df$Method[1:3] <- c("SKAT-O", "LRT-orig", "LRT-modified")
  p <- ggbarplot(melt(df[1:3, ]), x="variable", y="value", fill = "Method", xlab = "Setting", ylab="Power", position = position_dodge(0.8))
  p
}

# power analysis - different protective rate
{
  library(data.table)
  library(dplyr)
  library(ggpubr)
  library(Cairo)
  data.dir <- "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest"
  power = data.frame()
  for (a in c(0.3, 0.6, 0.9, 1.2, 1.5)) {
    for (r in c(0.03, 0.05, 0.07, "0.10")) {
      for (p in c(0, 0.25, 0.50, 0.75)) {
        causal = paste("causal", r, sep = "")
        eff = paste("a", a, sep = "")
        if (p == 0.50) {
          file.pattern = paste(causal, eff, "repeat1000", "\\d+\\.csv", sep = ".")
        } else {
          p.rate = paste("protective_rate", p, sep = "")
          file.pattern = paste(causal, eff, "repeat1000", p.rate, "\\d+\\.csv", sep = ".")
        }
        raw.data = rbindlist(lapply(dir(data.dir, pattern = file.pattern, full.names = T), fread))
        dfrow = c(colSums(raw.data < 0.05) / nrow(raw.data), a, r, p)
        names(dfrow) = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O", "LRT-q", "a",
                         "causal", "protective.rate")
        power = bind_rows(power, dfrow)
      }
    }
  }
  power.transform = data.frame()
  for (method in colnames(power)[1:6]) {
    sub.df = power[, c(method, "a", "causal", "protective.rate")]
    sub.df$method = method
    colnames(sub.df)[1] = "power"
    power.transform = rbind(power.transform, sub.df)
  }
  power.transform$power = as.numeric(power.transform$power)
  power.transform$method = factor(power.transform$method, 
                                  levels = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O",
                                             "LRT-q"))
  
  # protective rate = 25%
  CairoPDF("projects/rare_variant/materials/effect.size.4.95.protective25.pdf")
  power.transform %>% filter(a == 1.5 & protective.rate == 0.25) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=4.95,
                                    ", protective rate = 25%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/effect.size.0.99.protective25.pdf")
  power.transform %>% filter(a == 0.3 & protective.rate == 0.25) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=0.99,
                                    ", protective rate = 25%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/effect.size.2.97.protective25.pdf")
  power.transform %>% filter(a == 0.9 & protective.rate == 0.25) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=2.97,
                                    ", protective rate = 25%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/causal0.10.protective25.pdf")
  power.transform %>% filter(causal == "0.10" & protective.rate == 0.25) %>%
    mutate(effect = as.numeric(a) * 3.3) %>%
    ggline(x = "effect", y = "power", color = "method", 
           xlab = "Effect size upper bound", ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%",
                                    ", protective rate = 25%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  
  # protective rate = 75%
  CairoPDF("projects/rare_variant/materials/effect.size.4.95.protective75.pdf")
  power.transform %>% filter(a == 1.5 & protective.rate == 0.75) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=4.95,
                                    ", protective rate = 75%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/effect.size.0.99.protective75.pdf")
  power.transform %>% filter(a == 0.3 & protective.rate == 0.75) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=0.99,
                                    ", protective rate = 75%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/effect.size.2.97.protective75.pdf")
  power.transform %>% filter(a == 0.9 & protective.rate == 0.75) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=2.97,
                                    ", protective rate = 75%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/causal0.10.protective75.pdf")
  power.transform %>% filter(causal == "0.10" & protective.rate == 0.75) %>%
    mutate(effect = as.numeric(a) * 3.3) %>%
    ggline(x = "effect", y = "power", color = "method", 
           xlab = "Effect size upper bound", ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%",
                                    ", protective rate = 75%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  
  # protective rate = 0%
  CairoPDF("projects/rare_variant/materials/effect.size.4.95.protective0.pdf")
  power.transform %>% filter(a == 1.5 & protective.rate == 0) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=4.95,
                                    ", protective rate = 0"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/effect.size.0.99.protective0.pdf")
  power.transform %>% filter(a == 0.3 & protective.rate == 0) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=0.99,
                                    ", protective rate = 0"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/effect.size.2.97.protective0.pdf")
  power.transform %>% filter(a == 0.9 & protective.rate == 0) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, effect size"<=2.97,
                                    ", protective rate = 0"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/causal0.10.protective0.pdf")
  power.transform %>% filter(causal == "0.10" & protective.rate == 0) %>%
    mutate(effect = as.numeric(a) * 3.3) %>%
    ggline(x = "effect", y = "power", color = "method", 
           xlab = "Effect size upper bound", ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%",
                                    ", protective rate = 0"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
}

# power analysis - fixed effect + different protective rate
{
  library(data.table)
  library(dplyr)
  library(ggpubr)
  library(Cairo)
  data.dir <- "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/fixed_effect/"
  power = data.frame()
  
  for (a in c(0.99, 1.98, 2.97, 3.96, 4.95)) {
    for (r in c(0.03, 0.05, 0.07, "0.10")) {
      for (p in c(0, 0.25, "0.50", 0.75)) {
        causal = paste("causal", r, sep = "")
        eff = paste("fixed_effect", a, sep = "")
        p.rate = paste("protective_rate", p, sep = "")
        file.pattern = paste(causal, eff, "repeat1000", p.rate, "\\d+\\.csv", sep = ".")
        
        raw.data = rbindlist(lapply(dir(data.dir, pattern = file.pattern, full.names = T), fread))
        dfrow = c(colSums(raw.data < 0.05) / nrow(raw.data), a, r, p)
        names(dfrow) = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O", "LRT-q", 
                         "fixed.effect", "causal", "protective.rate")
        power = bind_rows(power, dfrow)
      }
    }
  }
  power.transform = data.frame()
  for (method in colnames(power)[1:6]) {
    sub.df = power[, c(method, "fixed.effect", "causal", "protective.rate")]
    sub.df$method = method
    colnames(sub.df)[1] = "power"
    power.transform = rbind(power.transform, sub.df)
  }
  power.transform$power = as.numeric(power.transform$power)
  power.transform$method = factor(power.transform$method, 
                                  levels = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O",
                                             "LRT-q"))
  
  # protective rate = 25%
  CairoPDF("projects/rare_variant/materials/fixed_effect.4.95.protective25.pdf")
  power.transform %>% filter(fixed.effect == 4.95 & protective.rate == 0.25) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 4.95",
                                    ", protective rate = 25%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/fixed_effect.0.99.protective25.pdf")
  power.transform %>% filter(fixed.effect == 0.99 & protective.rate == 0.25) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 0.99",
                                    ", protective rate = 25%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/fixed_effect.2.97.protective25.pdf")
  power.transform %>% filter(fixed.effect == 2.97 & protective.rate == 0.25) %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 2.97",
                                    ", protective rate = 25%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/causal0.10.fixed.effect.protective25.pdf")
  power.transform %>% filter(causal == "0.10" & protective.rate == 0.25) %>%
    ggline(x = "fixed.effect", y = "power", color = "method", 
           xlab = "Effect size", ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%",
                                    ", protective rate = 25%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  
  # protective rate = 50%
  CairoPDF("projects/rare_variant/materials/fixed_effect.4.95.protective50.pdf")
  power.transform %>% filter(fixed.effect == 4.95 & protective.rate == "0.50") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 4.95",
                                    ", protective rate = 50%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/fixed_effect.0.99.protective50.pdf")
  power.transform %>% filter(fixed.effect == 0.99 & protective.rate == "0.50") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 0.99",
                                    ", protective rate = 50%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/fixed_effect.2.97.protective50.pdf")
  power.transform %>% filter(fixed.effect == 2.97 & protective.rate == "0.50") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 2.97",
                                    ", protective rate = 50%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/causal0.10.fixed.effect.protective50.pdf")
  power.transform %>% filter(causal == "0.10" & protective.rate == "0.50") %>%
    ggline(x = "fixed.effect", y = "power", color = "method", 
           xlab = "Effect size", ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%",
                                    ", protective rate = 50%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  
  # protective rate = 75%
  CairoPDF("projects/rare_variant/materials/fixed_effect.4.95.protective75.pdf")
  power.transform %>% filter(fixed.effect == 4.95 & protective.rate == "0.75") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 4.95",
                                    ", protective rate = 75%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/fixed_effect.0.99.protective75.pdf")
  power.transform %>% filter(fixed.effect == 0.99 & protective.rate == "0.75") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 0.99",
                                    ", protective rate = 75%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/fixed_effect.2.97.protective75.pdf")
  power.transform %>% filter(fixed.effect == 2.97 & protective.rate == "0.75") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 2.97",
                                    ", protective rate = 75%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/causal0.10.fixed.effect.protective75.pdf")
  power.transform %>% filter(causal == "0.10" & protective.rate == "0.75") %>%
    ggline(x = "fixed.effect", y = "power", color = "method", 
           xlab = "Effect size", ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%",
                                    ", protective rate = 75%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  
  # protective rate = 0%
  CairoPDF("projects/rare_variant/materials/fixed_effect.4.95.protective0.pdf")
  power.transform %>% filter(fixed.effect == 4.95 & protective.rate == "0") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 4.95",
                                    ", protective rate = 0"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/fixed_effect.0.99.protective0.pdf")
  power.transform %>% filter(fixed.effect == 0.99 & protective.rate == "0") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 0.99",
                                    ", protective rate = 0"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/fixed_effect.2.97.protective0.pdf")
  power.transform %>% filter(fixed.effect == 2.97 & protective.rate == "0") %>% 
    ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
           ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, fixed effect size = 2.97",
                                    ", protective rate = 0"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
  CairoPDF("projects/rare_variant/materials/causal0.10.fixed.effect.protective0.pdf")
  power.transform %>% filter(causal == "0.10" & protective.rate == "0") %>%
    ggline(x = "fixed.effect", y = "power", color = "method", 
           xlab = "Effect size", ylab = "Power", legend.title = "Method",
           title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%",
                                    ", protective rate = 0"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13)
  dev.off()
}

# power analysis - different sample size
{
  library(data.table)
  library(dplyr)
  library(ggpubr)
  library(Cairo)
  data.dir <- "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/different_sample_size/"
  power = data.frame()
  for (a in c(0.3, 0.6, 0.9, 1.2, 1.5)) {
    for (r in c(0.03, 0.05, 0.07, "0.10")) {
      for (s in seq(100, 1000, 50)) {
        causal = paste("causal", r, sep = "")
        eff = paste("a", a, sep = "")
        sample.size = paste("sample_size", s, sep = "")
        file.pattern = paste(causal, eff, "repeat1000", sample.size, "\\d+\\.csv", sep = ".")
        
        raw.data = rbindlist(lapply(dir(data.dir, pattern = file.pattern, full.names = T), fread))
        dfrow = c(colSums(raw.data < 0.05, na.rm = T) / nrow(raw.data), a, r, s)
        names(dfrow) = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O", "LRT-q", "a",
                         "causal", "sample.size")
        power = bind_rows(power, dfrow)
      }
    }
  }
  power.transform = data.frame()
  for (method in colnames(power)[1:6]) {
    sub.df = power[, c(method, "a", "causal", "sample.size")]
    sub.df$method = method
    colnames(sub.df)[1] = "power"
    power.transform = rbind(power.transform, sub.df)
  }
  power.transform$power = as.numeric(power.transform$power)
  power.transform$method = factor(power.transform$method, 
                                  levels = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O",
                                             "LRT-q"))
  power.transform$sample.size = as.numeric(power.transform$sample.size)
  
  CairoPDF("projects/rare_variant/materials/different_sample_size.causal10.a1.5.pdf")
  power.transform %>% filter(a == "1.5" & causal == "0.10") %>% 
    ggline(x = "sample.size", y = "power", color = "method", xlab = "Sample size",
           ylab = "Power", legend.title = "Method", numeric.x.axis = T,
           title = expression(paste(alpha, "\ = 0.05, effect size"<=4.95,
                                    ", causal ratio = 10%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13, xticks.by = 100)
  dev.off()
  CairoPDF("projects/rare_variant/materials/different_sample_size.causal5.a1.5.pdf")
  power.transform %>% filter(a == "1.5" & causal == "0.05") %>% 
    ggline(x = "sample.size", y = "power", color = "method", xlab = "Sample size",
           ylab = "Power", legend.title = "Method", numeric.x.axis = T,
           title = expression(paste(alpha, "\ = 0.05, effect size"<=4.95,
                                    ", causal ratio = 5%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13, xticks.by = 100)
  dev.off()
  CairoPDF("projects/rare_variant/materials/different_sample_size.causal10.a0.6.pdf")
  power.transform %>% filter(a == "0.6" & causal == "0.10") %>% 
    ggline(x = "sample.size", y = "power", color = "method", xlab = "Sample size",
           ylab = "Power", legend.title = "Method", numeric.x.axis = T,
           title = expression(paste(alpha, "\ = 0.05, effect size"<=1.98,
                                    ", causal ratio = 10%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13, xticks.by = 100)
  dev.off()
  CairoPDF("projects/rare_variant/materials/different_sample_size.causal5.a0.6.pdf")
  power.transform %>% filter(a == "0.6" & causal == "0.05") %>% 
    ggline(x = "sample.size", y = "power", color = "method", xlab = "Sample size",
           ylab = "Power", legend.title = "Method", numeric.x.axis = T,
           title = expression(paste(alpha, "\ = 0.05, effect size"<=1.98,
                                    ", causal ratio = 5%"))) -> p
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
        font.y = 13, xticks.by = 100)
  dev.off()
  
}

# latest type I error rate
{
  require(data.table)
  typeIerror.check <- function(p.value) {
    fdr.table <- data.frame("alpha.level" = c(0.05, 0.01, 1e-3, 1e-4, 1e-5, 2.5e-6, 5e-6, 1e-6),
                            "fdr" = 0
    )
    for(i in 1:8) {
      fdr.table$fdr[i] <- sum(p.value < fdr.table$alpha.level[i]) / length(p.value)
    }
    return(fdr.table)
  }
  res = rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/type_I_error/latest/", full.names = T), fread))
  fdr = apply(res, 2, typeIerror.check)
  fdr.table = cbind(fdr[[1]]$alpha.level, fdr[[1]]$fdr, fdr[[2]]$fdr, fdr[[3]]$fdr, 
                    fdr[[4]]$fdr, fdr[[5]]$fdr, fdr[[6]]$fdr)
  colnames(fdr.table) = c("alpha", names(fdr))
  fdr.table
}

{
  require(data.table)
  require(ggpubr)
  require(Cairo)
  c03a3 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.03.a0.3", full.names = T), fread))
  c03a6 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.03.a0.6", full.names = T), fread))
  c03a9 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.03.a0.9", full.names = T), fread))
  c03a12 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.03.a1.2", full.names = T), fread))
  c03a15 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.03.a1.5", full.names = T), fread))
  
  c05a3 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.05.a0.3", full.names = T), fread))
  c05a6 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.05.a0.6", full.names = T), fread))
  c05a9 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.05.a0.9", full.names = T), fread))
  c05a12 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.05.a1.2", full.names = T), fread))
  c05a15 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.05.a1.5", full.names = T), fread))
  
  c07a3 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.07.a0.3", full.names = T), fread))
  c07a6 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.07.a0.6", full.names = T), fread))
  c07a9 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.07.a0.9", full.names = T), fread))
  c07a12 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.07.a1.2", full.names = T), fread))
  c07a15 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.07.a1.5", full.names = T), fread))
  
  c10a3 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.10.a0.3", full.names = T), fread))
  c10a6 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.10.a0.6", full.names = T), fread))
  c10a9 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.10.a0.9", full.names = T), fread))
  c10a12 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.10.a1.2", full.names = T), fread))
  c10a15 <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/latest/", pattern="causal0.10.a1.5", full.names = T), fread))
  
  df <- data.frame(colSums(c03a3 < 0.05) / dim(c03a3)[1], 
                   colSums(c03a6 < 0.05) / dim(c03a6)[1], 
                   colSums(c03a9 < 0.05) / dim(c03a9)[1], 
                   colSums(c03a12 < 0.05) / dim(c03a12)[1], 
                   colSums(c03a15 < 0.05) / dim(c03a15)[1], 
                   
                   colSums(c05a3 < 0.05) / dim(c05a3)[1], 
                   colSums(c05a6 < 0.05) / dim(c05a6)[1], 
                   colSums(c05a9 < 0.05) / dim(c05a9)[1], 
                   colSums(c05a12 < 0.05) / dim(c05a12)[1], 
                   colSums(c05a15 < 0.05) / dim(c05a15)[1], 
                   
                   colSums(c07a3 < 0.05) / dim(c07a3)[1], 
                   colSums(c07a6 < 0.05) / dim(c07a6)[1], 
                   colSums(c07a9 < 0.05) / dim(c07a9)[1], 
                   colSums(c07a12 < 0.05) / dim(c07a12)[1], 
                   colSums(c07a15 < 0.05) / dim(c07a15)[1], 
                   
                   colSums(c10a3 < 0.05) / dim(c10a3)[1],
                   colSums(c10a6 < 0.05) / dim(c10a6)[1],
                   colSums(c10a9 < 0.05) / dim(c10a9)[1],
                   colSums(c10a12 < 0.05) / dim(c10a12)[1], 
                   colSums(c10a15 < 0.05) / dim(c10a15)[1])
  colnames(df) <- c("causal ratio=0.03\na=0.3", "causal ratio=0.03\na=0.6", 
                    "causal ratio=0.03\na=0.9", "causal ratio=0.03\na=1.2",
                    "causal ratio=0.03\na=1.5",
                    
                    "causal ratio=0.05\na=0.3", "causal ratio=0.05\na=0.6", 
                    "causal ratio=0.05\na=0.9", "causal ratio=0.05\na=1.2",
                    "causal ratio=0.05\na=1.5", 
                    
                    "causal ratio=0.07\na=0.3", "causal ratio=0.07\na=0.6", 
                    "causal ratio=0.07\na=0.9", "causal ratio=0.07\na=1.2",
                    "causal ratio=0.07\na=1.5", 
                    
                    "causal ratio=0.10\na=0.3", "causal ratio=0.10\na=0.6",
                    "causal ratio=0.10\na=0.9", "causal ratio=0.10\na=1.2",
                    "causal ratio=0.10\na=1.5")
  df$Method <- c("CMC", "WSS", "BURDEN", "VT", "SKAT-O", "LRT-q")
  df$Method <- factor(df$Method, levels = c("CMC", "WSS", "BURDEN", "VT", "SKAT-O", "LRT-q"))
  # summary
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/power.summary.pdf")
  p <- ggbarplot(melt(df), x="variable", y="value", fill = "Method", xlab = "Setting", ylab="Power", position = position_dodge(0.8))
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  
  # break down - causal = 10%
  temp.df = df[, c(16:20, 21)]
  colnames(temp.df) = c(0.99, 1.98, 2.97, 3.96, 4.95, "Method")
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/causal0.10.pdf")
  p1 <- ggline(melt(temp.df), x="variable", y="value", 
               color = "Method", xlab = "Effect size upper bound", ylab="Power",
               title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%")))
  ggpar(p1, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  # effect size <= 0.99
  temp.df = df[, c(1, 6, 11, 16, 21)]
  colnames(temp.df) = c("3%", "5%", "7%", "10%", "Method")
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/effect.size.0.99.pdf")
  p2 <- ggline(melt(temp.df), x="variable", y="value", 
               color = "Method", xlab = "Causal ratio", ylab="Power",
               title = expression(paste(alpha, "\ = 0.05, effect size"<=0.99)))
  ggpar(p2, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  # effect size <= 1.98
  temp.df = df[, c(2, 7, 12, 17, 21)]
  colnames(temp.df) = c("3%", "5%", "7%", "10%", "Method")
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/effect.size.1.98.pdf")
  p3 <- ggline(melt(temp.df), x="variable", y="value", 
               color = "Method", xlab = "Causal ratio", ylab="Power",
               title = expression(paste(alpha, "\ = 0.05, effect size"<=1.98)))
  ggpar(p3, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  # effect size <= 2.97
  temp.df = df[, c(3, 8, 13, 18, 21)]
  colnames(temp.df) = c("3%", "5%", "7%", "10%", "Method")
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/effect.size.2.97.pdf")
  p4 <- ggline(melt(temp.df), x="variable", y="value", 
               color = "Method", xlab = "Causal ratio", ylab="Power",
               title = expression(paste(alpha, "\ = 0.05, effect size"<=2.97)))
  ggpar(p4, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  # effect size <= 3.96
  temp.df = df[, c(4, 9, 14, 19, 21)]
  colnames(temp.df) = c("3%", "5%", "7%", "10%", "Method")
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/effect.size.3.96.pdf")
  p5 <- ggline(melt(temp.df), x="variable", y="value", 
               color = "Method", xlab = "Causal ratio", ylab="Power",
               title = expression(paste(alpha, "\ = 0.05, effect size"<=3.96)))
  ggpar(p5, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  # effect size <= 4.95
  temp.df = df[, c(5, 10, 15, 20, 21)]
  colnames(temp.df) = c("3%", "5%", "7%", "10%", "Method")
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/effect.size.4.95.pdf")
  p6 <- ggline(melt(temp.df), x="variable", y="value", 
               color = "Method", xlab = "Causal ratio", ylab="Power",
               title = expression(paste(alpha, "\ = 0.05, effect size"<=4.95)))
  ggpar(p6, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, font.y = 13)
  dev.off()
}

# outliers - deprecated
{
  require(data.table)
  require(ggpubr)
  require(biomaRt)
  require(dplyr)
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  exac <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/exac/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt")
  ensembl.id.symbol.table <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                                   filters = "hgnc_symbol", values = exac$gene,
                                   mart = ensembl)
  exac <- cbind(exac, ensembl.id.symbol.table[match(exac$gene, ensembl.id.symbol.table$hgnc_symbol), ])
  fwrite(exac, "~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/exac/ensembl_id_added.exac.constraints.metrics.txt")
  
  # lung.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.enrichment.fdr10//outliers.rare_var.enrichment.Lung")
  # adi.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.enrichment.fdr10//outliers.rare_var.enrichment.Adipose_Visceral_Omentum")
  # art.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.enrichment.fdr10//outliers.rare_var.enrichment.Artery_Tibial")
  # em.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.enrichment.fdr10//outliers.rare_var.enrichment.Esophagus_Muscularis")
  # ms.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.enrichment.fdr10//outliers.rare_var.enrichment.Muscle_Skeletal")
  # wb.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.enrichment.fdr10//outliers.rare_var.enrichment.Whole_Blood")
  # 
  # enrichment of rare variants
  {
    calculate.enrichment.helper <- function(outliers) {
      num.rare.outliers <- sum(outliers$rare.outliers)
      num.total.outliers <- sum(outliers$total.outliers)
      num.rare.controls <- sum(outliers$rare.controls)
      num.total.controls <- sum(outliers$total.controls)
      estimate <- (num.rare.outliers / num.total.outliers) / (num.rare.controls / num.total.controls)
      log.se <- sqrt(1/num.rare.outliers - 1/num.total.outliers + 1/num.rare.controls - 1/num.total.controls)
      max.ci <- estimate * exp(1.96 * log.se)
      min.ci <- estimate * exp(-1.96 * log.se)
      rare.outlier.prop <- num.rare.outliers / (num.total.controls + num.total.outliers)
      return(c(num.rare.outliers, num.total.outliers, num.rare.controls, num.total.controls,
               estimate, log.se, max.ci, min.ci, rare.outlier.prop))
    }
    
    calculate.enrichment <- function(outliers, tissue) {
      vt.idx <- which(outliers$vt == 1)
      skat.idx <- which(outliers$skat == 1)
      lrt.idx <- which(outliers$lrt == 1)
      vt.stat <- calculate.enrichment.helper(outliers[vt.idx, ])
      skat.stat <- calculate.enrichment.helper(outliers[skat.idx, ])
      lrt.stat <- calculate.enrichment.helper(outliers[lrt.idx, ])
      df <- as.data.frame(rbind(vt.stat, skat.stat, lrt.stat))
      colnames(df) <- c("num.rare.outliers", "num.total.outliers", 
                        "num.rare.controls", "num.total.controls",
                        "estimate", "log.se", "max.ci", "min.ci", "rare.outlier.prop")
      df$Methods <- c("VT", "SKAT-O", "LRT-q")
      df$Tissue <- tissue
      return(df)
    }
    # lung.enrich <- calculate.enrichment(lung.outliers, "Lung")
    # adi.enrich <- calculate.enrichment(adi.outliers, "Adipose - Visceral (Omentum)")
    # art.enrich <- calculate.enrichment(art.outliers, "Artery - Tibial")
    # em.enrich <- calculate.enrichment(em.outliers, "Esophagus - Muscularis")
    # ms.enrich <- calculate.enrichment(ms.outliers, "Muscle - Skeletal")
    # wb.enrich <- calculate.enrichment(wb.outliers, "Whole Blood")
    # result <- rbind(lung.enrich, adi.enrich, art.enrich, em.enrich, ms.enrich, wb.enrich)
    # 
    # CairoPDF("projects/rare_variant/materials/egenes.enrichment.for.rare.var.pdf")
    # result$Tissues <- c(rep("Lung", 3), rep("AVO", 3), rep("AT", 3), rep("EM", 3),
    #                     rep("MS", 3), rep("WB", 3))
    # p = ggplot(result, aes(x = Tissues, y = estimate)) +
    #   theme_bw() + ylab('Enrichment') +
    #   geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
    #   geom_pointrange(aes(x = Tissues, ymin = max.ci, ymax = min.ci, 
    #                       colour = Methods), 
    #                   position = position_dodge(width = 0.6)) +
    #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(), 
    #         axis.text.x = element_text(size = 10),
    #         axis.line = element_line(colour = "black"),
    #         text = element_text(size = 15),
    #         legend.position = "top") + 
    #   ggtitle("Enrichment for rare variants in RV eGenes")
    # print(p)
    # dev.off()
    # 

    # all tissues
    # sum.results <- rbind(colSums(result[which(result$Methods == "VT"), c(1:4)]),
    #                      colSums(result[which(result$Methods == "SKAT-O"), c(1:4)]),
    #                      colSums(result[which(result$Methods == "LRT-q"), c(1:4)]))
    # sum.results <- as.data.frame(sum.results)
    # sum.results$Methods <- c("VT", "SKAT-O", "LRT-q")
    # sum.results <- sum.results %>% mutate(Enrichment = (num.rare.outliers / num.total.outliers) / (num.rare.controls / num.total.controls))
    # sum.results <- sum.results %>% mutate(rare.outliers.prop = num.rare.outliers / (num.total.outliers + num.total.controls))
    # sum.results <- sum.results %>% mutate(log.se = sqrt(1/num.rare.outliers - 1/num.total.outliers + 1/num.rare.controls - 1/num.total.controls))
    # sum.results <- sum.results %>% mutate(max.ci = Enrichment * exp(1.96 * log.se))
    # sum.results <- sum.results %>% mutate(min.ci = Enrichment * exp(-1.96 * log.se))
    # CairoPDF("projects/rare_variant/materials/all.tissues.egenes.enrichment.for.rare.var.pdf")
    # p = ggplot(sum.results, aes(x = Methods, y = Enrichment)) +
    #   theme_bw() + ylab('Enrichment') +
    #   geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
    #   geom_pointrange(aes(x = Methods, ymin = max.ci, ymax = min.ci, 
    #                       colour = Methods), 
    #                   position = position_dodge(width = 0.6)) +
    #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(), legend.position = "none",
    #         # axis.text.x = element_text(angle = 90, hjust = 1),
    #         axis.line = element_line(colour = "black"))
    # 
    # print(p)
    # dev.off()
    
  }
  
  # outliers per gene (quantile normalized)
  {
    # proportion of outliers
    result$outlier.prop <- result$num.rare.outliers / (result$num.rare.controls + result$num.rare.outliers)
    ggbarplot(result, x = "Tissue", y = "outlier.prop", color = "Methods", 
              position = position_dodge(0.9))
    
    # outliers per egene
    plot.outliers.per.gene <- function(outliers.df, tissues) {
      df <- data.frame()
      for (i in 1:length(outliers.df)) {
        skat.df <- outliers.df[[i]] %>% filter(skat == 1) %>% dplyr::select(rare.outliers) %>%
          mutate(Methods = "SKAT-O", Tissue = tissues[i])
        vt.df <- outliers.df[[i]] %>% filter(vt == 1) %>% dplyr::select(rare.outliers) %>%
          mutate(Methods = "VT", Tissue = tissues[i])
        lrt.df <- outliers.df[[i]] %>% filter(lrt == 1) %>% dplyr::select(rare.outliers) %>%
          mutate(Methods = "LRT-q", Tissue = tissues[i])
        df <- rbind(df, skat.df, vt.df, lrt.df)  
      }
      colnames(df) <- c("rare.outliers", "Methods", "Tissue")
      df$rare.outliers <- as.numeric(df$rare.outliers)
      # p <- ggboxplot(df, x = "Tissue", y = "`Rare outliers per eGene`", color = "Methods")
      # p
      return(df)
    }
    rare.outlier.per.gene <- plot.outliers.per.gene(list(lung.outliers, 
                                                         adi.outliers, art.outliers, em.outliers,
                                                         ms.outliers, wb.outliers), c("Lung", 
                                                                                      "Adipose - Visceral (Omentum)", "Artery - Tibial",
                                                                                      "Esophagus - Muscularis", "Muscle - Skeletal", 
                                                                                      "Whole Blood"))
    
    # wilcox test
    for (t in unique(rare.outlier.per.gene$Tissue)) {
      print(t)
      rare.outlier.per.gene %>% filter(Tissue == t) %>% 
        filter(Methods == "LRT-q") %>% dplyr::pull(rare.outliers) -> lrt.outliers
      rare.outlier.per.gene %>% filter(Tissue == t) %>% 
        filter(Methods == "VT") %>% dplyr::pull(rare.outliers) -> vt.outliers
      rare.outlier.per.gene %>% filter(Tissue == t) %>% 
        filter(Methods == "SKAT-O") %>% dplyr::pull(rare.outliers) -> skat.outliers
      print("LRT vs SKAT")
      print(wilcox.test(lrt.outliers, skat.outliers)$p.value)
      print("LRT vs VT")
      print(wilcox.test(lrt.outliers, vt.outliers)$p.value)
    }
    
    CairoPDF("projects/rare_variant/materials/outliers.per.egene.pdf")
    rare.outlier.per.gene %>% group_by(Tissue, Methods) %>% 
      dplyr::summarise(mean = mean(rare.outliers), n = n(), 
                       se = sd(rare.outliers) / sqrt(n)) -> df
    df <- as.data.frame(df)
    df$Tissues <- c(rep("AVO", 3), rep("AT", 3), rep("EM", 3), rep("Lung", 3),
                    rep("MS", 3), rep("WB", 3))
    ###
    df %>% filter(Methods != "VT") -> df
    ###
    p <- ggplot(df, aes(x=Tissues, y=mean, colour=Methods)) + 
      geom_pointrange(aes(x = Tissues, ymin = mean + 1.96 * se, 
                          ymax = mean - 1.96 * se, 
                          colour = Methods), 
                      position = position_dodge(width = 0.6)) +
      ylab("Number of outliers") + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.text.x = element_text(size = 10),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 15),
            legend.position = "top") + 
      ggtitle("Outliers carrying rare variants in RV eGenes")
    print(p)
    
    dev.off()
  }
  
  # rare outliers z-scores
  {rare.z <- lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/rare.outlier.z", full.name=T), fread)
    names(rare.z) <- c("lrt.adi", "lrt.art", "lrt.em", "lrt.lung", "lrt.ms", "lrt.wb",
                       "skat.adi", "skat.art", "skat.em", "skat.lung", "skat.ms", "skat.wb",
                       "vt.adi", "vt.art", "vt.em", "vt.lung", "vt.ms", "vt.wb")
    rare.z.abs <- lapply(rare.z, function(x) {x$rare.z = abs(x$rare.z); x})
    df <- melt(rare.z.abs)
    name.split <- tstrsplit(df$L1, "\\.")
    df$Method <- name.split[[1]]
    df$Tissue <- name.split[[2]]
    ggboxplot(df, x = "Tissue", y = "value", color = "Method")}
  
  # using gene lists from other databases
  {
    # gene list preprocessing
    {
      # overlap with clinvar
      clinvar <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/clinvar/unique_gene_symbolclinvar.vcf", header=F, data.table = F)
      # clinvar.id.table <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
      #                           filters = "hgnc_symbol", values = clinvar$V1,
      #                           mart = ensembl)
      # overlap.clinvar <- as.data.frame(row.names = c("Lung", 
      #                                                "Adipose - Visceral (Omentum)", 
      #                                                "Artery - Tibial",
      #                                                "Esophagus - Muscularis", 
      #                                                "Muscle - Skeletal", 
      #                                                "Whole Blood"))
      # colnames(overlap.clinvar) <- c("VT", "SKAT-O", "LRT-q")
      # gene2phenotype
      g2p <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gene2phenotype", 
                                  full.names = T), fread))
      # gwas catalog (genes in REPORTED GENE(S) or MAPPED_GENE)
      gwas.catalog <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gwas_catalog/reported.genes")
      gwas.catalog %>% dplyr::select("REPORTED GENE(S)") %>% distinct() -> gwas.catalog
      # omim
      omim <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/omim/mim2gene.txt")
      colnames(omim) <- c("MIM_Number", "MIM_Entry_Type", "Entrez_Gene_ID", "HGNC", "Ensembl_Gene_ID")
      # orphanet
      require(XML)
      orphanet <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/orphanet/orphanet.genes",
                        header = F)
    }
    
    # overlap table
    {
      human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
      annotate.outliers <- function(mart, outliers) {
        genes <- outliers$gene
        entrez.hgnc <- getBM(attributes = c("entrezgene_id", "hgnc_id", "hgnc_symbol",
                                            "ensembl_gene_id"), 
                             filters = "ensembl_gene_id", mart = mart, values = genes)
        entrez.hgnc$hgnc_id <- tstrsplit(entrez.hgnc$hgnc_id, ":", keep = 2)[[1]]
        outliers$hgnc_symbol <- entrez.hgnc$hgnc_symbol[match(genes, entrez.hgnc$ensembl_gene_id)]
        outliers$hgnc_id <- entrez.hgnc$hgnc_id[match(genes, entrez.hgnc$ensembl_gene_id)]
        outliers$entrezgene_id <- entrez.hgnc$entrezgene_id[match(genes, entrez.hgnc$ensembl_gene_id)]
        return(outliers)
      }
      
      lung.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.rare_var.enrichment.Lung")
      adi.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.rare_var.enrichment.Adipose_Visceral_Omentum")
      art.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.rare_var.enrichment.Artery_Tibial")
      em.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.rare_var.enrichment.Esophagus_Muscularis")
      ms.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.rare_var.enrichment.Muscle_Skeletal")
      wb.outliers <- fread("projects/rare_variant/gtex/outliers/outliers.rare_var.enrichment.Whole_Blood")
      
      lung.outliers <- annotate.outliers(human, lung.outliers)
      adi.outliers <- annotate.outliers(human, adi.outliers)
      art.outliers <- annotate.outliers(human, art.outliers)
      em.outliers <- annotate.outliers(human, em.outliers)
      ms.outliers <- annotate.outliers(human, ms.outliers)
      wb.outliers <- annotate.outliers(human, wb.outliers)
      
      fwrite(lung.outliers, "projects/rare_variant/gtex/outliers/lung.outliers.annotated.csv")
      fwrite(adi.outliers, "projects/rare_variant/gtex/outliers/adi.outliers.annotated.csv")
      fwrite(art.outliers, "projects/rare_variant/gtex/outliers/art.outliers.annotated.csv")
      fwrite(em.outliers, "projects/rare_variant/gtex/outliers/em.outliers.annotated.csv")
      fwrite(ms.outliers, "projects/rare_variant/gtex/outliers/ms.outliers.annotated.csv")
      fwrite(wb.outliers, "projects/rare_variant/gtex/outliers/wb.outliers.annotated.csv")
      
      get.overlap.rate <- function(dis.genes, outliers, info.col) {
        outliers %>% select(info.col) -> outliers.genes
        outliers.genes <- as.data.frame(outliers.genes)[, 1]
        
        skat.overlap.num <- sum(outliers.genes[outliers$skat == 1] %in% dis.genes)
        skat.overlap.rate <- skat.overlap.num / sum(outliers$skat == 1)
        
        lrt.overlap.num <- sum(outliers.genes[outliers$lrt == 1] %in% dis.genes)
        lrt.overlap.rate <- lrt.overlap.num / sum(outliers$lrt == 1)
        
        vt.overlap.num <- sum(outliers.genes[outliers$vt == 1] %in% dis.genes)
        vt.overlap.rate <- vt.overlap.num / sum(outliers$vt == 1)
        
        row <- list("skat.overlap.num" = skat.overlap.num, 
                    "skat.overlap.prop" = skat.overlap.rate,
                    "lrt.overlap.num" = lrt.overlap.num,
                    "lrt.overlap.prop" = lrt.overlap.rate,
                    "vt.overlap.num" = vt.overlap.num,
                    "vt.overlap.rate" = vt.overlap.rate)
        return(row)
      }
      get.overlap.summary.table <- function(dis.gene.list, database.names, info.cols, 
                                            outliers, tissue) 
      {
        df <- data.frame("skat.overlap.num" = numeric(), 
                         "skat.overlap.prop" = numeric(),
                         "lrt.overlap.num" = numeric(),
                         "lrt.overlap.prop" = numeric(),
                         "vt.overlap.num" = numeric(),
                         "vt.overlap.rate" = numeric())
        for (i in 1:length(dis.gene.list)) {
          df <- rbind(df, get.overlap.rate(dis.gene.list[[i]], outliers, info.cols[i]))
        }
        df$database <- database.names
        df$tissue <- tissue
        return(df)
      }
      lung.overlap <- get.overlap.summary.table(list(clinvar$V1, g2p$`gene symbol`, 
                                                     omim$Ensembl_Gene_ID, orphanet$V1,
                                                     gwas.catalog$`REPORTED GENE(S)`),
                                                c("ClinVar", "G2P", "OMIM", "OrphaNet",
                                                  "GWAS"), 
                                                c("hgnc_symbol", "hgnc_symbol",
                                                  "gene", "hgnc_symbol",
                                                  "hgnc_symbol"),
                                                lung.outliers, "Lung"
      )
      adi.overlap <- get.overlap.summary.table(list(clinvar$V1, g2p$`gene symbol`, 
                                                    omim$Ensembl_Gene_ID, orphanet$V1,
                                                    gwas.catalog$`REPORTED GENE(S)`),
                                               c("ClinVar", "G2P", "OMIM", "OrphaNet",
                                                 "GWAS"), 
                                               c("hgnc_symbol", "hgnc_symbol",
                                                 "gene", "hgnc_symbol",
                                                 "hgnc_symbol"),
                                               adi.outliers, "Adipose - Visceral (Omentum)"
      )
      art.overlap <- get.overlap.summary.table(list(clinvar$V1, g2p$`gene symbol`, 
                                                    omim$Ensembl_Gene_ID, orphanet$V1,
                                                    gwas.catalog$`REPORTED GENE(S)`),
                                               c("ClinVar", "G2P", "OMIM", "OrphaNet",
                                                 "GWAS"), 
                                               c("hgnc_symbol", "hgnc_symbol",
                                                 "gene", "hgnc_symbol",
                                                 "hgnc_symbol"),
                                               art.outliers, "Artery - Tibial")
      em.overlap <- get.overlap.summary.table(list(clinvar$V1, g2p$`gene symbol`, 
                                                   omim$Ensembl_Gene_ID, orphanet$V1,
                                                   gwas.catalog$`REPORTED GENE(S)`),
                                              c("ClinVar", "G2P", "OMIM", "OrphaNet",
                                                "GWAS"), 
                                              c("hgnc_symbol", "hgnc_symbol",
                                                "gene", "hgnc_symbol",
                                                "hgnc_symbol"),
                                              em.outliers, "Esophagus - Muscularis")
      ms.overlap <- get.overlap.summary.table(list(clinvar$V1, g2p$`gene symbol`, 
                                                   omim$Ensembl_Gene_ID, orphanet$V1,
                                                   gwas.catalog$`REPORTED GENE(S)`),
                                              c("ClinVar", "G2P", "OMIM", "OrphaNet",
                                                "GWAS"), 
                                              c("hgnc_symbol", "hgnc_symbol",
                                                "gene", "hgnc_symbol",
                                                "hgnc_symbol"),
                                              ms.outliers, "Muscle - Skeletal")
      wb.overlap <- get.overlap.summary.table(list(clinvar$V1, g2p$`gene symbol`, 
                                                   omim$Ensembl_Gene_ID, orphanet$V1,
                                                   gwas.catalog$`REPORTED GENE(S)`),
                                              c("ClinVar", "G2P", "OMIM", "OrphaNet",
                                                "GWAS"), 
                                              c("hgnc_symbol", "hgnc_symbol",
                                                "gene", "hgnc_symbol",
                                                "hgnc_symbol"),
                                              wb.outliers, "Whole Blood")
      
      all.outliers <- rbind(lung.outliers, adi.outliers, art.outliers, em.outliers,
                            ms.outliers, wb.outliers)
      all.outliers %>% filter(skat == 1) %>% select(gene, hgnc_symbol, hgnc_id, entrezgene_id) %>%
        distinct() -> skat.all.egenes
      all.outliers %>% filter(vt == 1) %>% select(gene, hgnc_symbol, hgnc_id, entrezgene_id) %>%
        distinct() -> vt.all.egenes
      all.outliers %>% filter(lrt == 1) %>% select(gene, hgnc_symbol, hgnc_id, entrezgene_id) %>%
        distinct() -> lrt.all.egenes
      clinvar.overlap <- sapply(list(vt.all.egenes, skat.all.egenes, lrt.all.egenes), 
                                function(x) {num <- sum(x$hgnc_symbol %in% clinvar$V1); c(num, num/nrow(x))})
      g2p.overlap <- sapply(list(vt.all.egenes, skat.all.egenes, lrt.all.egenes), 
                            function(x) {num <- sum(x$hgnc_symbol %in% g2p$`gene symbol`); c(num, num/nrow(x))})
      omim.overlap <- sapply(list(vt.all.egenes, skat.all.egenes, lrt.all.egenes), 
                             function(x) {num <- sum((x$hgnc_symbol %in% omim$HGNC) | (x$gene %in% omim$Ensembl_Gene_ID) | (x$entrezgene_id %in% omim$Entrez_Gene_ID)); 
                             c(num, num/nrow(x))})
      orphanet.overlap <- sapply(list(vt.all.egenes, skat.all.egenes, lrt.all.egenes), 
                                 function(x) {num <- sum(x$hgnc_symbol %in% orphanet$V1); c(num, num/nrow(x))})
      gwas.overlap <- sapply(list(vt.all.egenes, skat.all.egenes, lrt.all.egenes), 
                             function(x) {num <- sum(x$hgnc_symbol %in% gwas.catalog$`REPORTED GENE(S)`); c(num, num/nrow(x))})
      
    }
    
    # gene list enrichment
    {
      gene.list.enrichment <- function(tissue, egenes = NULL, gene.list.attributes = NULL) 
      {
        if (is.null(gene.list.attributes)) {
          require(biomaRt)
          human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
          gene.list <- fread(paste("projects/rare_variant/gtex/gene.names/", tissue,
                                   ".gene.names", sep = ""), header = F)
          gene.list.attributes <- getBM(attributes = c("entrezgene_id", "hgnc_id", 
                                                       "hgnc_symbol", "ensembl_gene_id"), 
                                        filters = "ensembl_gene_id", mart = human, 
                                        values = gene.list$V1)
        }
        
        if (grepl(":", gene.list.attributes$hgnc_id[[1]])) {
          gene.list.attributes$hgnc_id <- tstrsplit(gene.list.attributes$hgnc_id, ":", keep = 2)[[1]]
        }
        
        if (is.null(egenes)) {
          egenes <- vector("list", 3)
          names(egenes) <- c("VT", "SKAT-O", "LRT-q")
          data <- fread(
            paste("projects/rare_variant/gtex/outliers/outliers.rare_var.enrichment",
                  tissue, sep = ".")
          )
          egenes[[1]] <- data %>% filter(vt == 1)
          egenes[[2]] <- data %>% filter(skat == 1)
          egenes[[3]] <- data %>% filter(lrt == 1)
          print(lapply(egenes, nrow))
        }
        
        gene.list.attributes$vt <- gene.list.attributes$ensembl_gene_id %in% egenes[[1]]$gene
        gene.list.attributes$skat <- gene.list.attributes$ensembl_gene_id %in% egenes[[2]]$gene
        gene.list.attributes$lrt <- gene.list.attributes$ensembl_gene_id %in% egenes[[3]]$gene
        gene.list.attributes$clinvar <- gene.list.attributes$hgnc_symbol %in% clinvar$V1
        gene.list.attributes$g2p <- gene.list.attributes$hgnc_symbol %in% g2p$`gene symbol`
        gene.list.attributes$omim <- (gene.list.attributes$hgnc_symbol %in% omim$HGNC) | 
          (gene.list.attributes$ensembl_gene_id %in% omim$Ensembl_Gene_ID) | 
          (gene.list.attributes$entrezgene_id %in% omim$Entrez_Gene_ID)
        gene.list.attributes$orphanet <- gene.list.attributes$hgnc_symbol %in% orphanet$V1
        gene.list.attributes$gwas <- gene.list.attributes$hgnc_symbol %in% gwas.catalog$`REPORTED GENE(S)`
        gene.list.attributes$all_databases <- 
          gene.list.attributes$clinvar | gene.list.attributes$g2p | gene.list.attributes$omim | gene.list.attributes$orphanet | gene.list.attributes$gwas
        
        res <- data.frame(odd.ratio = numeric(), min.conf = numeric(),
                          max.conf = numeric(), pval = numeric(), database = character(),
                          method = character(), tissue = character(), stringsAsFactors = F)
        options(stringsAsFactors = FALSE)
        for (i in c("vt", "skat", "lrt")) {
          for (t in c("clinvar", "g2p", "omim", "orphanet", "gwas", "all_databases")) {
            test <- fisher.test(table(gene.list.attributes[, c(i, t), with = F]))
            dfrow <- list(odd.ratio = as.numeric(test$estimate), 
                          max.conf = test$conf.int[2], min.conf = test$conf.int[1], 
                          pval = test$p.value, database = t, method = i, 
                          tissue = tissue)
            res <- rbind(res, dfrow)
          }
        }
        fwrite(gene.list.attributes, 
               paste("projects/rare_variant/gtex/outliers/annotated_genes/", 
                     tissue, ".attributes",
                     sep = ""))
        fwrite(res, 
               paste("projects/rare_variant/gtex/outliers/gene_set_enrichment/", 
                     tissue, ".enrichment",
                     sep = ""))
        return(res)
      }
      
      {
        lung.gene.names <- fread("projects/rare_variant/gtex/gene.names/lung.gene.names",
                                 header = F)
        lung.attributes <- getBM(attributes = c("entrezgene_id", "hgnc_id", 
                                                "hgnc_symbol", "ensembl_gene_id"), 
                                 filters = "ensembl_gene_id", mart = human, 
                                 values = lung.gene.names$V1)
        
        adi.gene.names <- fread("projects/rare_variant/gtex/gene.names/adipose.gene.names",
                                header = F)
        adi.attributes <- getBM(attributes = c("entrezgene_id", "hgnc_id", 
                                               "hgnc_symbol", "ensembl_gene_id"), 
                                filters = "ensembl_gene_id", mart = human, 
                                values = adi.gene.names$V1)
        
        art.gene.names <- fread("projects/rare_variant/gtex/gene.names/artery.gene.names",
                                header = F)
        art.attributes <- getBM(attributes = c("entrezgene_id", "hgnc_id", 
                                               "hgnc_symbol", "ensembl_gene_id"), 
                                filters = "ensembl_gene_id", mart = human, 
                                values = art.gene.names$V1)
        
        em.gene.names <- fread("projects/rare_variant/gtex/gene.names/em.gene.names",
                               header = F)
        em.attributes <- getBM(attributes = c("entrezgene_id", "hgnc_id", 
                                              "hgnc_symbol", "ensembl_gene_id"), 
                               filters = "ensembl_gene_id", mart = human, 
                               values = em.gene.names$V1)
        
        ms.gene.names <- fread("projects/rare_variant/gtex/gene.names/ms.gene.names",
                               header = F)
        ms.attributes <- getBM(attributes = c("entrezgene_id", "hgnc_id", 
                                              "hgnc_symbol", "ensembl_gene_id"), 
                               filters = "ensembl_gene_id", mart = human, 
                               values = ms.gene.names$V1)
        
        wb.gene.names <- fread("projects/rare_variant/gtex/gene.names/wb.gene.names",
                               header = F)
        wb.attributes <- getBM(attributes = c("entrezgene_id", "hgnc_id", 
                                              "hgnc_symbol", "ensembl_gene_id"), 
                               filters = "ensembl_gene_id", mart = human, 
                               values = wb.gene.names$V1)
      }
      
      {
        all.tissue.attributes <- rbind(lung.attributes, adi.attributes, art.attributes, 
                                       em.attributes, ms.attributes, wb.attributes)
        all.tissue.attributes <- 
          all.tissue.attributes[!duplicated(all.tissue.attributes$ensembl_gene_id), ]
        
        all.tissue.gene.list.enrichment <- 
          gene.list.enrichment("all.tissue", 
                               list(vt.all.egenes, skat.all.egenes, lrt.all.egenes),
                               all.tissue.attributes)
        {
          # lung.attributes <- fread("projects/rare_variant/gtex/outliers/annotated_genes/propably_wrong/lung.attributes")
          # adi.attributes <- fread("projects/rare_variant/gtex/outliers/annotated_genes/propably_wrong/adipose.attributes")
          # art.attributes <- fread("projects/rare_variant/gtex/outliers/annotated_genes/propably_wrong/artery.attributes")
          # em.attributes <- fread("projects/rare_variant/gtex/outliers/annotated_genes/propably_wrong/esophagus.attributes")
          # ms.attributes <- fread("projects/rare_variant/gtex/outliers/annotated_genes/propably_wrong/muscle.attributes")
          # wb.attributes <- fread("projects/rare_variant/gtex/outliers/annotated_genes/propably_wrong/blood.attributes")
          }
        # lung.gene.list.enrichment <- 
        #   gene.list.enrichment(tissue = "Lung", 
        #                        gene.list.attributes = lung.attributes)
        # adi.gene.list.enrichment <- 
        #   gene.list.enrichment(tissue = "Adipose_Visceral_Omentum", 
        #                        gene.list.attributes = adi.attributes)
        # art.gene.list.enrichment <- 
        #   gene.list.enrichment(tissue = "Artery_Tibial", 
        #                        gene.list.attributes = art.attributes)
        # em.gene.list.enrichment <- 
        #   gene.list.enrichment(tissue = "Esophagus_Muscularis", 
        #                        gene.list.attributes = em.attributes)
        # ms.gene.list.enrichment <- 
        #   gene.list.enrichment(tissue = "Muscle_Skeletal", 
        #                        gene.list.attributes = ms.attributes)
        # wb.gene.list.enrichment <- 
        #   gene.list.enrichment(tissue = "Whole_Blood", 
        #                        gene.list.attributes = wb.attributes)
      }
      
      # use this
      CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/gene.set.enrichment.pdf")
      all.tissue.gene.list.enrichment <- 
        fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gene_set_enrichment/should_be_correct/all.tissue.enrichment")
      all.tissue.gene.list.enrichment$database =
        rep(c("ClinVar", "G2P", "OMIM", "OrphaNet", "GWAS Catalog"), 3)
      all.tissue.gene.list.enrichment$method =
        c(rep("VT", 5), rep("SKAT-O", 5), rep("LRT-q", 5))
      p = ggplot(all.tissue.gene.list.enrichment %>% filter(method == "LRT-q"), aes(x = database, y = odd.ratio)) +
        theme_bw() + ylab('Odd ratio') + xlab("Database") +
        geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
        geom_pointrange(aes(x = database, ymin = min.conf, ymax = max.conf, colour = database), #colour = method), 
                        position = position_dodge(width = 0.6)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              # axis.text.x = element_text(angle = 90, hjust = 1),
              axis.line = element_line(colour = "black"), 
              text = element_text(size = 15),
              legend.position = "none") + coord_flip() +
        geom_text(aes(label=formatC(pval, digits = 2, format = "e")), 
                  vjust = -1.5, hjust = -0.05)
      # ggtitle("Enrichment for disease-associated genes detected by LRT-q") + coord_flip()
      print(p)
      dev.off() 
      
      CairoPDF("projects/rare_variant/materials/gene.set.enrichment.each.tissue.pdf")
      enrichment.table <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gene_set_enrichment/should_be_correct/", full.names = T), fread))
      enrichment.table %>% filter(database == "all_databases") %>% 
        filter(tissue != "all.tissue") -> enrichment.table
      enrichment.table %>% 
        mutate(
          method = gsub("vt", "VT", method), 
          method = gsub("skat", "SKAT-O", method), 
          method = gsub("lrt", "LRT-q", method)
        ) -> enrichment.table
      enrichment.table %>% 
        mutate(
          tissue = gsub("Adipose_Visceral_Omentum", "Adipose - Visceral (Omentum)", tissue), 
          tissue = gsub("Artery_Tibial", "Artery - Tibial", tissue), 
          tissue = gsub("Esophagus_Muscularis", "Esophagus - Muscularis", tissue),
          tissue = gsub("Lung", "Lung", tissue),
          tissue = gsub("Muscle_Skeletal", "Muscle - Skeletal", tissue),
          tissue = gsub("Whole_Blood", "Whole Blood", tissue)
        ) -> enrichment.table
      
      p = ggplot(enrichment.table, aes(x = tissue, y = odd.ratio)) +
        theme_bw() + ylab('Odd ratio') + xlab("Database") +
        geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
        geom_pointrange(aes(x = tissue, ymin = min.conf, ymax = max.conf, 
                            colour = method), 
                        position = position_dodge(width = 0.6)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              # axis.text.x = element_text(angle = 90, hjust = 1),
              axis.line = element_line(colour = "black"))
      print(p)
      dev.off() 
    }
    
    # novel genes other than those in databases
    {
      novel.egenes <- function(path) {
        res <- data.frame(vt = numeric(), skat = numeric(), lrt = numeric(),
                          vt.prop = numeric(), skat.prop = numeric(), 
                          lrt.prop = numeric(), tissue = character(), 
                          stringsAsFactors = F)
        options(stringsAsFactors = FALSE)
        data <- lapply(dir(path, full.names = T), fread)
        names(data) <- 
          tstrsplit(basename(dir(path)), 
                    "\\.", keep = 1)[[1]]
        for (i in 1:length(data)) {
          data[[i]] %>% 
            mutate(novel = !(clinvar | g2p | omim | orphanet | gwas)) -> 
            data[[i]]
          data[[i]] %>% filter(vt & novel) %>% count() %>% as.numeric() -> vt
          data[[i]] %>% filter(skat & novel) %>% count() %>% as.numeric() -> skat
          data[[i]] %>% filter(lrt & novel) %>% count() %>% as.numeric() -> lrt
          lrt.prop <- lrt / sum(data[[i]]$lrt)
          skat.prop <- skat / sum(data[[i]]$skat)
          vt.prop <- vt / sum(data[[i]]$vt)
          res <- rbind(res, list(vt = vt, skat = skat, lrt = lrt, vt.prop = vt.prop,
                                 skat.prop = skat.prop, lrt.prop = lrt.prop,
                                 tissue = names(data)[i]))
        }
        return(res)
      }
      tissues.novel.egenes <- 
        novel.egenes("projects/rare_variant/gtex/outliers/annotated_genes/should_be_correct/")
      fwrite(tissues.novel.egenes, 
             "projects/rare_variant/gtex/outliers/novel_genes/novel.genes.by.tissue.csv")
    }
  }
  
  # unnormalized data (avoid shrinkage of outlier expression values caused by quantile normalization)
  # use TPM that are corrected with covariates
  {
    require(data.table)
    require(ggpubr)
    require(tidyr)
    setwd("projects/rare_variant/gtex/outliers/")
    indiv.with.rare.var <- lapply(dir("indiv.with.rare.var.idx/", full.names = T), fread)
    vt.rare.z <- lapply(dir("tpm/", pattern = "vt", full.names = T), fread)
    skat.rare.z <- lapply(dir("tpm/", pattern = "skat", full.names = T), fread)
    lrt.rare.z <- lapply(dir("tpm/", pattern = "lrt", full.names = T), fread)
    names(vt.rare.z) <- c("AVO", "AT", "EM", "Lung", "MS", "WB")
    names(skat.rare.z) <- c("AVO", "AT", "EM", "Lung", "MS", "WB")
    names(lrt.rare.z) <- c("AVO", "AT", "EM", "Lung", "MS", "WB")
    names(indiv.with.rare.var) <- c("AVO", "AT", "EM", "Lung", "MS", "WB")
    
    compute.outliers.enrichment <- function(sample.tpm, sample.rare) {
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
        estimate <- (num.rare.outliers / num.total.outliers) / (num.rare.controls / num.total.controls)
        log.se <- sqrt(1/num.rare.outliers - 1/num.total.outliers + 1/num.rare.controls - 1/num.total.controls)
        max.ci <- estimate * exp(1.96 * log.se)
        min.ci <- estimate * exp(-1.96 * log.se)
        result = rbind(result, c(num.total.outliers, num.total.controls, num.rare.outliers, num.rare.controls, 
                                 num.common.outliers, num.common.controls, estimate, log.se, max.ci, min.ci, 
                                 threshold))
      }
      colnames(result) = c("total.outliers", "total.controls", "rare.outliers", "rare.controls", "common.outliers", 
                           "common.controls", "estimate", "log.se", "max.ci", "min.ci", "threshold")
      return(result)
    }
    
    standard.z <- function(x) {
      return((x - mean(x)) / sd(x))
    }
    outliers.enrichment = vector("list", length(indiv.with.rare.var))
    names(outliers.enrichment) <- names(indiv.with.rare.var)
    for (tissue in names(indiv.with.rare.var)) {
      outliers.enrichment[[tissue]] = compute.outliers.enrichment(lrt.rare.z[[tissue]], 
                                                                  indiv.with.rare.var[[tissue]])
    }
  }
  
  # should use this!
  # transformation: log2(TPM + 1) -> standardization -> corrected with covariates
  {
    require(data.table)
    require(ggpubr)
    require(tidyr)
    require(dplyr)
    require(Cairo)
    setwd("/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/")
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
      names(indiv.with.rare.var) <- tstrsplit(dir("indiv.with.rare.var.idx/"), "\\.")[[5]]
      names(lrt.rare.z) <- tstrsplit(dir("log2.standardized.corrected.tpm/", pattern = "lrt"), "\\.")[[6]]
      names(vt.rare.z) <- tstrsplit(dir("log2.standardized.corrected.tpm/", pattern = "vt"), "\\.")[[6]]
      names(skat.rare.z) <- tstrsplit(dir("log2.standardized.corrected.tpm/", pattern = "skat"), "\\.")[[6]]
      
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
    
    # lrt enrichment of rare variants per tissue - plot
    {
      df = rbindlist(lapply(outliers.enrichment, function(x) x[2, ]))
      df$tissue = names(outliers.enrichment)
      CairoPDF("../../materials/log_tpm.egenes.enrichment.for.rare.var.per.tissue.lrt.only.pdf", 
               width = 12, height = 8)
      p = ggplot(df, aes(x = tissue, y = relative.risk)) +
        theme_bw() + ylab('Enrichment of rare variants') + xlab("Tissue") +
        geom_abline(intercept = 1, slope = 0, linetype = "dashed") + #ylim(0.95, 1.15) +
        geom_pointrange(aes(x = tissue, ymin = min.ci, ymax = max.ci, colour = tissue)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), legend.position = "none",
              axis.text.x = element_text(size = 10), text = element_text(size = 15),
              axis.line = element_line(colour = "black")) + coord_flip() #+ 
        # geom_text(aes(label=formatC(pval, digits = 2, format = "e")), vjust = -1.5)
      print(p)
      dev.off()
    }
    
    # using gene lists from other databases
    {
      # gene list preprocessing
      {
        # overlap with clinvar
        clinvar <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/clinvar/unique_gene_symbolclinvar.vcf", header=F, data.table = F)
        g2p <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gene2phenotype", 
                                    full.names = T), fread))
        # gwas catalog (genes in REPORTED GENE(S) or MAPPED_GENE)
        gwas.catalog <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gwas_catalog/reported.genes")
        gwas.catalog %>% dplyr::select("REPORTED GENE(S)") %>% distinct() -> gwas.catalog
        # omim
        omim <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/omim/mim2gene.txt")
        colnames(omim) <- c("MIM_Number", "MIM_Entry_Type", "Entrez_Gene_ID", "HGNC", "Ensembl_Gene_ID")
        # orphanet
        require(XML)
        orphanet <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/orphanet/orphanet.genes",
                          header = F)
      }
      
      # annotate outliers
      {
        human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
        annotate.outliers <- function(mart, outliers) {
          genes <- outliers$gene
          entrez.hgnc <- getBM(attributes = c("entrezgene_id", "hgnc_id", "hgnc_symbol",
                                              "ensembl_gene_id"), 
                               filters = "ensembl_gene_id", mart = mart, values = genes)
          entrez.hgnc$hgnc_id <- tstrsplit(entrez.hgnc$hgnc_id, ":", keep = 2)[[1]]
          outliers$hgnc_symbol <- entrez.hgnc$hgnc_symbol[match(genes, entrez.hgnc$ensembl_gene_id)]
          outliers$hgnc_id <- entrez.hgnc$hgnc_id[match(genes, entrez.hgnc$ensembl_gene_id)]
          outliers$entrezgene_id <- entrez.hgnc$entrezgene_id[match(genes, entrez.hgnc$ensembl_gene_id)]
          return(outliers)
        }
        
        annotated.outliers <- lapply(outliers.per.gene$lrt, function(x) annotate.outliers(human, x))
        
        saveRDS(annotated.outliers, "~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/all_tissues.outliers.annotated.lrt_only.rds")
      }
      
      # gene list enrichment
      {
        gene.list.enrichment <- function(tissue, egenes, gene.list.attributes, early.output=F) 
        {
          gwas.catalog %>% filter(!(`REPORTED GENE(S)` == "" | is.na(`REPORTED GENE(S)`))) -> gwas.catalog
          clinvar %>% filter(!(V1 == "" | is.na(V1))) -> clinvar
          g2p %>% filter(!(`gene symbol` == "" | is.na(`gene symbol`))) -> g2p
          orphanet %>% filter(!(V1 == "" | is.na(V1))) -> orphanet
          omim$Entrez_Gene_ID[is.na(omim$Entrez_Gene_ID)] = 0
          omim$HGNC[omim$HGNC == ""] = "_"
          omim$Ensembl_Gene_ID[omim$Ensembl_Gene_ID == ""] = "_"
          
          if (grepl(":", gene.list.attributes$hgnc_id[[1]])) {
            gene.list.attributes$hgnc_id <- tstrsplit(gene.list.attributes$hgnc_id, ":", keep = 2)[[1]]
          }
          
          gene.list.attributes$lrt <- gene.list.attributes$gene %in% egenes
          gene.list.attributes$non.lrt <- 1 - gene.list.attributes$lrt
          gene.list.attributes$clinvar <- gene.list.attributes$hgnc_symbol %in% clinvar$V1
          gene.list.attributes$g2p <- gene.list.attributes$hgnc_symbol %in% g2p$`gene symbol`
          gene.list.attributes$omim <- (gene.list.attributes$hgnc_symbol %in% omim$HGNC) | 
            (gene.list.attributes$gene %in% omim$Ensembl_Gene_ID) | 
            (gene.list.attributes$entrezgene_id %in% omim$Entrez_Gene_ID)
          gene.list.attributes$orphanet <- gene.list.attributes$hgnc_symbol %in% orphanet$V1
          gene.list.attributes$gwas <- gene.list.attributes$hgnc_symbol %in% gwas.catalog$`REPORTED GENE(S)`
          gene.list.attributes$all_databases <- 
            gene.list.attributes$clinvar | gene.list.attributes$g2p | gene.list.attributes$omim | gene.list.attributes$orphanet | gene.list.attributes$gwas
          
          if (early.output) {
            fwrite(gene.list.attributes, 
                   "~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/all.tissues.gene.list.attributes")
          }
          
          res <- data.frame(odd.ratio = numeric(), min.conf = numeric(),
                            max.conf = numeric(), pval = numeric(), database = character(),
                            method = character(), tissue = character(), 
                            odd.ratio.non.lrt = numeric(),
                            max.conf.non.lrt = numeric(),
                            min.conf.non.lrt = numeric(),
                            pval.non.lrt = numeric(),
                            stringsAsFactors = F)
          options(stringsAsFactors = FALSE)
          for (t in c("clinvar", "g2p", "omim", "orphanet", "gwas", "all_databases")) {
            test <- fisher.test(table(gene.list.attributes[, c("lrt", t)]))
            test2 <- fisher.test(table(gene.list.attributes[, c("non.lrt", t)]))
            dfrow <- list(odd.ratio = as.numeric(test$estimate), 
                          max.conf = test$conf.int[2], min.conf = test$conf.int[1], 
                          pval = test$p.value, database = t, method = "lrt", 
                          tissue = tissue, 
                          odd.ratio.non.lrt = as.numeric(test2$estimate),
                          max.conf.non.lrt = test2$conf.int[2],
                          min.conf.non.lrt = test2$conf.int[1],
                          pval.non.lrt = test2$p.value)
            res <- rbind(res, dfrow)
          }
          return(res)
        }
        
        {
          library(biomaRt)
          human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
          annotate.outliers <- function(mart, outliers) {
            genes <- outliers$gene
            entrez.hgnc <- getBM(attributes = c("entrezgene_id", "hgnc_id", "hgnc_symbol",
                                                "ensembl_gene_id"), 
                                 filters = "ensembl_gene_id", mart = mart, values = genes)
            entrez.hgnc$hgnc_id <- tstrsplit(entrez.hgnc$hgnc_id, ":", keep = 2)[[1]]
            outliers$hgnc_symbol <- entrez.hgnc$hgnc_symbol[match(genes, entrez.hgnc$ensembl_gene_id)]
            outliers$hgnc_id <- entrez.hgnc$hgnc_id[match(genes, entrez.hgnc$ensembl_gene_id)]
            outliers$entrezgene_id <- entrez.hgnc$entrezgene_id[match(genes, entrez.hgnc$ensembl_gene_id)]
            return(outliers)
          }
          total.expressed.genes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/lrt.all_tissues.expresed_genes", 
                                        data.table=F, header=F)
          lrt.egenes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/lrt.all_tissues.egenes", 
                             data.table=F, header=F)
          colnames(total.expressed.genes) = "gene"
          colnames(lrt.egenes) = "gene"
          lrt.egenes = annotate.outliers(human, lrt.egenes)
          total.expressed.genes = annotate.outliers(human, total.expressed.genes)
          saveRDS(total.expressed.genes, "~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/total.annotated.genes.rds")
          total.expressed.genes = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/total.annotated.genes.rds")
          all.tissue.gene.list.enrichment = gene.list.enrichment("all.tissues", lrt.egenes$gene, total.expressed.genes)
        }
        
        # use this
        CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/gene.set.enrichment.pdf")
        all.tissue.gene.list.enrichment$database = c("ClinVar", "G2P", "OMIM", "OrphaNet", "GWAS Catalog", "All Databases")
        all.tissue.gene.list.enrichment$method = rep("LRT-q", 6)
        df = all.tissue.gene.list.enrichment[1:5, 1:5]
        df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[1, c(8:11, 1)]))
        df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[2, c(8:11, 1)]))
        df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[3, c(8:11, 1)]))
        df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[4, c(8:11, 1)]))
        df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[5, c(8:11, 1)]))
        df$database[6:10] = all.tissue.gene.list.enrichment$database[1:5]
        df$type = c(rep("RV eGenes", 5), rep("Non-RV eGenes", 5))
        p = ggplot(df %>% filter(type == "RV eGenes"), aes(x = database, y = odd.ratio)) +
          theme_bw() + ylab('Odd ratio') + xlab("Database") +
          geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
          geom_pointrange(aes(x = database, y=odd.ratio, ymin = min.conf, ymax = max.conf, color = database), #colour = method), 
                          position = position_dodge(width = 0.6)) +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                # axis.text.x = element_text(angle = 90, hjust = 1),
                axis.line = element_line(colour = "black"), 
                text = element_text(size = 15),
                legend.position = "none") + coord_flip() + labs(color=NULL) + ylim(0, 3.5) +
          geom_text(aes(label=formatC(pval, digits = 2, format = "e")),
                    vjust = -1.5, hjust = -0.05)
        # ggtitle("Enrichment for disease-associated genes detected by LRT-q") + coord_flip()
        print(p)
        dev.off() 
        
        # gene set enrichment per tissue - not useful, temporarily
        {
          CairoPDF("projects/rare_variant/materials/gene.set.enrichment.each.tissue.pdf")
          enrichment.table <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gene_set_enrichment/should_be_correct/", full.names = T), fread))
          enrichment.table %>% filter(database == "all_databases") %>% 
            filter(tissue != "all.tissue") -> enrichment.table
          enrichment.table %>% 
            mutate(
              method = gsub("vt", "VT", method), 
              method = gsub("skat", "SKAT-O", method), 
              method = gsub("lrt", "LRT-q", method)
            ) -> enrichment.table
          enrichment.table %>% 
            mutate(
              tissue = gsub("Adipose_Visceral_Omentum", "Adipose - Visceral (Omentum)", tissue), 
              tissue = gsub("Artery_Tibial", "Artery - Tibial", tissue), 
              tissue = gsub("Esophagus_Muscularis", "Esophagus - Muscularis", tissue),
              tissue = gsub("Lung", "Lung", tissue),
              tissue = gsub("Muscle_Skeletal", "Muscle - Skeletal", tissue),
              tissue = gsub("Whole_Blood", "Whole Blood", tissue)
            ) -> enrichment.table
          
          p = ggplot(enrichment.table, aes(x = tissue, y = odd.ratio)) +
            theme_bw() + ylab('Odd ratio') + xlab("Database") +
            geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
            geom_pointrange(aes(x = tissue, ymin = min.conf, ymax = max.conf, 
                                colour = method), 
                            position = position_dodge(width = 0.6)) +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  # axis.text.x = element_text(angle = 90, hjust = 1),
                  axis.line = element_line(colour = "black"))
          print(p)
          dev.off() 
        }
      }
      
    }
    
    # sum over all tissues
    {
      # summary.outliers.enrichment = outliers.enrichment$AVO[, 1:6] + outliers.enrichment$AT[, 1:6] +
        # outliers.enrichment$EM[, 1:6] + outliers.enrichment$Lung[, 1:6] + outliers.enrichment$MS[, 1:6] +
        # outliers.enrichment$WB[, 1:6]
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
      
      CairoPDF("../../materials/log_tpm.egenes.enrichment.for.rare.var.pdf")
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
    
    {
      # for each gene
      {
        get.outliers.per.gene <- function(sample.tpm, sample.rare, threshold = 2) {
          if (!any(class(sample.tpm) == "data.frame")) {
            return()
          }
          # not compute for each gene, but for the whole dataset
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
          result = data.frame(
            "gene" = selected.gene,
            "total.outliers" = num.total.outliers,
            "total.controls" = num.total.controls,
            "rare.outliers" = num.rare.outliers,
            "rare.controls" = num.rare.controls,
            "common.outliers" = num.common.outliers,
            "common.controls" = num.common.controls
          )
          return(result)
        }
        {
          outliers.per.gene <- vector("list", 3)
          names(outliers.per.gene) <- c("vt", "skat", "lrt")
          outliers.per.gene[["vt"]] <- vector("list", 48)
          outliers.per.gene[["lrt"]] <- vector("list", 48)
          outliers.per.gene[["skat"]] <- vector("list", 48)
          names(outliers.per.gene[["vt"]]) <- names(indiv.with.rare.var)
          names(outliers.per.gene[["lrt"]]) <- names(indiv.with.rare.var)
          names(outliers.per.gene[["skat"]]) <- names(indiv.with.rare.var)
          for (tissue in names(indiv.with.rare.var)) {
            outliers.per.gene[["vt"]][[tissue]] <- get.outliers.per.gene(vt.rare.z[[tissue]],
                                                                         indiv.with.rare.var[[tissue]])
            outliers.per.gene[["lrt"]][[tissue]] <- get.outliers.per.gene(lrt.rare.z[[tissue]], 
                                                                          indiv.with.rare.var[[tissue]])
            outliers.per.gene[["skat"]][[tissue]] <- get.outliers.per.gene(skat.rare.z[[tissue]],
                                                                           indiv.with.rare.var[[tissue]])
          }
        }
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
      
      # outlier proportion - not useful
      {
        result$outlier.prop <- result$mean.rare.outliers / (result$mean.rare.controls + result$mean.rare.outliers)
        ggbarplot(na.omit(result), x = "tissue", y = "outlier.prop", color = "method", 
                  position = position_dodge(0.9))
      }
      
      # wilcox test - not useful
      for (t in unique(names(indiv.with.rare.var))) {
        print(t)
        if (t %in% names(outliers.per.gene$lrt)) {
          outliers.per.gene$lrt[[t]] %>% dplyr::pull(rare.outliers) -> lrt.outliers
        } else {
          lrt.outliers = NULL
        }
        if (t %in% names(outliers.per.gene$vt)) {
          outliers.per.gene$vt[[t]] %>% dplyr::pull(rare.outliers) -> vt.outliers
        } else {
          vt.outliers = NULL
        }
        if (t %in% names(outliers.per.gene$skat)) {
          outliers.per.gene$skat[[t]] %>% dplyr::pull(rare.outliers) -> skat.outliers
        } else {
          skat.outliers = NULL
        }
        
        print("LRT vs SKAT")
        if (!is.null(lrt.outliers) & !is.null(skat.outliers)) {
          print(wilcox.test(lrt.outliers, skat.outliers)$p.value)
        } else {
          print(NULL)
        }
        print("LRT vs VT")
        if (!is.null(lrt.outliers) & !is.null(vt.outliers)) {
          print(wilcox.test(lrt.outliers, vt.outliers)$p.value)
        } else {
          print(NULL)
        }
      }
      
      # outliers per tissue
      {
        CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/outliers.per.tissue.pdf")
        p <- ggplot(na.omit(result), aes(x=tissue, y=mean.rare.outliers, colour=method)) + 
          geom_pointrange(aes(x = tissue, ymin = mean.rare.outliers + 1.96 * se.rare.outliers, 
                              ymax = mean.rare.outliers - 1.96 * se.rare.outliers, 
                              colour = method), 
                          position = position_dodge(width = 0.6)) +
          ylab("Average number of outliers") + theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.text.x = element_text(size = 10),
                axis.line = element_line(colour = "black"),
                text = element_text(size = 15),
                legend.position = "top") + xlab("Tissue") + labs(colour = "Methods") + coord_flip()
        # ggtitle("Outliers carrying rare variants in RV eGenes")
        print(p)
        dev.off()
        
        CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/outliers.per.tissue.lrt.only.pdf")
        res = na.omit(result)
        p <- ggplot(res %>% filter(method == "LRT-q") %>% arrange(mean.rare.outliers) %>% 
                      mutate(tissue=factor(tissue, levels=tissue)), 
                    aes(x=tissue, y=mean.rare.outliers, colour=tissue)) + 
          geom_pointrange(aes(x = tissue, ymin = mean.rare.outliers + 1.96 * se.rare.outliers, 
                              ymax = mean.rare.outliers - 1.96 * se.rare.outliers, 
                              colour = tissue), 
                          position = position_dodge(width = 0.6)) +
          ylab("Average number of outliers") + theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.text.x = element_text(size = 10),
                axis.line = element_line(colour = "black"),
                text = element_text(size = 15), 
                legend.position = "none") + xlab("Tissue") + coord_flip()
        # ggtitle("Outliers carrying rare variants in RV eGenes")
        print(p)
        dev.off()
        
        # by proportion of total samples
        result$mean.rare.outliers.prop = 
          result$mean.rare.outliers / (result$mean.total.outliers + result$mean.total.controls)
        result$se.rare.outliers.prop = 
          result$se.rare.outliers / (result$mean.total.outliers + result$mean.total.controls)
        CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/outliers.prop.by.total.samples.per.tissue.lrt.only.pdf",
                 width = 12, height = 8)
        p <- ggplot(na.omit(result) %>% filter(method == "LRT-q") %>% arrange(mean.rare.outliers.prop) %>% 
                      mutate(tissue=factor(tissue, levels=tissue)), 
                    aes(x=tissue, y=mean.rare.outliers.prop, colour=tissue)) + 
          geom_pointrange(aes(x = tissue, ymin = mean.rare.outliers.prop + 1.96 * se.rare.outliers.prop, 
                              ymax = mean.rare.outliers.prop - 1.96 * se.rare.outliers.prop, 
                              colour = tissue), 
                          position = position_dodge(width = 0.6)) +
          ylab("Proportion") + theme_bw() +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.text.x = element_text(size = 10),
                axis.line = element_line(colour = "black"),
                text = element_text(size = 15), 
                legend.position = "none") + xlab("Tissue") + coord_flip()
        # ggtitle("Outliers carrying rare variants in RV eGenes")
        print(p)
        dev.off()
        
        # by proportion of total outliers
        # result$mean.rare.outliers.prop = 
        #   result$mean.rare.outliers / result$mean.total.outliers
        # result$se.rare.outliers.prop = 
        #   result$se.rare.outliers / result$mean.total.outliers
        CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/outliers.prop.by.total.outliers.per.tissue.lrt.only.pdf",
                 width = 12, height = 8)
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
  }
}

# outliers analysis
# should use this!
# transformation: log2(TPM + 1) -> standardization -> corrected with covariates
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
  
  # lrt enrichment of rare variants per tissue - plot
  {
    df = rbindlist(lapply(outliers.enrichment, function(x) x[2, ]))
    df$tissue = names(outliers.enrichment)
    CairoPDF("../../materials/log_tpm.egenes.enrichment.for.rare.var.per.tissue.lrt.only.pdf", 
             width = 12, height = 8)
    p = ggplot(df, aes(x = tissue, y = relative.risk)) +
      theme_bw() + ylab('Enrichment of rare variants') + xlab("Tissue") +
      geom_abline(intercept = 1, slope = 0, linetype = "dashed") + #ylim(0.95, 1.15) +
      geom_pointrange(aes(x = tissue, ymin = min.ci, ymax = max.ci, colour = tissue)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), legend.position = "none",
            axis.text.x = element_text(size = 10), text = element_text(size = 15),
            axis.line = element_line(colour = "black")) + coord_flip() #+ 
    # geom_text(aes(label=formatC(pval, digits = 2, format = "e")), vjust = -1.5)
    print(p)
    dev.off()
  }
  
  # using gene lists from other databases
  {
    # gene list preprocessing
    {
      # overlap with clinvar
      clinvar <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/clinvar/unique_gene_symbolclinvar.vcf", header=F, data.table = F)
      g2p <- rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gene2phenotype", 
                                  full.names = T), fread))
      # gwas catalog (genes in REPORTED GENE(S) or MAPPED_GENE)
      gwas.catalog <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/gwas_catalog/reported.genes")
      gwas.catalog %>% dplyr::select("REPORTED GENE(S)") %>% distinct() -> gwas.catalog
      # omim
      omim <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/omim/mim2gene.txt")
      colnames(omim) <- c("MIM_Number", "MIM_Entry_Type", "Entrez_Gene_ID", "HGNC", "Ensembl_Gene_ID")
      # orphanet
      require(XML)
      orphanet <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/orphanet/orphanet.genes",
                        header = F)
      cancer.genes <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/GeneListsOther/cancer.genes.gold.standard.csv")
      cardio.genes <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/GeneListsOther/cardio.genes.gold.standard.csv")
      acmg <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/acmg/genes.list")
      height <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/other_genesets/generif.height.geneset.txt", header = F)
      bmi <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/other_genesets/generif.bmi.geneset.txt", header = F)
    }
    
    # annotate outliers
    {
      human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
      annotate.outliers <- function(mart, outliers) {
        genes <- outliers$gene
        entrez.hgnc <- getBM(attributes = c("entrezgene_id", "hgnc_id", "hgnc_symbol",
                                            "ensembl_gene_id"), 
                             filters = "ensembl_gene_id", mart = mart, values = genes)
        entrez.hgnc$hgnc_id <- tstrsplit(entrez.hgnc$hgnc_id, ":", keep = 2)[[1]]
        outliers$hgnc_symbol <- entrez.hgnc$hgnc_symbol[match(genes, entrez.hgnc$ensembl_gene_id)]
        outliers$hgnc_id <- entrez.hgnc$hgnc_id[match(genes, entrez.hgnc$ensembl_gene_id)]
        outliers$entrezgene_id <- entrez.hgnc$entrezgene_id[match(genes, entrez.hgnc$ensembl_gene_id)]
        return(outliers)
      }
      
      annotated.outliers <- lapply(outliers.per.gene$lrt, function(x) annotate.outliers(human, x))
      
      saveRDS(annotated.outliers, "~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/outliers/all_tissues.outliers.annotated.lrt_only.rds")
    }
    
    # gene list enrichment
    {
      gene.list.enrichment <- function(tissue, egenes, gene.list.attributes, early.output=F) 
      {
        gwas.catalog %>% filter(!(`REPORTED GENE(S)` == "" | is.na(`REPORTED GENE(S)`))) -> gwas.catalog
        clinvar %>% filter(!(V1 == "" | is.na(V1))) -> clinvar
        g2p %>% filter(!(`gene symbol` == "" | is.na(`gene symbol`))) -> g2p
        orphanet %>% filter(!(V1 == "" | is.na(V1))) -> orphanet
        omim$Entrez_Gene_ID[is.na(omim$Entrez_Gene_ID)] = 0
        omim$HGNC[omim$HGNC == ""] = "_"
        omim$Ensembl_Gene_ID[omim$Ensembl_Gene_ID == ""] = "_"
        
        # test
        uniq.gwas.catalog <- unique(gwas.catalog$`REPORTED GENE(S)`)
        uniq.clinvar <- unique(clinvar$V1)
        uniq.g2p <- unique(g2p$`gene symbol`)
        uniq.orpha <- unique(orphanet$V1)
        uniq.omim.entrez <- unique(omim$Entrez_Gene_ID)
        uniq.omim.hgnc <- unique(omim$HGNC)
        uniq.omim.ensembl <- unique(omim$Ensembl_Gene_ID)
        uniq.cancer.genes <- unique(cancer.genes$Gene)
        uniq.cardio.genes <- unique(cardio.genes$Gene)
        uniq.acmg <- unique(acmg$gene_name)
        uniq.height <- unique(height$V1)
        uniq.bmi <- unique(bmi$V1)
        
        
        if (grepl(":", gene.list.attributes$hgnc_id[[1]])) {
          gene.list.attributes$hgnc_id <- tstrsplit(gene.list.attributes$hgnc_id, ":", keep = 2)[[1]]
        }
        
        gene.list.attributes$lrt <- gene.list.attributes$gene %in% egenes
        gene.list.attributes$non.lrt <- 1 - gene.list.attributes$lrt
        
        gene.list.attributes$clinvar <- gene.list.attributes$hgnc_symbol %in% uniq.clinvar
        gene.list.attributes$g2p <- gene.list.attributes$hgnc_symbol %in% uniq.g2p
        gene.list.attributes$omim <- (gene.list.attributes$hgnc_symbol %in% uniq.omim.hgnc) | 
          (gene.list.attributes$gene %in% uniq.omim.ensembl) | 
          (gene.list.attributes$entrezgene_id %in% uniq.omim.entrez)
        gene.list.attributes$orphanet <- gene.list.attributes$hgnc_symbol %in% uniq.orpha
        gene.list.attributes$gwas <- gene.list.attributes$hgnc_symbol %in% uniq.gwas.catalog
        gene.list.attributes$cancer <- gene.list.attributes$hgnc_symbol %in% uniq.cancer.genes
        gene.list.attributes$cardio <- gene.list.attributes$hgnc_symbol %in% uniq.cardio.genes
        gene.list.attributes$acmg <- gene.list.attributes$hgnc_symbol %in% uniq.acmg
        gene.list.attributes$height <- gene.list.attributes$hgnc_symbol %in% uniq.height
        gene.list.attributes$bmi <- gene.list.attributes$hgnc_symbol %in% uniq.bmi
        ###########################################################
        
        # gene.list.attributes$lrt <- gene.list.attributes$gene %in% egenes
        # gene.list.attributes$non.lrt <- 1 - gene.list.attributes$lrt
        # gene.list.attributes$clinvar <- gene.list.attributes$hgnc_symbol %in% clinvar$V1
        # gene.list.attributes$g2p <- gene.list.attributes$hgnc_symbol %in% g2p$`gene symbol`
        # gene.list.attributes$omim <- (gene.list.attributes$hgnc_symbol %in% omim$HGNC) |
        #   (gene.list.attributes$gene %in% omim$Ensembl_Gene_ID) |
        #   (gene.list.attributes$entrezgene_id %in% omim$Entrez_Gene_ID)
        # gene.list.attributes$orphanet <- gene.list.attributes$hgnc_symbol %in% orphanet$V1
        # gene.list.attributes$gwas <- gene.list.attributes$hgnc_symbol %in% gwas.catalog$`REPORTED GENE(S)`
        gene.list.attributes$all_databases <-
          gene.list.attributes$clinvar | gene.list.attributes$g2p | gene.list.attributes$omim | gene.list.attributes$orphanet | gene.list.attributes$gwas | gene.list.attributes$cancer | gene.list.attributes$cardio | gene.list.attributes$acmg | gene.list.attributes$height
        
        if (early.output) {
          fwrite(gene.list.attributes, 
                 "~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/all.tissues.gene.list.attributes")
        }
        
        res <- data.frame(odd.ratio = numeric(), min.conf = numeric(),
                          max.conf = numeric(), pval = numeric(), database = character(),
                          method = character(), tissue = character(), 
                          odd.ratio.non.lrt = numeric(),
                          max.conf.non.lrt = numeric(),
                          min.conf.non.lrt = numeric(),
                          pval.non.lrt = numeric(),
                          stringsAsFactors = F)
        options(stringsAsFactors = FALSE)
        for (t in c("clinvar", "g2p", "omim", "orphanet", "gwas", "cancer", "cardio", "acmg", "height", "bmi", "all_databases")) {
          test <- fisher.test(table(gene.list.attributes[, c("lrt", t)]))
          test2 <- fisher.test(table(gene.list.attributes[, c("non.lrt", t)]))
          dfrow <- list(odd.ratio = as.numeric(test$estimate), 
                        max.conf = test$conf.int[2], min.conf = test$conf.int[1], 
                        pval = test$p.value, database = t, method = "lrt", 
                        tissue = tissue, 
                        odd.ratio.non.lrt = as.numeric(test2$estimate),
                        max.conf.non.lrt = test2$conf.int[2],
                        min.conf.non.lrt = test2$conf.int[1],
                        pval.non.lrt = test2$p.value)
          res <- rbind(res, dfrow)
        }
        return(res)
      }
      
      {
        library(biomaRt)
        human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
        annotate.outliers <- function(mart, outliers) {
          genes <- outliers$gene
          genes <- gsub("\\..*", "", genes)
          entrez.hgnc <- getBM(attributes = c("entrezgene_id", "hgnc_id", "hgnc_symbol",
                                              "ensembl_gene_id"), 
                               filters = "ensembl_gene_id", mart = mart, values = genes)
          entrez.hgnc$hgnc_id <- tstrsplit(entrez.hgnc$hgnc_id, ":", keep = 2)[[1]]
          outliers$hgnc_symbol <- entrez.hgnc$hgnc_symbol[match(genes, entrez.hgnc$ensembl_gene_id)]
          outliers$hgnc_id <- entrez.hgnc$hgnc_id[match(genes, entrez.hgnc$ensembl_gene_id)]
          outliers$entrezgene_id <- entrez.hgnc$entrezgene_id[match(genes, entrez.hgnc$ensembl_gene_id)]
          return(outliers)
        }
        total.expressed.genes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8//lrt.all_tissues.expresed_genes", 
                                      data.table=F, header=F)
        lrt.egenes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8//lrt.all_tissues.egenes", 
                           data.table=F, header=F)
        colnames(total.expressed.genes) = "gene"
        colnames(lrt.egenes) = "gene"
        lrt.egenes = annotate.outliers(human, lrt.egenes)
        total.expressed.genes = annotate.outliers(human, total.expressed.genes)
        saveRDS(total.expressed.genes, "~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8//total.annotated.genes.rds")
        total.expressed.genes = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8//total.annotated.genes.rds")
        all.tissue.gene.list.enrichment = gene.list.enrichment("all.tissues", lrt.egenes$gene, total.expressed.genes)
      }
      
      # use this
      CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/gene.set.enrichment.pdf")
      all.tissue.gene.list.enrichment$database = c("ClinVar", "G2P", "OMIM", "OrphaNet", "GWAS Catalog", 
                                                   "Cancer", "Cardio", "ACMG", "Height", "BMI", "All Databases")
      all.tissue.gene.list.enrichment$method = rep("LRT-q", nrow(all.tissue.gene.list.enrichment))
      all.tissue.gene.list.enrichment <- all.tissue.gene.list.enrichment[order(all.tissue.gene.list.enrichment$odd.ratio), ]
      all.tissue.gene.list.enrichment$database <- factor(all.tissue.gene.list.enrichment$database, 
                                                         levels = all.tissue.gene.list.enrichment$database)
      df = all.tissue.gene.list.enrichment[1:5, 1:5]
      df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[1, c(8:11, 1)]))
      df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[2, c(8:11, 1)]))
      df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[3, c(8:11, 1)]))
      df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[4, c(8:11, 1)]))
      df = rbind(df, as.numeric(all.tissue.gene.list.enrichment[5, c(8:11, 1)]))
      df$database[6:10] = all.tissue.gene.list.enrichment$database[1:5]
      df$type = c(rep("RV eGenes", 5), rep("Non-RV eGenes", 5))
      # p = ggplot(df %>% filter(type == "RV eGenes"), aes(x = database, y = odd.ratio)) +
      p = ggplot(all.tissue.gene.list.enrichment[c(2,3,5,6,7,8,10), ], aes(x = database, y = odd.ratio)) +
        theme_bw() + ylab('Odds ratio') + xlab("Database") +
        geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
        geom_pointrange(aes(x = database, y=odd.ratio, ymin = min.conf, ymax = max.conf, color = database), #colour = method), 
                        position = position_dodge(width = 0.6)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              # axis.text.x = element_text(angle = 90, hjust = 1),
              axis.line = element_line(colour = "black"), 
              text = element_text(size = 15),
              legend.position = "none") + coord_flip() + labs(color=NULL) + ylim(0, 5) +
        geom_text(aes(label=formatC(pval, digits = 2, format = "e")),
                  vjust = -1.5, hjust = -0.05)
      # ggtitle("Enrichment for disease-associated genes detected by LRT-q") + coord_flip()
      print(p)
      dev.off() 
      
    }
    
  }
  
  # sum over all tissues
  {
    # summary.outliers.enrichment = outliers.enrichment$AVO[, 1:6] + outliers.enrichment$AT[, 1:6] +
    # outliers.enrichment$EM[, 1:6] + outliers.enrichment$Lung[, 1:6] + outliers.enrichment$MS[, 1:6] +
    # outliers.enrichment$WB[, 1:6]
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
  
  {
    # for each gene
    {
      get.outliers.per.gene <- function(sample.tpm, sample.rare, threshold = 2) {
        if (!any(class(sample.tpm) == "data.frame")) {
          return()
        }
        # not compute for each gene, but for the whole dataset
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
        result = data.frame(
          "gene" = selected.gene,
          "total.outliers" = num.total.outliers,
          "total.controls" = num.total.controls,
          "rare.outliers" = num.rare.outliers,
          "rare.controls" = num.rare.controls,
          "common.outliers" = num.common.outliers,
          "common.controls" = num.common.controls
        )
        return(result)
      }
      {
        outliers.per.gene <- vector("list", 3)
        names(outliers.per.gene) <- c("vt", "skat", "lrt")
        outliers.per.gene[["vt"]] <- vector("list", 49)
        outliers.per.gene[["lrt"]] <- vector("list", 49)
        outliers.per.gene[["skat"]] <- vector("list", 49)
        names(outliers.per.gene[["vt"]]) <- names(indiv.with.rare.var)
        names(outliers.per.gene[["lrt"]]) <- names(indiv.with.rare.var)
        names(outliers.per.gene[["skat"]]) <- names(indiv.with.rare.var)
        for (tissue in names(indiv.with.rare.var)) {
          outliers.per.gene[["vt"]][[tissue]] <- get.outliers.per.gene(vt.rare.z[[tissue]],
                                                                       indiv.with.rare.var[[tissue]])
          outliers.per.gene[["lrt"]][[tissue]] <- get.outliers.per.gene(lrt.rare.z[[tissue]], 
                                                                        indiv.with.rare.var[[tissue]])
          outliers.per.gene[["skat"]][[tissue]] <- get.outliers.per.gene(skat.rare.z[[tissue]],
                                                                         indiv.with.rare.var[[tissue]])
        }
      }
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
    
    # outliers per tissue
    {
      CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/outliers.per.tissue.pdf")
      p <- ggplot(na.omit(result), aes(x=tissue, y=mean.rare.outliers, colour=method)) + 
        geom_pointrange(aes(x = tissue, ymin = mean.rare.outliers + 1.96 * se.rare.outliers, 
                            ymax = mean.rare.outliers - 1.96 * se.rare.outliers, 
                            colour = method), 
                        position = position_dodge(width = 0.6)) +
        ylab("Average number of outliers") + theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size = 10),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 15),
              legend.position = "top") + xlab("Tissue") + labs(colour = "Methods") + coord_flip()
      # ggtitle("Outliers carrying rare variants in RV eGenes")
      print(p)
      dev.off()
      
      CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/outliers.per.tissue.lrt.only.pdf")
      res = na.omit(result)
      p <- ggplot(res %>% filter(method == "LRT-q") %>% arrange(mean.rare.outliers) %>% 
                    mutate(tissue=factor(tissue, levels=tissue)), 
                  aes(x=tissue, y=mean.rare.outliers, colour=tissue)) + 
        geom_pointrange(aes(x = tissue, ymin = mean.rare.outliers + 1.96 * se.rare.outliers, 
                            ymax = mean.rare.outliers - 1.96 * se.rare.outliers, 
                            colour = tissue), 
                        position = position_dodge(width = 0.6)) +
        ylab("Average number of outliers") + theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size = 10),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 15), 
              legend.position = "none") + xlab("Tissue") + coord_flip()
      # ggtitle("Outliers carrying rare variants in RV eGenes")
      print(p)
      dev.off()
      
      # by proportion of total samples
      result$mean.rare.outliers.prop = 
        result$mean.rare.outliers / (result$mean.total.outliers + result$mean.total.controls)
      result$se.rare.outliers.prop = 
        result$se.rare.outliers / (result$mean.total.outliers + result$mean.total.controls)
      CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/outliers.prop.by.total.samples.per.tissue.lrt.only.pdf",
               width = 12, height = 8)
      p <- ggplot(na.omit(result) %>% filter(method == "LRT-q") %>% arrange(mean.rare.outliers.prop) %>% 
                    mutate(tissue=factor(tissue, levels=tissue)), 
                  aes(x=tissue, y=mean.rare.outliers.prop, colour=tissue)) + 
        geom_pointrange(aes(x = tissue, ymin = mean.rare.outliers.prop + 1.96 * se.rare.outliers.prop, 
                            ymax = mean.rare.outliers.prop - 1.96 * se.rare.outliers.prop, 
                            colour = tissue), 
                        position = position_dodge(width = 0.6)) +
        ylab("Proportion") + theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size = 10),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 15), 
              legend.position = "none") + xlab("Tissue") + coord_flip()
      # ggtitle("Outliers carrying rare variants in RV eGenes")
      print(p)
      dev.off()
      
      # by proportion of total outliers
      # result$mean.rare.outliers.prop = 
      #   result$mean.rare.outliers / result$mean.total.outliers
      # result$se.rare.outliers.prop = 
      #   result$se.rare.outliers / result$mean.total.outliers
      CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/outliers.prop.by.total.outliers.per.tissue.lrt.only.pdf",
               width = 12, height = 8)
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
}

# generate curve for effect sizes
{
  library(Cairo)
  library(ggpubr)
  library(data.table)
  x = seq(0.0005, 0.05, 0.00001)
  a = seq(0.3, 1.5, 0.3)
  y = as.matrix(a) %*% abs(log10(x))
  data = data.frame("Effect size" = c(y[1, ], y[2, ], y[3, ], y[4, ], y[5, ]),
                    "MAF" = rep(x, 5),
                    "a" = c(rep("a = 0.3", 4951), rep("a = 0.6", 4951),
                            rep("a = 0.9", 4951), rep("a = 1.2", 4951),
                            rep("a = 1.5", 4951)))
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/fig.s.maf.vs.effect.size.pdf")
  p = ggline(data, x = "MAF", y = "Effect.size", color = "a", ylab = "Effect size",
             plot_type = "l", numeric.x.axis = T) + rremove("legend.title")
  print(p)
  dev.off()
}

# decision boundary
{
  library(data.table)
  library(ggpubr)
  library(Cairo)
  library(dplyr)
  setwd("/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/power/decision_boundary")
  
  c05 = rbindlist(lapply(dir(pattern = "c05", full.names = T), fread))
  c1 = rbindlist(lapply(dir(pattern = "c1", full.names = T), fread))
  
  c05 = na.omit(as.data.frame(c05))
  c05[, 7:14] = if_else(c05[, 7:14] < 0.05, "Significant", "Not significant")
  colnames(c05)[11:14] = c("SKATO", "ACATV", "ACATO", "LRTq")
  c05$truth = if_else(c05$causal1 == 0 & c05$causal2 == 0, 
                      "CE = 0", "CE > 0")
  
  # use transparency to avoid overplotting
  ggplot(c05, aes(x=reg.t.stat1, y=reg.t.stat2, color=LRTq, shape=truth)) +
    geom_point(alpha=0.3)
  # use sampling to avoid overplotting
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.LRTq.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", 
                   color = "LRTq", shape = "truth", ellipse.alpha = 0.3)
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = quote("LRT-q ("~c[i]~"= 0.5)"))
    dev.off()
    }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.VT.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", 
                   color = "VT", shape = "truth")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = quote("VT ("~c[i]~"= 0.5)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.SKATO.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", 
                   color = "SKATO", shape = "truth")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("SKAT-O ("~c[i]~"= 0.5)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.WSS.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:10000, ], x = "reg.t.stat1", y = "reg.t.stat2", color = "WSS")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("WSS ("~c[i]~"= 0.5)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.BURDEN.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:10000, ], x = "reg.t.stat1", y = "reg.t.stat2", color = "BURDEN")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("BURDEN ("~c[i]~"= 0.5)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.CMC.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", 
                   color = "CMC", shape = "truth")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("CMC ("~c[i]~"= 0.5)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.ACATV.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", color = "ACATV", shape = "truth")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("ACAT-V ("~c[i]~"= 0.5)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.ACATO.c05.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c05[1:2500, ], x = "reg.t.stat1", y = "reg.t.stat2", color = "ACATO", shape = "truth")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("ACAT-O ("~c[i]~"= 0.5)"))
    dev.off()
  }
  
  c1 = na.omit(as.data.frame(c1))
  c1[, 7:ncol(c1)] = if_else(c1[, 7:ncol(c1)] < 0.05, "Significant", "Not significant")
  colnames(c1)[11:14] = c("SKATO", "ACATV", "ACATO", "LRTq")
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.LRTq.c1.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c1[1:10000, ], x = "reg.t.stat1", y = "reg.t.stat2", color = "LRTq")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("LRT-q ("~c[i]~"= 1.0)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.VT.c1.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c1[1:10000, ], x = "reg.t.stat1", y = "reg.t.stat2", color = "VT")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("VT ("~c[i]~"= 1.0)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.SKATO.c1.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c1[1:10000, ], x = "reg.t.stat1", y = "reg.t.stat2", color = "SKATO")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("SKAT-O ("~c[i]~"= 1.0)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.CMC.c1.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c1[1:10000, ], x = "reg.t.stat1", y = "reg.t.stat2", color = "CMC")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("CMC ("~c[i]~"= 1.0)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.ACATV.c1.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c1, x = "var1", y = "var2", color = "ACATV")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("ACAT-V ("~c[i]~"= 1.0)"))
    dev.off()
  }
  {
    folder = "/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/materials/"
    pic.name = paste(folder, "decision_boundary.ACATO.c1.pdf", sep = "")
    CairoPDF(pic.name)
    p <- ggscatter(c1, x = "var1", y = "var2", color = "ACATO")
    ggpar(p, legend.title = "", font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
          font.y = 13, xlab = "t-score of variant 1", ylab = "t-score of variant 2",
          title = bquote("ACAT-O ("~c[i]~"= 1.0)"))
    dev.off()
  }
  
}

# egenes vs sample size
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
  
  p = ggscatter(egene.count, x="sample.size", y="lrt.egenes", color="tissue")
  CairoPDF("projects/rare_variant/materials/egenes.vs.sample.size.pdf")
  ggpar(p, legend="none", ylab="Number of eGenes detected by LRT-q",
        xlab="Number of samples", font.tickslab = 13, font.legend = 13, 
        font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  
  p = ggscatter(egene.count, x="sample.size", y="lrt.egenes.per.gene", color="tissue")
  CairoPDF("projects/rare_variant/materials/egenes.per.gene.vs.sample.size.pdf")
  ggpar(p, legend="none", ylab="Number of eGenes detected by LRT-q / total expressed genes",
        xlab="Number of samples", font.tickslab = 13, font.legend = 13, 
        font.main = 13, font.x = 13, font.y = 13)
  dev.off()
  
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
  
  CairoPDF("projects/rare_variant/materials/novel.rv.egenes.each.tissue.pdf", 
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

# tissue sharing
{
  library(cowplot)
  library(ggplot2)
  library(ggpubr)
  library(ggcorrplot)
  library(reshape2)
  library(Cairo)
  library(dplyr)
  # shared.q = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/qvals.shared.by.weights.tissues.rds")
  reorder_cormat <- function(cormat) {
    # Use correlation between variables as distance
    na.idx <- apply(cormat, 1, function(x) all(is.na(x)))
    cormat <- cormat[!na.idx, !na.idx]
    cormat[is.na(cormat)] <- 0.20
    
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
    return(cormat)
  }
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  plot.sharing.mat <- function(cormat) {
    # cormat <- reorder_cormat(cormat)
    upper_tri = get_lower_tri(cormat)
    melted_cormat <- melt(as.matrix(upper_tri), na.rm = TRUE)
    melted_cormat <- melt(as.matrix(cormat), na.rm = T)
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+ 
      scale_fill_distiller(palette = "YlOrBr", direction = 1) +
      # scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.20,
      #                      limit = c(0,1), space = "Lab",
      #                      name="Sharing\nfraction") +
      theme_minimal() + # minimal theme
      #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       #     panel.background = element_blank(), axis.line = element_line(colour = "black"))
      theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                       size = 5, hjust = 1),
            axis.text.y = element_text(vjust = 1, 
                                       size = 5, hjust = 1),
            axis.title.x=element_blank(), axis.title.y=element_blank())+
      coord_fixed()
    return(ggheatmap)
  }
  
  # ggcorrplot(cor(shared.q$MAFxCADD, method = "spearman"), hc.order = T, 
  #            outline.color = "white", type = "upper", lab = F, 
  #            title = "Tissue correlation", 
  #            legend.title = "Spearman's coefficient", tl.cex = 8)
  
  # RV eGenes
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/tissue.sharing.corr.matrix.pdf")
  corr.matrix <- readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing//tissue_sharing.fraction_shared_egenes.lrt.0.10.rds")
  # plot.sharing.mat(reorder_cormat(corr.matrix))
  plot.sharing.mat(corr.matrix)
  dev.off()
  
  # clustered
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/tissue.sharing.corr.matrix.clustered.pdf")
  corr.matrix <- readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing//tissue_sharing.fraction_shared_egenes.lrt.0.10.rds")
  plot.sharing.mat(reorder_cormat(corr.matrix))
  # plot.sharing.mat(corr.matrix)
  dev.off()
  
  # CV eGenes
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/tissue.sharing.corr.matrix.cv.egenes.pdf")
  # corr.matrix <- readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.fraction_shared_egenes.cv.egenes.q.thres.0.1.rds")
  corr.matrix <- readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.fraction_shared_egenes.cv.egenes.first.200.genes.rds")
  # plot.sharing.mat(reorder_cormat(corr.matrix))
  plot.sharing.mat(corr.matrix)
  dev.off()
  
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/tissue.sharing.corr.matrix.cv.egenes.clustered.pdf")
  # corr.matrix <- readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.fraction_shared_egenes.cv.egenes.q.thres.0.1.rds")
  corr.matrix <- readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.fraction_shared_egenes.cv.egenes.first.200.genes.rds")
  plot.sharing.mat(reorder_cormat(corr.matrix))
  # plot.sharing.mat(corr.matrix)
  dev.off()
  
  
  # ggcorrplot(corr.matrix, lab = F, outline.color = "white", type="upper", title = "Tissue correlation",
             # legend.title = "Correlation", tl.cex = 8)
  
  # shared egenes by number of tissues
  shared.matrix <- readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.rv.egenes.rds")
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.shared.rv.egenes.num.tissues.pdf")
  p = ggbarplot(shared.matrix, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
            ylab = "Proportion of shared RV eGenes", fill = "black")
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13)
  dev.off()
  
  # fraction of tissue-specific egenes and universal egenes
  frac.matrix <- readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.fraction.universal.tissue_specific.egenes.lrt.fdr.0.05.rds")
  frac.matrix$tissue = rownames(frac.matrix)
  frac.matrix$short.tissue = gsub("Brain.*", "Brain tissues", frac.matrix$tissue)
  frac.matrix$short.tissue = gsub("Skin.*", "Skin tissues", frac.matrix$short.tissue)
  frac.matrix$short.tissue = gsub("Artery.*", "Artery tissues", frac.matrix$short.tissue)
  frac.matrix$short.tissue = gsub("Adipose.*", "Adipose tissues", frac.matrix$short.tissue)
  frac.matrix$short.tissue = gsub("Heart.*", "Heart tissues", frac.matrix$short.tissue)
  frac.matrix$short.tissue = gsub("Colon.*", "Colon tissues", frac.matrix$short.tissue)
  frac.matrix$short.tissue = gsub("Esophagus.*", "Esophagus tissues", frac.matrix$short.tissue)
  frac.matrix$short.tissue[nrow(frac.matrix)] = "Universal"
  frac.matrix$type = c(rep("Tissue-specific", nrow(frac.matrix)-1), "Universal")
  
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.tissue.specific.universal.egenes.by.tissue.pdf")
  frac.matrix %>% group_by(short.tissue, type) %>% dplyr::summarise(frac=sum(prop)) %>% 
    arrange(desc(frac)) %>%
    ggbarplot(x="short.tissue", y="frac", fill="type") + rotate_x_text(45) -> p
  ggpar(p, font.tickslab = 11, font.legend = 13, font.main = 13,
        ylab = "Fraction of total RV eGenes", xlab = "Tissue", legend = "none")
  dev.off()
  
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.tissue.specific.universal.egenes.pdf")
  simple.frac.matrix = data.frame(
    "Type" = c("Tissue-specific RV eGenes", "Universal RV eGenes"),
    "Fraction" = c(sum(frac.matrix$prop[-nrow(frac.matrix)]), frac.matrix$prop[nrow(frac.matrix)])
    )
  p = ggbarplot(simple.frac.matrix, x="Type", y="Fraction", fill = "Type", 
                legend="none", ylab = "Fraction of toal RV eGenes")
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13)
  dev.off()
  
  # number of tissue-specific cv and rv egenes
  cv.ts.egenes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/all_tissues.tissue.specific.cv.egenes", data.table = F)
  rv.ts.egenes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/all_tissues.tissue.specific.rv.egenes", data.table = F)
  ts.egenes = cv.ts.egenes[, c(1, 3, 4)]
  colnames(ts.egenes) = c("Tissue", "cv.ts", "cv.ts.prop")
  ts.egenes$lrt.ts = rv.ts.egenes$`LRT-q`
  ts.egenes$lrt.ts.prop = ts.egenes$lrt.ts / egene.count$lrt.egenes
  df = melt(ts.egenes, id.vars = "Tissue", 
            measure.vars = c("cv.ts.prop", "lrt.ts.prop"))
  df$value = log2(df$value)
  library(plyr)
  df$variable = revalue(df$variable, 
                        c("cv.ts.prop"="CV eGenes", "lrt.ts.prop"="RV eGenes"))
  colnames(df)[2] = "Type"
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/tissue_specific.rv.cv.egenes.pdf")
  p = ggboxplot(df, x="Type", y="value", fill="Type",
            ylab="Proportion of tissue-specific eGenes (Log-transformed)")
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, 
        font.x = 13, font.y = 13, legend = "none")
  dev.off()
}

# tissue sharing - combine RV eGenes and CV eGenes
{
  library(ggpubr)
  library(data.table)
  library(Cairo)
  
  # test - not useful
  {
    rv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/combine_rv_cv_egenes/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.rv.egenes.rds")
    cv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/combine_rv_cv_egenes/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.cv.egenes.rds")
    df = rbind(rv, cv)
    df$type = c(rep("RV eGenes", nrow(rv)), rep("CV eGenes", nrow(cv)))
    
    CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.shared.cv.egenes.rv.egenes.num.tissues.pdf")
    p = ggbarplot(df, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
                  ylab = "Proportion of shared RV eGenes", fill = "type", 
                  position = position_dodge())
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, legend.title = "",
          font.x = 13, font.y = 13, title = "Tissues with more than 300 RV eGenes")
    dev.off()
    
    group.limit = 4
    rv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/combine_rv_cv_egenes/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.rv.egenes.rds")
    cv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/combine_rv_cv_egenes/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.cv.egenes.rds")
    df2 = rbind(rv[1:group.limit, ], cv[1:group.limit, ])
    df2$num.tissues = as.character(df2$num.tissues)
    df2$type = c(rep("RV eGenes", group.limit), rep("CV eGenes", group.limit))
    df2 = rbind(df2, c(paste(">", group.limit, sep=" "), 
                       1-sum(rv$prop.shared[1:group.limit]), "RV eGenes"))
    df2 = rbind(df2, c(paste(">", group.limit, sep=" "), 
                       1-sum(cv$prop.shared[1:group.limit]), "CV eGenes"))
    df2$prop.shared = as.numeric(df2$prop.shared)
    CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.shared.cv.egenes.rv.egenes.num.tissues.grouped.200.pdf")
    p = ggbarplot(df2, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
                  ylab = "Proportion of shared RV eGenes", fill = "type", 
                  position = position_dodge())
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, legend.title = "",
          font.x = 13, font.y = 13, title = "Tissues with more than 200 RV eGenes",
          font.title = 13)
    dev.off()
    
    group.limit = 4
    rv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/combine_rv_cv_egenes/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.100.rv.egenes.rds")
    cv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/combine_rv_cv_egenes/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.100.cv.egenes.rds")
    df2 = rbind(rv[1:group.limit, ], cv[1:group.limit, ])
    df2$num.tissues = as.character(df2$num.tissues)
    df2$type = c(rep("RV eGenes", group.limit), rep("CV eGenes", group.limit))
    df2 = rbind(df2, c(paste(">", group.limit, sep=" "), 
                       1-sum(rv$prop.shared[1:group.limit]), "RV eGenes"))
    df2 = rbind(df2, c(paste(">", group.limit, sep=" "), 
                       1-sum(cv$prop.shared[1:group.limit]), "CV eGenes"))
    df2$prop.shared = as.numeric(df2$prop.shared)
    CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.shared.cv.egenes.rv.egenes.num.tissues.grouped.100.pdf")
    p = ggbarplot(df2, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
                  ylab = "Proportion of shared RV eGenes", fill = "type", 
                  position = position_dodge())
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, legend.title = "",
          font.x = 13, font.y = 13, title = "Tissues with more than 100 RV eGenes",
          font.title = 13)
    dev.off()
  }
  
  group.limit = 4
  rv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.rv.egenes.rds")
  cv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing//tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.cv.egenes.rds")
  df3 = rbind(rv[1, ], cv[1, ])
  df3$num.tissues = as.character(df3$num.tissues)
  df3$type = c(rep("RV eGenes", 1), rep("CV eGenes", 1))
  df3 = rbind(df3, c(paste("2 -", group.limit, sep=" "), 
                     sum(rv$prop.shared[2:group.limit]), "RV eGenes"))
  df3 = rbind(df3, c(paste("2 -", group.limit, sep=" "), 
                     sum(cv$prop.shared[2:group.limit]), "CV eGenes"))
  df3 = rbind(df3, c(paste(">", group.limit, sep=" "), 
                     1-sum(rv$prop.shared[1:group.limit]), "RV eGenes"))
  df3 = rbind(df3, c(paste(">", group.limit, sep=" "), 
                     1-sum(cv$prop.shared[1:group.limit]), "CV eGenes"))
  df3$prop.shared = as.numeric(df3$prop.shared)
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.shared.cv.egenes.rv.egenes.num.tissues.grouped.200.pdf")
  p = ggbarplot(df3, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
                ylab = "Proportion of shared eGenes", fill = "type", 
                position = position_dodge())
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, legend.title = "",
        font.x = 13, font.y = 13, title = "Tissues with more than 200 RV eGenes",
        font.title = 13)
  dev.off()
  
  group.limit = 4
  rv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.100.rv.egenes.rds")
  cv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.100.cv.egenes.rds")
  df3 = rbind(rv[1, ], cv[1, ])
  df3$num.tissues = as.character(df3$num.tissues)
  df3$type = c(rep("RV eGenes", 1), rep("CV eGenes", 1))
  df3 = rbind(df3, c(paste("2 -", group.limit, sep=" "), 
                     sum(rv$prop.shared[2:group.limit]), "RV eGenes"))
  df3 = rbind(df3, c(paste("2 -", group.limit, sep=" "), 
                     sum(cv$prop.shared[2:group.limit]), "CV eGenes"))
  df3 = rbind(df3, c(paste(">", group.limit, sep=" "), 
                     1-sum(rv$prop.shared[1:group.limit]), "RV eGenes"))
  df3 = rbind(df3, c(paste(">", group.limit, sep=" "), 
                     1-sum(cv$prop.shared[1:group.limit]), "CV eGenes"))
  df3$prop.shared = as.numeric(df3$prop.shared)
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.shared.cv.egenes.rv.egenes.num.tissues.grouped.100.pdf")
  p = ggbarplot(df3, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
                ylab = "Proportion of shared eGenes", fill = "type", 
                position = position_dodge())
  ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, legend.title = "",
        font.x = 13, font.y = 13, title = "Tissues with more than 100 RV eGenes",
        font.title = 13)
  dev.off()
  
  # plot all numbers of tissues
  {
    # tissues with >= 200 RV eGenes
    rv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.rv.egenes.rds")
    cv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing//tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.200.cv.egenes.rds")
    df3 = rbind(rv, cv)
    df3$type = c(rep("RV eGenes", nrow(rv)), rep("CV eGenes", nrow(cv)))
    df3$prop.shared = as.numeric(df3$prop.shared)
    CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.shared.cv.egenes.rv.egenes.num.tissues.grouped.200.all_num_tissues.pdf")
    p = ggbarplot(df3, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
                  ylab = "Proportion of shared eGenes", fill = "type", 
                  position = position_dodge())
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, legend.title = "",
          font.x = 13, font.y = 13, title = "Tissues with more than 200 RV eGenes",
          font.title = 13)
    dev.off()
    
    # tissues with >= 100 RV eGenes
    rv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.100.rv.egenes.rds")
    cv = readRDS("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/tissue_sharing/tissue_sharing.prop_shared_egenes.by.num_tissues.lrt.q.thres.0.05.num.egene.thres.100.cv.egenes.rds")
    df3 = rbind(rv, cv)
    df3$type = c(rep("RV eGenes", nrow(rv)), rep("CV eGenes", nrow(cv)))
    df3$prop.shared = as.numeric(df3$prop.shared)
    CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/prop.shared.cv.egenes.rv.egenes.num.tissues.grouped.100.all_num_tissues.pdf")
    p = ggbarplot(df3, x="num.tissues", y="prop.shared", xlab = "Number of tissues", 
                  ylab = "Proportion of shared eGenes", fill = "type", 
                  position = position_dodge())
    ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, legend.title = "",
          font.x = 13, font.y = 13, title = "Tissues with more than 100 RV eGenes",
          font.title = 13)
    dev.off()
    
  }
}

# number of rare variants within 1mb and 20kb from TSS
{
  library(data.table)
  library(Cairo)
  library(ggpubr)
  one_mb = rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/rare_var_per_gene/v8/1mb/", full.names = T), fread))
  twenty_kb = rbindlist(lapply(dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/rare_var_per_gene/v8/20kb/", full.names = T), fread))
  one_mb$dis = "1 Mb"
  twenty_kb$dis = "20 Kb"
  CairoPDF("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/rare_var_per_gene.pdf")
  p = ggboxplot(rbind(one_mb, twenty_kb), x = "dis", y = "V2", fill = "dis",
            xlab = "Distance from TSS", ylab = "Number of rare variants per gene")
  ggpar(p, legend.title = "Distance from TSS")
  dev.off()
}

# GTEx v7 v8 sample sizes
{
  require(data.table)
  require(Cairo)
  require(ggpubr)
  sample.sizes = fread("/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/sample_sizes/v7.v8.sample.size.comparison.list.txt")
  colnames(sample.sizes) = c("Tissue", "GTEx v7", "GTEx v8")
  df = melt(sample.sizes, id.vars = "Tissue")
  colnames(df) = c("Tissue", "Dataset", "Sample.size")
  ggbarplot(df, fill = "Dataset", position = position_dodge(0.9), 
            x = "Tissue", y = "Sample.size", ylab = "Sample size") + rotate_x_text(90)
}

# GTEx v8 tissues
{
  require(data.table)
  require(Cairo)
  require(ggpubr)
  egene.count1 = fread("/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/gtex8/egene.counts.p1.csv")
  egene.count2 = fread("/Users/avallonking/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/gtex8/egene.counts.csv")
  egene.count = rbind(egene.count1, egene.count2)
  
}

# Nahyun genes
{
  library(data.table)
  library(dplyr)
  
  files = dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/nahyun/send", recursive = T,
              pattern = "disease_gene.txt", full.names = T)
  gene.v7 = lapply(files, fread)
  names(gene.v7) = tstrsplit(basename(files), "\\.")[[1]]
  
  files2 = dir("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/nahyun/replicated", full.names = T)
  replicated.genes = lapply(files2, fread)
  names(replicated.genes) = tstrsplit(basename(files2), "\\.")[[1]]
  
  gene.v8 = gene.v7
  for (i in names(gene.v7)) {
    gene.v8[[i]] = gene.v7[[i]] %>% filter(EnsemblID %in% replicated.genes[[i]]$replicated.disease.genes)
  }
  
  sort(table(unlist(lapply(gene.v8, function(x) unique(x$EnsemblID)))))
}

# 100 perm vs 100k perm
{
  library(data.table)
  library(Cairo)
  library(ggpubr)
  library(VennDiagram)
  library(RColorBrewer)
  setwd("~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/")
  perm100.egenes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.wb.perm100.egenes.list")
  perm10k.egenes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.wb.perm10k.egenes.list")
  perm1k.egenes = fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.wb.perm1k.egenes.list")
  myCol <- brewer.pal(3, "Pastel2")
  venn.diagram(list(perm10k.egenes$perm10k.egenes, perm1k.egenes$perm1k.egenes, perm100.egenes$perm100.egenes), 
               filename = "adaptive.perm.wb.png", 
               category.names = c("100k permutations", "1k permutations", "100 permutations"), 
               output = T, imagetype = "png", fill = myCol,
               height = 480, width = 480, lwd = 2, lty = 'blank',
               # number
               cex = .3,
               fontface = "bold",
               fontfamily = "sans",
               # label
               cat.cex = 0.2,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.pos = c(0, 0, 0),
               cat.dist = c(-0.05, 0.03, 0.01),
               cat.fontfamily = "sans"
               )
}

# Venn diagram
{
  require(data.table)
  require(VennDiagram)
  library(RColorBrewer)

  lrt.egenes <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.all_tissues.egenes", header = F)$V1
  vt.egenes <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/vt.all_tissues.egenes")$V1
  skat.egenes <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/skat.all_tissues.egenes")$V1
  acat.egenes <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/acat.all_tissues.egenes")$V1
  
  myCol <- brewer.pal(4, "Pastel2")
  venn.diagram(
    x = list(lrt.egenes, vt.egenes, skat.egenes, acat.egenes),
    category.names = c("LRT-q" , "VT", "SKAT-O", "ACAT-O"),
    filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.all.tissues.egenes.png",
    output=T,
    
    # Output features
    imagetype="png" ,
    height = 480 ,
    width = 480 ,
    resolution = 300,
    compression = "lzw",

    # Circles
    lwd = 1,
    col = myCol,
    fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),

    # Numbers
    cex = .4,
    # fontface = "bold",
    fontfamily = "sans",

    # Set names
    cat.cex = 0.4,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    # cat.pos = c(-27, 23, 135, 13),
    # cat.dist = c(0.055, 0.035, 0.085, 0.02),
    cat.fontfamily = "sans",
    cat.col = myCol
    # rotation = 1
  )
  
  # tissues with > 100 egenes
  {
    rv.egenes.shared.count <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.rv.egenes.shared.count")
    cv.egenes.shared.count <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.cv.egenes.shared.count")
    myCol <- brewer.pal(6, "Pastel2")[1:2]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq==1) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq==1) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.1tissue.png",
      output = T,
      main = "Tissue-specific eGenes",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
    
    myCol <- brewer.pal(6, "Pastel2")[3:4]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq >= 2 & Freq <= 4) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq >= 2 & Freq <= 4) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.2-4tissues.png",
      output = T,
      main = "eGenes shared in 2 - 4 tissues",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
    
    myCol <- brewer.pal(8, "Set1")[c(1,2)]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq > 4) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq > 4) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.4_or_more_tissues.png",
      output = T,
      main = "eGenes shared in more than 4 tissues",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(-180, 180),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
  }
  
  # tissues with > 200 egenes
  {
    rv.egenes.shared.count <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.rv.egenes.shared.count.tissue_200egenes")
    cv.egenes.shared.count <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.cv.egenes.shared.count.tissue_200egenes")
    myCol <- brewer.pal(6, "Pastel2")[1:2]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq==1) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq==1) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.1tissue.tissue_200egenes.png",
      output = T,
      main = "Tissue-specific eGenes",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
    
    myCol <- brewer.pal(6, "Pastel2")[3:4]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq >= 2 & Freq <= 4) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq >= 2 & Freq <= 4) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.2-4tissues.tissue_200egenes.png",
      output = T,
      main = "eGenes shared in 2 - 4 tissues",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
    
    myCol <- brewer.pal(8, "Set1")[c(1,2)]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq > 4) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq > 4) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.4_or_more_tissues.tissue_200egenes.png",
      output = T,
      main = "eGenes shared in more than 4 tissues",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol,
      inverted = T
    )
  }
  
  # tissues with > 100 egenes, MAF < 1%
  {
    rv.egenes.shared.count <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.rv.egenes.shared.count.tissue_100egenes.maf01")
    cv.egenes.shared.count <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.cv.egenes.shared.count.tissue_100egenes.maf01")
    myCol <- brewer.pal(6, "Pastel2")[1:2]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq==1) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq==1) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.1tissue.maf01.png",
      output = T,
      main = "Tissue-specific eGenes",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
    
    myCol <- brewer.pal(6, "Pastel2")[3:4]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq >= 2 & Freq <= 4) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq >= 2 & Freq <= 4) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.2-4tissues.maf01.png",
      output = T,
      main = "eGenes shared in 2 - 4 tissues",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
    
    myCol <- brewer.pal(8, "Set1")[c(1,2)]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq > 4) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq > 4) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.4_or_more_tissues.maf01.png",
      output = T,
      main = "eGenes shared in more than 4 tissues",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(-180, 180),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
  }
  
  # tissues with > 200 egenes, MAF < 1%
  {
    rv.egenes.shared.count <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.rv.egenes.shared.count.tissue_200egenes.maf01")
    cv.egenes.shared.count <- fread("~/Documents/scires/Jae-Hoon/projects/rare_variant/gtex/all_tissues/gtex8/lrt.cv.egenes.shared.count.tissue_200egenes.maf01")
    myCol <- brewer.pal(6, "Pastel2")[1:2]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq==1) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq==1) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.1tissue.tissue_200egenes.maf01.png",
      output = T,
      main = "Tissue-specific eGenes",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
    
    myCol <- brewer.pal(6, "Pastel2")[3:4]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq >= 2 & Freq <= 4) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq >= 2 & Freq <= 4) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.2-4tissues.tissue_200egenes.maf01.png",
      output = T,
      main = "eGenes shared in 2 - 4 tissues",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol
      # rotation = 1
    )
    
    myCol <- brewer.pal(8, "Set1")[c(1,2)]
    venn.diagram(
      x = list(rv.egenes.shared.count %>% filter(Freq > 4) %>% pull(gene), 
               cv.egenes.shared.count %>% filter(Freq > 4) %>% pull(gene)
      ),
      category.names = c("RV eGenes", "CV eGenes"),
      filename = "~/Documents/scires/Jae-Hoon/projects/rare_variant/materials/venn.diagram.egenes.tissue.sharing.4_or_more_tissues.tissue_200egenes.maf01.png",
      output = T,
      main = "eGenes shared in more than 4 tissues",
      main.cex = 0.5,
      main.fontfamily = "sans",
      
      # Output features
      imagetype="png" ,
      height = 480 ,
      width = 480 ,
      resolution = 300,
      compression = "lzw",
      
      # Circles
      lwd = 1,
      col = myCol,
      fill = unlist(lapply(myCol, function(x) alpha(x, 0.3))),
      
      # Numbers
      cex = .4,
      # fontface = "bold",
      fontfamily = "sans",
      
      # Set names
      cat.cex = 0.4,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.pos = c(0, 0),
      cat.dist = c(0.05, 0.05),
      cat.fontfamily = "sans",
      cat.col = myCol,
      inverted = T
    )
  }
}