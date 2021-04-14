{
    require(data.table)
    require(dplyr)
    require(tidyr)
    require(ggpubr)
    require(Cairo)
    require(ACAT)
    setwd("~/Documents/scires/Jae-Hoon/projects/rare_variant/simulation/")

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
        # Fig 1A
        CairoPDF("../materials/acat.causal0.10.protective50.pdf")
        power.transform %>% filter(causal == "0.10" & protective.rate == "0.50") %>%
            mutate(effect = as.numeric(a) * 3.3) %>%
            ggline(x = "effect", y = "power", color = "method",  shape = "method", 
                   xlab = "Effect size upper bound", ylab = "Power", legend.title = "Method",
                   title = expression(paste(alpha, "\ = 0.05, causal ratio = 10%"))) -> p
            ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
                  font.y = 13)
            dev.off()

        # Fig 1B
        CairoPDF("../materials/acat.effect.size.0.99.protective50.pdf")
        power.transform %>% filter(a == "0.3" & protective.rate == "0.50") %>% 
            ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
                   ylab = "Power", legend.title = "Method", shape = "method", 
                   title = expression(paste(alpha, "\ = 0.05, effect size"<=0.99))) -> p
        ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
              font.y = 13)
        dev.off()

        # Fig 1C
        CairoPDF("../materials/acat.effect.size.2.97.protective50.pdf")
        power.transform %>% filter(a == "0.9" & protective.rate == "0.50") %>% 
            ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
                   ylab = "Power", legend.title = "Method", shape = "method", 
                   title = expression(paste(alpha, "\ = 0.05, effect size"<=2.97))) -> p
        ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
              font.y = 13)
        dev.off()

        # Fig 1D
        CairoPDF("../materials/acat.effect.size.4.95.protective50.pdf")
        power.transform %>% filter(a == "1.5" & protective.rate == "0.50") %>% 
            ggline(x = "causal", y = "power", color = "method", xlab = "Causal ratio",
                   ylab = "Power", legend.title = "Method", shape = "method", 
                   title = expression(paste(alpha, "\ = 0.05, effect size"<=4.95))) -> p
        ggpar(p, font.tickslab = 13, font.legend = 13, font.main = 13, font.x = 13, 
              font.y = 13)
        dev.off()
    }
}
