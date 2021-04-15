{
    require(data.table)
    require(dplyr)
    require(tidyr)
    require(ggpubr)
    require(Cairo)
    require(ACAT)
    setwd("../data/simulation/")

    # type I error
    # Table 1
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
        res = rbindlist(lapply(dir("./type_I_error/", full.names = T), fread))
        res$`ACAT-O` = apply(res, 1, function(x) ACAT(c(x[4], x[5], x[6])))
        fdr = apply(res, 2, typeIerror.check)
        fdr.table = cbind(fdr[[1]]$alpha.level, fdr[[1]]$fdr, fdr[[2]]$fdr, fdr[[3]]$fdr, 
                          fdr[[4]]$fdr, fdr[[5]]$fdr, fdr[[6]]$fdr, fdr[[7]]$fdr, 
                          fdr[[8]]$fdr)
        colnames(fdr.table) = c("alpha", names(fdr))
        fdr.table
    }
}

