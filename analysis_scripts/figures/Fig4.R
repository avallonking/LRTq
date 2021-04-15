# Fig 4A: line 77-93
# Fig 4B: line 260-291 
{
  require(data.table)
  require(ggpubr)
  require(tidyr)
  require(dplyr)
  require(Cairo)
  setwd("./data/")
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
  
  # Fig 4A
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
      clinvar <- fread("./data/gene_sets/clinvar/unique_gene_symbolclinvar.vcf", header=F, data.table = F)
      g2p <- rbindlist(lapply(dir("./data/gene_sets/gene2phenotype", 
                                  full.names = T), fread))
      # gwas catalog (genes in REPORTED GENE(S) or MAPPED_GENE)
      gwas.catalog <- fread("./data/gene_sets/gwas_catalog/reported.genes")
      gwas.catalog %>% dplyr::select("REPORTED GENE(S)") %>% distinct() -> gwas.catalog
      # omim
      omim <- fread("./data/gene_sets/omim/mim2gene.txt")
      colnames(omim) <- c("MIM_Number", "MIM_Entry_Type", "Entrez_Gene_ID", "HGNC", "Ensembl_Gene_ID")
      # orphanet
      require(XML)
      orphanet <- fread("./data/gene_sets/orphanet/orphanet.genes",
                        header = F)
      cancer.genes <- fread("./data/gene_sets/GeneListsOther/cancer.genes.gold.standard.csv")
      cardio.genes <- fread("./data/gene_sets/GeneListsOther/cardio.genes.gold.standard.csv")
      acmg <- fread("./data/gene_sets/acmg/genes.list")
      height <- fread("./data/gene_sets/other_genesets/generif.height.geneset.txt", header = F)
      bmi <- fread("./data/gene_sets/other_genesets/generif.bmi.geneset.txt", header = F)
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
      
      saveRDS(annotated.outliers, "./data/gene_sets/all_tissues.outliers.annotated.lrt_only.rds")
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
        
        gene.list.attributes$all_databases <-
          gene.list.attributes$clinvar | gene.list.attributes$g2p | gene.list.attributes$omim | gene.list.attributes$orphanet | gene.list.attributes$gwas | gene.list.attributes$cancer | gene.list.attributes$cardio | gene.list.attributes$acmg | gene.list.attributes$height
        
        if (early.output) {
          fwrite(gene.list.attributes, 
                 "./data/all.tissues.gene.list.attributes")
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
        total.expressed.genes = fread("./data/lrt.all_tissues.expresed_genes", 
                                      data.table=F, header=F)
        lrt.egenes = fread("./data/lrt.all_tissues.egenes", 
                           data.table=F, header=F)
        colnames(total.expressed.genes) = "gene"
        colnames(lrt.egenes) = "gene"
        lrt.egenes = annotate.outliers(human, lrt.egenes)
        total.expressed.genes = annotate.outliers(human, total.expressed.genes)
        saveRDS(total.expressed.genes, "./data/total.annotated.genes.rds")
        total.expressed.genes = readRDS("./data/total.annotated.genes.rds")
        all.tissue.gene.list.enrichment = gene.list.enrichment("all.tissues", lrt.egenes$gene, total.expressed.genes)
      }
      
      # Fig 4B
      CairoPDF("../../materials/gene.set.enrichment.pdf")
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
}
