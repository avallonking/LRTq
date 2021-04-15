# calculate standard deviation using n
std <- function(x) {
  # calculate standard deviation with denominator of n
  return(sqrt(mean((x-mean(x, na.rm = T))^2, na.rm = T)))
}

# normalize CADD scores
simple.normalize <- function(x) {
  norm.val <- (x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
  if (all(is.na(norm.val))) {
      return(rep(0.5, length(x)))
  } else {
      return(norm.val)
  }
}

norm.l2 <- function(x) {
  res <- abs(x) / sqrt(sum(x^2, na.rm=T))
  #res[res == 1] = res[res == 1] - 0.01
  #res[res == 0] = res[res == 0] + 0.01
  return(res)
}

# LRT
# quantitative trait
my_lrt_method <- function(E, G, idv.rare, e = exp(1), c = NULL, grid.search = F) {
  require(plyr)
  # each row represents each individual, each column represents each variant

  # extract expression levels of idv with rare variants or without rare variants
  E.y = alply(idv.rare, 2, function(x) E[x])  # with rare variants
  E.x = alply(idv.rare, 2, function(x) E[!x]) # without rare variants

  K = ncol(G)
  N = nrow(G)
  n = colSums(idv.rare)
  m = colSums(!idv.rare)
  
  # set causal.rate
  if (is.null(c)) {
    c = runif(K)
  } else if (length(c) == 1) {
    c = rep(c, K)
  } else if (length(c) > 1) {
    c = c
  }
  
  # parameters
  exp.x.constant <- unlist(lapply(E.x, function(x) {sum( (x-mean(x))^2 )}))
  exp.y.constant <- unlist(lapply(E.y, function(y) {sum( (y-mean(y))^2 )}))
  sigma.all = std(E)
  sigma.var <- sqrt( ( exp.x.constant + exp.y.constant ) / (m+n) )
  
  
  if (grid.search) {
    space <- seq(-1, 1, 0.01)
    sigma.var.range <- t(matrix(0, nrow=length(sigma.var), ncol = length(space)) + sigma.var) + space
    sigma.var.range <- t(sigma.var.range)
    
    ln.b.range <- log(c, base = e) - ((m+n) / 2) * log(2 * pi * sigma.var.range ^ 2, base = e) -
        1 / (2 * sigma.var.range ^ 2) * (exp.x.constant + exp.y.constant)
    
    ln.b <- apply(ln.b.range, 1, max)
  } else {
    ln.b = log(c, base = e) - ((m+n) / 2) * log(2 * pi * sigma.var ^ 2, base = e) - N / 2
  }
      
  ln.a = log(1 - c, base = e) - ((m + n) / 2) * log(2 * pi * sigma.all ^ 2, base = e) - N / 2
  l0 = sum(ln.a)
  l0l1 = sum(ln.a) + sum(log(1 + e ^ (ln.b - ln.a), base = e))
  l1 = l0l1 + log(1 - e ^ (l0 - l0l1))
  lr = l1 - l0

  return(lr)
}

LRT <- function(y, X, perm = 1000, e = exp(1), c = NULL, maf = 0.05, grid.search = F) {
  # get minor allele frequencies
  MAF = colMeans(X, na.rm=TRUE) / 2   
  # how many variants < maf
  rare.maf = MAF < maf & MAF > 0
  rare = sum(rare.maf)
  # expression with rare variants and genotypes of rare variants
  X = as.matrix(X[, rare.maf])
  idv.rare = as.matrix(X > 0)
  
  # permutation
  # unpermuted LRT statistics
  # changed
  if (is.null(c)) {
    c = 1 - MAF[rare.maf]
    #print(paste("LRT causal ratio:", c))
  } else if (length(c) > 1) {
    c = c[rare.maf]
  }
  lrt.stat = my_lrt_method(y, X, idv.rare, e, c, grid.search = grid.search)
  
  # For a specific score, returns how often permuted data has a higher score than unpermuted data.
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_lrt_method(y[perm.sample], X, idv.rare, e, c, grid.search = grid.search)
    }
    # p-value
    sig.perm = sum(x.perm >= lrt.stat)
    perm.pval = (sig.perm + 1) / (perm + 1)
  }
  
  # results
  name = "LRT: Likelihood Ratio Test"
  arg.spec = c(length(y), ncol(X), perm)
  names(arg.spec) = c("samples", "variants", "n.perms")	
  res = list(lrt.stat = lrt.stat, 
             perm.pval = perm.pval,
             sig.perm = sig.perm,
             args = arg.spec, 
             name = name)
  return(res)
}

# VT
compute.vt <- function(y, z.num.pre, z.denom) {
  y.new = y - mean(y)
  z.num = unlist(lapply(z.num.pre, function(x) {sum(x * y.new, na.rm = T)}))
  z.scores = z.num / z.denom
  vt.stat = max(z.scores, na.rm = T)
  return(vt.stat)
}

VT <- function(y, X, perm=1000, ksi=NULL) {
  ## running vt method
  mafs = (1 + colSums(X, na.rm=TRUE)) / (2 + 2*nrow(X))
  h.maf = sort(unique(mafs))
  z.num.pre = vector("list", length(h.maf)-1)
  z.denom = rep(0, length(h.maf)-1)

  if (length(h.maf) == 1) {
    name = "VT: Variable Threshold"
    arg.spec = c(length(y), ncol(X), perm)
    names(arg.spec) = c("samples", "variants", "n.perms")
    res = list(vt.stat = NA,
               perm.pval = NA,
               sig.perm = NA,
               args = arg.spec,
               name = name)
    return(res)
  }

  # precompute unchanged variables
  for (i in 1:(length(h.maf)-1))
  {
    if (is.null(ksi)) {
      selected.var = mafs < h.maf[i+1]
      z.num.pre[[i]] = X[, selected.var]
      z.denom[i] = sqrt(sum((X[, selected.var]) ^ 2, na.rm=TRUE))
    } else {
      selected.var <- mafs < h.maf[i+1]
      z.num.pre[[i]] = t(ksi * t(X))[, selected.var]
      z.denom[i] = sqrt(sum((t(ksi * t(X)))[, selected.var] ^ 2, na.rm=TRUE))
    }
  }
  
  ## Unpermuted VT statistics
  vt.stat = compute.vt(y, z.num.pre, z.denom)
  ## permutations
  ## For a specific score, returns how often permuted data has a higher score than unpermuted data.
  perm.pval = NA
  sig.perm = 0
  if (perm > 0)
  {
    y.perm = replicate(perm, sample(y))
    x.perm = apply(y.perm, 2, function(y) compute.vt(y, z.num.pre, z.denom))
    sig.perm = sum(x.perm >= vt.stat)
  }
  ## p-value
  perm.pval = (sig.perm + 1) / (perm + 1)
  
  ## results
  name = "VT: Variable Threshold"
  arg.spec = c(length(y), ncol(X), perm)
  names(arg.spec) = c("samples", "variants", "n.perms")	
  res = list(vt.stat = vt.stat, 
             perm.pval = perm.pval, 
             sig.perm = sig.perm,
             args = arg.spec, 
             name = name)
  return(res)
}

# Weighted methods - WSS
# use linear regression + Wald statistic
WSS.agg <- function(y, X, threshold = 1.64, covariates = NULL, maf = 0.05) {
  require(car)
  # each row represents each individual, each column represents each variant
  # assume that all variants are in one small-size gene, so we do need to 
  # separate rare variants into different groups
  # get minor allele frequencies
  MAF = colMeans(X, na.rm=TRUE) / 2   
  # how many variants < maf
  rare.maf = MAF < maf & MAF > 0
  rare = sum(rare.maf)
  
  # threshold is set based on normal distribution. define top 5% expression as high
  unaffected <- y <= threshold
  # calculate WSS weight
  mU <- colSums(as.matrix(X[unaffected, ]), na.rm = T)
  nU <- colSums(as.matrix(!is.na(X)[unaffected, ]), na.rm = T)
  q = (mU+1) / (2*nU+2)
  n = colSums(as.matrix(!is.na(X)))
  w = sqrt(n * q * (1 - q))
  # WSS aggregation
  # multiply the weights
  X <- t(t(X) * w)
  if (rare <= 1) 
  {   
    # if rare variants <= 1, then NO aggregation is needed
    X.new = X
  } else {
    # aggregate rare variants into one column
    X.agg = rowSums(as.matrix(X[, rare.maf]), na.rm=T)
    # joining aggregated rare variants to common variants
    X.new = cbind(X[, !rare.maf], X.agg)
  }
  # linear regression and Wald statistic
  if (is.null(covariates)) {
    model <- lm(y ~ ., data = as.data.frame(X.new))
    if (sum(is.na(model$coefficients)) > 0) {
        model <- lm(y ~ ., data = as.data.frame(X.new[, !is.na(model$coefficients)[-1] ]))
    }
  } else {
    model <- lm(y ~ covariates + ., data = as.data.frame(X.new))
    if (sum(is.na(model$coefficients)) > 0) {
        model <- lm(y ~ covariates + ., data = as.data.frame(X.new[, !is.na(model$coefficients)[-1] ]))
    }
  }
  wald <- Anova(model, type = "II", test.statistic = "Wald")
  ## results
  name = "Weighted Aggregation: WSS Aggregation Method (Linear regression + Wald statistic)"
  arg.spec = c(length(y), ncol(X), rare,  maf)
  arg.spec = as.character(arg.spec)
  names(arg.spec) = c("samples", "variants", "rarevar", "maf")	
  res = list(wss.stat = rev(wald$`F value`)[2], 
             asym.pval = rev(wald$`Pr(>F)`)[2], 
             args = arg.spec, 
             name = name)
  return(res)
}

# Can treat low expression level as control, while high expression level as case
my_wss_method <- function(casecon, gen) {	
  controls = !casecon
  n.loc = ncol(gen)
  
  ## calculate weights
  w <- rep(0, n.loc)
  for (j in 1:n.loc)
  {
    nNA = is.na(gen[,j])
    mU = sum(gen[controls,j], na.rm=TRUE)
    nU = sum(controls[!nNA])		
    q = (mU+1) / (2*nU+2)
    n = sum(!nNA)
    w[j] = sqrt(n * q * (1-q))
  }
  
  ## calculate genetic score	
  score = rowSums(gen %*% diag(1/w), na.rm=TRUE) 
  # rank.score = order(score)
  rank.score = rank(score)
  
  ## sum of ranks of cases
  x = sum(rank.score[casecon==1])
  return(x)
}

WSS <- function(y, X, perm=1000, threshold = 1.64) {
  # threshold is set based on normal distribution. assume top 5% expression levels are high (cases)
  case = y > threshold
  y[case] = 1
  y[!case] = 0
  
  ## running wss method
  wss.stat = my_wss_method(y, X)  
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_wss_method(y[perm.sample], X) 
    }
    # p-value 
    perm.pval = (sum(x.perm >= wss.stat) + 1) / (perm + 1)
  }
  
  ## results
  name = "WSS: Weighted Sum Statistic"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm)
  names(arg.spec) = c("cases", "controls", "variants", "n.perms")	
  res = list(wss.stat = wss.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  return(res)
}

# Aggregation methods - BURDEN
BURDEN <- function(y, X, maf=0.05, covariates = NULL) {
  require(car)
  # each row represents each individual, each column represents each variant
  # assume that all variants are in one small-size gene, so we do need to 
  # separate rare variants into different groups
  
  # get minor allele frequencies
  MAF = colMeans(X, na.rm=TRUE) / 2   
  # how many variants < maf
  rare.maf = MAF < maf & MAF > 0
  rare = sum(rare.maf)
  # aggregation
  if (rare <= 1) 
  {   
    # if rare variants <= 1, then NO aggregation is needed
    X.new = X
  } else {
    # aggregate rare variants into one column
    X.agg = rowSums(as.matrix(X[,rare.maf]), na.rm=TRUE)
    # joining aggregated rare variants to common variants
    X.new = cbind(X[, !rare.maf], X.agg)
  }
  # linear regression and Wald statistic
  if (is.null(covariates)) {
    model <- lm(y ~ ., data = as.data.frame(X.new))
    if (sum(is.na(model$coefficients)) > 0) {
        model <- lm(y ~ ., data = as.data.frame(X.new[, !is.na(model$coefficients)[-1] ]))
    }
  } else {
    model <- lm(y ~ covariates + ., data = as.data.frame(X.new))
    if (sum(is.na(model$coefficients)) > 0) {
        model <- lm(y ~ covariates + ., data = as.data.frame(X.new[, !is.na(model$coefficients)[-1] ]))
    }
  }
  wald <- Anova(model, type = "II", test.statistic = "Wald")
  ## results
  name = "BURDEN: Aggregation Method (Linear regression + Wald statistic)"
  arg.spec = c(length(y), ncol(X), rare,  maf)
  arg.spec = as.character(arg.spec)
  names(arg.spec) = c("samples", "variants", "rarevar", "maf")	
  res = list(burden.stat = rev(wald$`F value`)[2], 
             asym.pval = rev(wald$`Pr(>F)`)[2], 
             args = arg.spec, 
             name = name)
  return(res)
}

# Collapsing methods - CMC
CMC <- function(y, X, maf=0.05, covariates = NULL) {
  require(car)
  # each row represents each individual, each column represents each variant
  # assume that all variants are in one small-size gene, so we do need to 
  # separate rare variants into different groups
  
  # get minor allele frequencies
  MAF = colMeans(X, na.rm=TRUE) / 2   
  # how many variants < maf
  rare.maf = MAF < maf & MAF > 0
  rare = sum(rare.maf)
  # collapsing
  if (rare <= 1) 
  {   
    # if rare variants <= 1, then NO collapse is needed
    X.new = X
  } else {
    # collapsing rare variants into one column
    X.collaps = rowSums(as.matrix(X[,rare.maf]), na.rm=TRUE)
    X.collaps[X.collaps != 0] = 1
    if (all(X.collaps == 1)) {
        ## results
        name = "CMC: Combined Multivariate and Collapsing Method (Linear regression + Wald statistic)"
        arg.spec = c(length(y), ncol(X), rare,  maf)
        arg.spec = as.character(arg.spec)
        names(arg.spec) = c("samples", "variants", "rarevar", "maf")	
        res = list(cmc.stat = NA, 
                   asym.pval = NA, 
                   args = arg.spec, 
                   name = name)
        return(res)
    }
    # joining collapsed to common variants
    X.new = cbind(X[, !rare.maf], X.collaps)
  }

  # linear regression and Wald statistic
  if (is.null(covariates)) {
    model <- lm(y ~ ., data = as.data.frame(X.new))
    if (sum(is.na(model$coefficients)) > 0) {
        model <- lm(y ~ ., data = as.data.frame(X.new[, !is.na(model$coefficients)[-1] ]))
    }
  } else {
    model <- lm(y ~ covariates + ., data = as.data.frame(X.new))
    if (sum(is.na(model$coefficients)) > 0) {
        model <- lm(y ~ covariates + ., data = as.data.frame(X.new[, !is.na(model$coefficients)[-1] ]))
    }
  }

  wald <- Anova(model, type = "II", test.statistic = "Wald")
  ## results
  name = "CMC: Combined Multivariate and Collapsing Method (Linear regression + Wald statistic)"
  arg.spec = c(length(y), ncol(X), rare,  maf)
  arg.spec = as.character(arg.spec)
  names(arg.spec) = c("samples", "variants", "rarevar", "maf")	
  res = list(cmc.stat = rev(wald$`F value`)[2], 
             asym.pval = rev(wald$`Pr(>F)`)[2], 
             args = arg.spec, 
             name = name)
  return(res)
}

# SKAT
SKAT_ <- function(y, X, covariates = NULL, weight = NULL) {
  require(SKAT)
  if (is.null(covariates)) {
    obj <- SKAT_Null_Model(y ~ 1, out_type = "C")
    skat.result <- SKAT(X, obj, method="SKATO", weights = weight)
    p <- skat.result$p.value
  } else {
    obj <- SKAT_Null_Model(y ~ covariates, out_type = "C")
    skat.result <- SKAT(X, obj, method="SKATO", weights = weight)
    p <- skat.result$p.value
  }
  ## results
  name = "SKAT-O: SNP-set (Sequence) Kernel Association Test"
  arg.spec = c(length(y), ncol(X))
  names(arg.spec) = c("samples", "variants")	
  res = list(p.value = p, 
             skat.obj = obj,
             skat.result = skat.result,
             args = arg.spec, 
             name = name)
  return(res)
}

##### apply to GTEx
remove.missing.impute <- function(geno.matrix) {
  # remove missing rate >= 0.05
  if (dim(geno.matrix)[1] == 0) {
      return(NULL)
  }
  good.var.idx <- rowSums(geno.matrix == -1) < 0.05 * (dim(geno.matrix)[2] - 1)
  if (sum(good.var.idx) == 0) {
      return(NULL)
  } else {
      geno.matrix <- geno.matrix[good.var.idx, ]
      geno.matrix[geno.matrix == -1] <- 0
      return(geno.matrix)
  }
}

get.weights <- function(weights, selected.snp, rare.var.idx, norm="none", impute="median") {
     weights.idx <- match(selected.snp, weights$V1)
     selected.weights <- weights$V2[weights.idx]

     if (all(is.na(selected.weights))) {
         selected.weights[is.na(selected.weights)] <- median(weights$V2, na.rm=T)
     } else if (impute == "median") {
         selected.weights[is.na(selected.weights)] <- median(selected.weights, na.rm=T)
     } else if (impute == "zero") {
         selected.weights[is.na(selected.weights)] <- 0
     }

     if (norm == "l2") {
         selected.weights <- norm.l2(selected.weights)
     } else if (norm == "20k") {
         selected.weights <- 1 - selected.weights / 20000
     }
     selected.weights[selected.weights == 1] <- selected.weights[selected.weights == 1] - 0.01
     selected.weights[selected.weights == 0] <- selected.weights[selected.weights == 0] + 0.01
     selected.weights <- selected.weights[rare.var.idx]
     return(selected.weights)
}

prepare.vt.ind <- function(E) {
    data.frame("indid" = c(0:(length(E) - 1)), "pheno" = E)
}

prepare.vt.csnp <- function(weights) {
    data.frame("snpid" = c(0:(length(weights) - 1)), "polyphen" = weights)
}

prepare.vt.cgeno <- function(G) {
    a = as.data.frame(which(G == 1, arr.ind = T) - 1)
    b = as.data.frame(which(G == 2, arr.ind = T) - 1)
    colnames(a) = c("indid", "snpid")
    colnames(b) = c("indid", "snpid")
    if (nrow(a) != 0 & nrow(b) != 0) {
        a$count = 1
        b$count = 2
        return(rbind(a, b))
    } else if (nrow(a) == 0) {
        b$count = 2
        return(b)
    } else {
        a$count = 1
        return(a)
    }
}

perform.test <- function(tissue.expr, geno.matrix, gene.set, gene.list, covar, start, end, weights.file.prefix, common.eqtl, common.geno.matrix, maf.cutoff = 0.05, perm = 100000, nodes=16) {
    #nodes = detectCores()
    require(stringr)
    require(RNOmni)
    require(data.table)
    require(doMC)
    p.values.vt <- as.data.frame(matrix(nrow = end-start+1, ncol = 9))
    p.values.lrt <- as.data.frame(matrix(nrow = end-start+1, ncol = 9))
    p.values.skat <- as.data.frame(matrix(nrow = end-start+1, ncol = 9))
    colnames(p.values.vt) <- c("Gene_ID", "Distance", "MAFxCADD", "MAF", "MAFxDistance", "DistancexCADDxMAF", "Unweighted", "LINSIGHT", "CADD")
    colnames(p.values.lrt) <- c("Gene_ID", "Distance", "MAFxCADD", "MAF", "MAFxDistance", "DistancexCADDxMAF", "Unweighted", "LINSIGHT", "CADD")
    colnames(p.values.skat) <- c("Gene_ID", "Distance", "MAFxCADD", "MAF", "MAFxDistance", "DistancexCADDxMAF", "Unweighted", "LINSIGHT", "CADD")
    selected.idv.idx <- which(colnames(tissue.expr) %in% colnames(geno.matrix))
    selected.idv.idx1 <- which(colnames(geno.matrix) %in% colnames(tissue.expr)[selected.idv.idx])
    geno.matrix <- geno.matrix[, c(1, selected.idv.idx1)]
    
    selected.idv.idx2 <- which(covar$ID %in% colnames(geno.matrix))
    covar <- covar[selected.idv.idx2, -1]
    
    selected.idv.idx3 <- which(colnames(common.geno.matrix) %in% colnames(geno.matrix))
    common.geno.matrix <- common.geno.matrix[, selected.idv.idx3]

    weights.files <- paste(weights.file.prefix, c("maf", "cadd_scaled", "dis", "linsight"), sep=".")
    weights <- lapply(weights.files, fread)

    for (i in start:end) {
         if (i > nrow(gene.list)) {
             break
         }
         selected.gene <- gene.list$V1[i]

         if (!selected.gene %in% common.eqtl$gene_id) {
             i <- i-start+1
             p.values.vt[i, 1] = selected.gene
             p.values.lrt[i, 1] = selected.gene
             p.values.skat[i, 1] = selected.gene
             print("VT")
             print(p.values.vt[i, ])
             print("LRT")
             print(p.values.lrt[i, ])
             print("SKAT")
             print(p.values.skat[i, ])
             next
         }

         common.snp <- common.eqtl$variant_id[common.eqtl$gene_id == selected.gene]
         common.snp.genotype <- common.geno.matrix[common.geno.matrix$ID == common.snp, -1]
         common.snp.genotype <- remove.missing.impute(common.snp.genotype)

         included.snp <- gene.set$V2[gene.set$V1 == selected.gene]

         G <- geno.matrix[geno.matrix$ID %in% included.snp, ]
         G <- remove.missing.impute(G)

         if (is.null(G) | is.null(common.snp.genotype) | selected.gene == "ENSG00000183117") {
            i <- i-start+1
            p.values.vt[i, 1] = selected.gene
            p.values.lrt[i, 1] = selected.gene
            p.values.skat[i, 1] = selected.gene
            print("VT")
            print(p.values.vt[i, ])
            print("LRT")
            print(p.values.lrt[i, ])
            print("SKAT")
            print(p.values.skat[i, ])
            next
         }

         new.covar <- cbind(covar, as.numeric(common.snp.genotype))

         selected.snp <- G[, 1]
         G <- G[, -1]
         G <- t(G)
         maf.list <- colMeans(G, na.rm = T) / 2
         G[, maf.list >= 0.50] = 2 - G[, maf.list >= 0.50]
         maf.list <- colMeans(G, na.rm = T) / 2

         rare.var.idx <- as.logical(maf.list > 0 & maf.list < maf.cutoff)

         # changed
         print(selected.gene)
         print(i)
         print(paste("rare variants", sum(rare.var.idx)))

         if (sum(rare.var.idx) < 1) {
             i <- i-start+1
             p.values.vt[i, 1] = selected.gene
             p.values.lrt[i, 1] = selected.gene
             p.values.skat[i, 1] = selected.gene
             print("VT")
             print(p.values.vt[i, ])
             print("LRT")
             print(p.values.lrt[i, ])
             print("SKAT")
             print(p.values.skat[i, ])
             next
         }

         G <- as.matrix(G[, rare.var.idx])

         expr <- as.numeric(tissue.expr[tissue.expr$gene_id==selected.gene, selected.idv.idx])
         E <- rankNorm(lm(expr ~ as.matrix(new.covar))$residuals)

         selected.weights.unweighted <- rep(0.30, ncol(G))
         selected.weights.maf <- get.weights(weights[[1]], selected.snp, rare.var.idx)
         selected.weights.cadd <- get.weights(weights[[2]], selected.snp, rare.var.idx, norm="l2")
         selected.weights.dis <- get.weights(weights[[3]], selected.snp, rare.var.idx, norm="20k", impute="zero")
         selected.weights.linsight <- get.weights(weights[[4]], selected.snp, rare.var.idx)

         i <- i-start+1

         p.values.vt[i, 1] = selected.gene

         ind.vt = prepare.vt.ind(E)
         cgeno.vt = prepare.vt.cgeno(G)

         p.values.vt[i, 2] = vt.fast(ind.vt, prepare.vt.csnp(selected.weights.dis), cgeno.vt, permutations = perm, nodes = nodes)
         p.values.vt[i, 3] = vt.fast(ind.vt, prepare.vt.csnp(selected.weights.maf * selected.weights.cadd), cgeno.vt, permutations = perm, nodes = nodes)
         p.values.vt[i, 4] = vt.fast(ind.vt, prepare.vt.csnp(selected.weights.maf), cgeno.vt, permutations = perm, nodes = nodes)
         p.values.vt[i, 5] = vt.fast(ind.vt, prepare.vt.csnp(selected.weights.maf * selected.weights.dis), cgeno.vt, permutations = perm, nodes = nodes)
         p.values.vt[i, 6] = vt.fast(ind.vt, prepare.vt.csnp(selected.weights.dis * selected.weights.cadd * selected.weights.maf), cgeno.vt, permutations = perm, nodes = nodes)
         p.values.vt[i, 7] = vt.fast(ind.vt, prepare.vt.csnp(selected.weights.unweighted), cgeno.vt, permutations = perm, nodes = nodes)
         p.values.vt[i, 8] = vt.fast(ind.vt, prepare.vt.csnp(selected.weights.linsight), cgeno.vt, permutations = perm, nodes = nodes)
         p.values.vt[i, 9] = vt.fast(ind.vt, prepare.vt.csnp(selected.weights.cadd), cgeno.vt, permutations = perm, nodes = nodes)

         p.values.lrt[i, 1] = selected.gene
         registerDoMC()
         p.values.lrt[i, 2] = ((foreach(i = 1:nodes, .combine = '+') %dopar% lrt_perm(expr = E, geno = G, causal_ratio = selected.weights.dis, perm = perm / nodes)) + 1) / (perm + 1)
         p.values.lrt[i, 3] = ((foreach(i = 1:nodes, .combine = '+') %dopar% lrt_perm(expr = E, geno = G, causal_ratio = selected.weights.maf * selected.weights.cadd, perm = perm / nodes)) + 1) / (perm + 1)
         p.values.lrt[i, 4] = ((foreach(i = 1:nodes, .combine = '+') %dopar% lrt_perm(expr = E, geno = G, causal_ratio = selected.weights.maf, perm = perm / nodes)) + 1) / (perm + 1)
         p.values.lrt[i, 5] = ((foreach(i = 1:nodes, .combine = '+') %dopar% lrt_perm(expr = E, geno = G, causal_ratio = selected.weights.maf * selected.weights.dis, perm = perm / nodes)) + 1) / (perm + 1)
         p.values.lrt[i, 6] = ((foreach(i = 1:nodes, .combine = '+') %dopar% lrt_perm(expr = E, geno = G, causal_ratio = selected.weights.dis * selected.weights.cadd * selected.weights.maf, perm = perm / nodes)) + 1) / (perm + 1)
         p.values.lrt[i, 7] = ((foreach(i = 1:nodes, .combine = '+') %dopar% lrt_perm(expr = E, geno = G, causal_ratio = selected.weights.unweighted, perm = perm / nodes)) + 1) / (perm + 1)
         p.values.lrt[i, 8] = ((foreach(i = 1:nodes, .combine = '+') %dopar% lrt_perm(expr = E, geno = G, causal_ratio = selected.weights.linsight, perm = perm / nodes)) + 1) / (perm + 1)
         p.values.lrt[i, 9] = ((foreach(i = 1:nodes, .combine = '+') %dopar% lrt_perm(expr = E, geno = G, causal_ratio = selected.weights.cadd, perm = perm / nodes)) + 1) / (perm + 1)
         
         p.values.skat[i, 1] = selected.gene
         p.values.skat[i, 2] = SKAT_(E, G, covariates = NULL, weight = selected.weights.dis)$p.value
         p.values.skat[i, 3] = SKAT_(E, G, covariates = NULL, weight = selected.weights.maf * selected.weights.cadd)$p.value
         p.values.skat[i, 4] = SKAT_(E, G, covariates = NULL, weight = selected.weights.maf)$p.value
         p.values.skat[i, 5] = SKAT_(E, G, covariates = NULL, weight = selected.weights.maf * selected.weights.dis)$p.value
         p.values.skat[i, 6] = SKAT_(E, G, covariates = NULL, weight = selected.weights.dis * selected.weights.cadd * selected.weights.maf)$p.value
         p.values.skat[i, 7] = SKAT_(E, G, covariates = NULL, weight = selected.weights.unweighted)$p.value
         p.values.skat[i, 8] = SKAT_(E, G, covariates = NULL, weight = selected.weights.linsight)$p.value
         p.values.skat[i, 9] = SKAT_(E, G, covariates = NULL, weight = selected.weights.cadd)$p.value
         
         print("VT")
         print(p.values.vt[i, ])
         print("LRT")
         print(p.values.lrt[i, ])
         print("SKAT")
         print(p.values.skat[i, ])
    }
    return(list(p.values.vt, p.values.lrt, p.values.skat))
    #return(list(p.values.vt, p.values.lrt))
}

# main
args <- commandArgs(trailingOnly=T)
tissue.expr.file <- args[1]
geno.matrix.file <- args[2]
gene.set.file <- args[3]
covariates.file <- args[4]
gene.list.file <- args[5]
start <- as.numeric(args[6])
end <- as.numeric(args[7])
result.file.prefix <- args[8]
weights.file.prefix <- args[9]
common.eqtl.file <- args[10]
common.geno.matrix.file <- args[11]

print(args)

require(Rcpp)
require(parallel)
sourceCpp("./methods_source/lrt.cpp")
source("./methods_source/rareVariantTests.functions_only.R")

require(data.table)
geno.matrix <- fread(geno.matrix.file, data.table=F)
tissue.expr <- fread(cmd=paste("zcat ", tissue.expr.file), data.table=F)
# tissue.expr$gene_id <- gsub("\\..*", "", tissue.expr$gene_id)
gene.set <- fread(gene.set.file, header=F, data.table=F)
gene.list <- fread(gene.list.file, header=F, data.table=F)
covar <- fread(covariates.file, data.table=F)
result.file <- paste(result.file.prefix, c("vt.csv", "lrt.csv", "skat.csv"), sep=".")
common.eqtl <- fread(common.eqtl.file, data.table=F)
common.geno.matrix <- fread(common.geno.matrix.file, data.table=F)

nodes <- parallel::detectCores(all.tests=T)
nodes <- 4
p.values <- perform.test(tissue.expr=tissue.expr, geno.matrix=geno.matrix, gene.set=gene.set, covar=covar, gene.list=gene.list, start=start, end=end, weights.file.prefix=weights.file.prefix, nodes=nodes, common.eqtl=common.eqtl, common.geno.matrix=common.geno.matrix)
write.csv(p.values[[1]], result.file[1], quote=F, row.names=F)
write.csv(p.values[[2]], result.file[2], quote=F, row.names=F)
write.csv(p.values[[3]], result.file[3], quote=F, row.names=F)
print(warnings())
