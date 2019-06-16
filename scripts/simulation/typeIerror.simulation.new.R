# calculate standard deviation using n
std <- function(x) {
  # calculate standard deviation with denominator of n
  return(sqrt(mean((x-mean(x, na.rm = T))^2, na.rm = T)))
}

# LRT
# quantitative trait
my_lrt_method <- function(E, G, idv.rare, e = exp(1), c = NULL, grid.search = T, amp.factor = 1) {
  # each row represents each individual, each column represents each variant

  # extract expression levels of idv with rare variants or without rare variants
  E.y = apply(idv.rare, 2, function(x) E[x])  # with rare variants
  E.x = apply(idv.rare, 2, function(x) E[!x]) # without rare variants

  K = ncol(G)
  N = nrow(G)
  n = colSums(idv.rare)
  m = colSums(!idv.rare)
  
  # set causal.rate
  if (is.null(c)) {
    c = runif(K)
  } else {
    c = rep(c, K)
  }
  
  # parameters
  mu.all = mean(E, na.rm = TRUE)
  sigma.all = std(E)
  
  mu.x = as.numeric(lapply(E.x, mean))
  sigma.x = as.numeric(lapply(E.x, std))
  sigma.x[sigma.x==0] <- sigma.all
  
  mu.y = as.numeric(lapply(E.y, mean))
  sigma.y = as.numeric(lapply(E.y, std))
  sigma.y[sigma.y==0] <- sigma.all
  
  # calculate log-transformed L_0
  l0 = sum(log(1 - c, base = e)) - (K * N / 2) * log(2 * pi * sigma.all ^ 2, base = e) - 
      K * N / 2
  
  # calculate log-transformed L0 + L1 
  # we use sum here because ln.a and ln.b are vectors
  ln.a = log(1 - c, base = e) - ((m + n) / 2) * log(2 * pi * sigma.all ^ 2, base = e) - N / 2

  if (grid.search) {
    exp.y.constant <- numeric(length(mu.y))
    exp.x.constant <- numeric(length(mu.x))
    for (i in 1:length(mu.y)) {
      exp.y.constant[i] <- sum(amp.factor ^ 2 * (E.y[[i]] - mu.y[i]) ^ 2)
      exp.x.constant[i] <- sum((E.x[[i]] - mu.x[i]) ^ 2)
    }
    
    space <- seq(-2, 2, 0.01)

    # changed
    sigma.x.range <- t(matrix(0, nrow=length(sigma.x), ncol = length(space)) + sigma.x) + space

    sigma.x.range <- t(sigma.x.range)
    sigma.x.range[sigma.x.range <= 0] = 0.01

    ln.bx.range <- matrix(0, ncol = length(space), nrow = length(c))
    #ln.b.max <- numeric(length(c))
    sigma.x.max <- numeric(length(c))
 
    sigma.y.range <- t(matrix(0, nrow=length(sigma.y), ncol = length(space)) + sigma.y) + space

    sigma.y.range <- t(sigma.y.range)
    sigma.y.range[sigma.y.range <= 0] = 0.01

    ln.by.range <- matrix(0, ncol = length(space), nrow = length(c))
    #ln.b.max <- numeric(length(c))
    sigma.y.max <- numeric(length(c))
    
    ln.bx.range <- log(c, base = e) - 
      (m / 2) * log(2 * pi * sigma.x.range ^ 2, base = e) - 
      1 / 2 * (exp.x.constant / sigma.x.range ^ 2) - 
      (n / 2) * log(2 * pi * sigma.y ^ 2, base = e) - 
      1 / 2 * (exp.y.constant / sigma.y ^ 2)
    
    ln.by.range <- log(c, base = e) - 
      (m / 2) * log(2 * pi * sigma.x ^ 2, base = e) - 
      1 / 2 * (exp.x.constant / sigma.x ^ 2) - 
      (n / 2) * log(2 * pi * sigma.y.range ^ 2, base = e) - 
      1 / 2 * (exp.y.constant / sigma.y.range ^ 2)
    
    max.y.idx <- apply(ln.by.range, 1, function(x) which(x == max(x))[1])
    max.x.idx <- apply(ln.bx.range, 1, function(x) which(x == max(x))[1])
    
    for (i in 1:length(c)) {
      #ln.b.max[i] <- ln.b.range[i, max.idx[i]]
      sigma.y.max[i] <- sigma.y.range[i, max.y.idx[i]]
      sigma.x.max[i] <- sigma.x.range[i, max.x.idx[i]]
    }
    #ln.b = ln.b.max
    ln.b <- log(c, base = e) - 
      (m / 2) * log(2 * pi * sigma.x.max ^ 2, base = e) - 
      1 / 2 * (exp.x.constant / sigma.x.max ^ 2) - 
      (n / 2) * log(2 * pi * sigma.y.max ^ 2, base = e) - 
      1 / 2 * (exp.y.constant / sigma.y.max ^ 2)
  } else {
    ln.b = log(c, base = e) - (m / 2) * log(2 * pi * sigma.x ^ 2, base = e) - 
      (n / 2) * log(2 * pi * sigma.y ^ 2, base = e) - N / 2 
  } 
      
  # changed
  l0l1 = sum(ln.a) + sum(log(1 + e ^ (ln.b - ln.a), base = e))
  l1 = l0l1 + log(1 - e ^ (l0 - l0l1))
  
  lr = l0 / l1
  
  return(lr)
}

LRT <- function(y, X, maf = 0.05, perm = 1000, e = exp(1), c = NULL, grid.search = T) {
  # get minor allele frequencies
  MAF = colMeans(X, na.rm=TRUE) / 2   
  # how many variants < maf
  rare.maf = MAF < maf & MAF > 0
  rare = sum(rare.maf)
  # expression with rare variants and genotypes of rare variants
  X = X[, rare.maf]
  idv.rare = X > 0
  
  # permutation
  # unpermuted LRT statistics
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
    perm.pval = (sum(x.perm >= lrt.stat) + 1) / (perm + 1)
  }
  
  # results
  name = "LRT: Likelihood Ratio Test"
  arg.spec = c(length(y), ncol(X), maf, perm)
  names(arg.spec) = c("samples", "variants", "maf", "n.perms")	
  res = list(lrt.stat = lrt.stat, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  return(res)
}

# binary trait
lrt1 <- function(n = 1000, m = 10, p.i = 0.05, p.plus = 0.50) {
  set.seed(0)
  c = runif(m)
  
  # x = log(1-c) + (2*n*p.i) * log(p.i) + (2*n - 2*n*p.i) * log(1-p.i)
  # k = log(c) + n * log(1-p.plus) + n * log(1-p.i) + (2*n*p.i) * (log(p.i) - log(1-p.i))
  # y = log(p.plus) - log(1-p.plus) + log(1-p.i) - log(p.i)
  # 
  # l0 = sum(x)
  # l0l1 = sum(x) + sum(log(1 + exp(n * p.plus * (k + y) - x)))
  # l1 = log(exp(l0l1) - exp(l0))
  
  x = (1-c) * p.i ^ (2 * n * p.i) * (1-p.i) ^ (2 * n - 2 * n * p.i)
  k = c * (1-p.plus) ^ n * (1-p.i) ^ n * (p.i / (1-p.i)) ^ (2 * n * p.i)
  y = (p.plus / (1-p.plus)) * ((1-p.i) / p.i)

  l0 = sum(log(x))
  l0l1 = sum(log(x + k * y ^ (n * p.plus)))
  l1 = log(exp(l0l1) - exp(l0))
  
  print(x)
  print(k)
  print(y)
  print(k * y ^ (n * p.plus))
  
  print(l0)
  print(l0l1)
  print(l1)
  print(l0-l1)
}

lrt2 <- function(n = 1000, m = 10, p.i = 0.05, p.plus = 0.50, e=exp(1)) {
  set.seed(0)
  c = runif(m)
  
  x = log(1-c, base = e) + (2*n*p.i) * log(p.i, base = e) + (2*n - 2*n*p.i) * log(1-p.i, base = e)
  k = log(c, base = e) + n * log(1-p.plus, base = e) + n * log(1-p.i, base = e) + (2*n*p.i) * (log(p.i, base = e) - log(1-p.i, base = e))
  y = log(p.plus, base = e) - log(1-p.plus, base = e) + log(1-p.i, base = e) - log(p.i, base = e)
  
  l0 = sum(x)
  l0l1 = sum(x) + sum(log(1 + e ^ (n * p.plus * y + k - x), base = e))
  # l0l1 = sum(log(x + exp(k) * exp(y) ^ (n * p.plus)))
  lr = l0l1 - l0
  
  # x = (1-c) * p.i ^ (2 * n * p.i) * (1-p.i) ^ (2 * n - 2 * n * p.i)
  # k = c * (1-p.plus) ^ n * (1-p.i) ^ n * (p.i / (1-p.i)) ^ (2 * n * p.i)
  # y = (p.plus / (1-p.plus)) * ((1-p.i) / p.i)
  
  # l0 = sum(log(x))
  # l0l1 = log(x + k * y ^ (n * p.plus))
  # l1 = log(exp(l0l1) - exp(l0))
  
  print(x)
  print(k)
  print(y)
  print(exp(k) * exp(y) ^ (n * p.plus))
  
  print(l0)
  print(l0l1)
  print(lr)
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
  mU <- colSums(X[unaffected, ], na.rm = T)
  nU <- colSums(!is.na(X)[unaffected, ], na.rm = T)
  q = (mU+1) / (2*nU+2)
  n = colSums(!is.na(X))
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
    X.agg = rowSums(X[, rare.maf], na.rm=T)
    # joining aggregated rare variants to common variants
    X.new = cbind(X[, !rare.maf], X.agg)
  }
  # linear regression and Wald statistic
  if (is.null(covariates)) {
    model <- lm(y ~ ., data = as.data.frame(X.new))
  } else {
    model <- lm(y ~ covariates + ., data = as.data.frame(X.new))
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
    X.agg = rowSums(X[,rare.maf], na.rm=TRUE)
    # joining aggregated rare variants to common variants
    X.new = cbind(X[, !rare.maf], X.agg)
  }
  # linear regression and Wald statistic
  if (is.null(covariates)) {
    model <- lm(y ~ ., data = as.data.frame(X.new))
  } else {
    model <- lm(y ~ covariates + ., data = as.data.frame(X.new))
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
my_cmc_method <- function(casecon, X.new) {
  # Internal function for CMAT method
  ## number of individuals N, cases nA, controls nU
  N = nrow(X.new)
  nA = sum(casecon)
  nU = N - nA
  ## matrix of genotypes in cases
  Xx = X.new[casecon==1,]  
  ## matrix of genotypes in controls  
  Yy = X.new[casecon==0,] 
  ## get means
  Xx.mean = colMeans(Xx, na.rm=TRUE)
  Yy.mean = colMeans(Yy, na.rm=TRUE)
  ## center matrices Xx and Yy
  Dx = sweep(Xx, 2, Xx.mean)
  Dy = sweep(Yy, 2, Yy.mean)
  
  ## pooled covariance matrix
  if (sum(complete.cases(X.new)) == N)  # no missing values
  {  
    COV = (t(Dx) %*% Dx + t(Dy) %*% Dy) / (N-2)
  } else {  # with missing values
    ## covariance matrix of cases
    tDx = t(Dx)
    Sx = matrix(0, ncol(X.new), ncol(X.new))
    for (i in 1:nrow(tDx))
    {
      for (j in i:ncol(Dx))
      {
        Sx[i,j] = sum(tDx[i,] * Dx[,j], na.rm=TRUE)
      }
    }
    sx.diag = diag(Sx)
    Sx = Sx + t(Sx)
    diag(Sx) = sx.diag
    ## covariance matrix of controls
    tDy = t(Dy)
    Sy = matrix(0, ncol(X.new), ncol(X.new))
    for (i in 1:nrow(tDy))
    {
      for (j in i:ncol(Dy))
      {
        Sy[i,j] = sum(tDy[i,] * Dy[,j], na.rm=TRUE)
      }
    }
    sy.diag = diag(Sy)
    Sy = Sy + t(Sy)
    diag(Sy) = sy.diag
    ## pooled covariance matrix
    COV = (1/(N-2)) * (Sx + Sy)	
  }
  
  ## general inverse
  if (nrow(COV) == 1) # only one variant
  { 
    if (COV < 1e-8) COV = 1e-8
    COV.inv = 1 / COV
  } else {
    COV.eigen = eigen(COV)
    eig.vals = COV.eigen$values  
    inv.vals = ifelse(abs(eig.vals) <= 1e-8, 0, 1/eig.vals)
    EV = solve(COV.eigen$vectors)
    COV.inv = t(EV) %*% diag(inv.vals) %*% EV
  }	
  
  ## Hotellings T2 statistic
  stat = t(Xx.mean - Yy.mean) %*% COV.inv %*% (Xx.mean - Yy.mean) * nA * nU / N
  as.numeric(stat)
}

CMC.casecon <- function(y, X, maf=0.05, perm=1000, threshold = 1.64) {
  # threshold is set based on normal distribution. assume top 5% expression levels are high (cases)
  case = y > threshold
  y[case] = 1
  y[!case] = 0
  ## number of individuals N
  N = nrow(X)
  ## get minor allele frequencies
  MAF = colMeans(X, na.rm=TRUE) / 2   
  ## how many variants < maf
  rare.maf = MAF < maf & MAF > 0
  rare = sum(rare.maf)
  ## collapsing
  if (rare <= 1) 
  {   
    # if rare variants <= 1, then NO collapse is needed
    X.new = X
  } else {
    # collapsing rare variants into one column
    X.collaps = rowSums(X[,rare.maf], na.rm=TRUE)
    X.collaps[X.collaps != 0] = 1
    # joining collapsed to common variants
    X.new = cbind(X[,!rare.maf], X.collaps)	   
  }
  ## change values to -1, 0, 1
  X.new = X.new - 1
  ## number of new variants
  M = ncol(X.new)
  ## Hotellings T2 statistic
  cmc.stat = my_cmc_method(y, X.new)
  
  ## Asymptotic p-values
  # under the null hypothesis T2 follows an F distribution 
  f.stat = cmc.stat * (N-M-1)/(M*(N-2))
  df1 = M          # degrees of freedom  
  df2 = N - M - 1  # degrees of freedom  
  asym.pval = 1 - pf(f.stat, df1, df2)
  
  ## under the alternative hyposthesis T2 follows a chi-square distr
  # pval = 1 - pchisq(cmc.stat, df=M)
  
  ## permutations
  perm.pval = NA
  if (perm > 0)
  {
    x.perm = rep(0, perm)
    for (i in 1:perm)
    {
      perm.sample = sample(1:length(y))
      x.perm[i] = my_cmc_method(y[perm.sample], X.new) 
    }
    # p-value 
    perm.pval = (sum(x.perm >= cmc.stat) + 1) / (perm + 1)
  }
  
  ## results
  name = "CMC: Combined Multivariate and Collapsing Method"
  arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare,  maf, perm)
  arg.spec = as.character(arg.spec)
  names(arg.spec) = c("cases", "controls", "variants", "rarevar", "maf", "perm")	
  res = list(cmc.stat = cmc.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  return(res)
}

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
    X.collaps = rowSums(X[,rare.maf], na.rm=TRUE)
    X.collaps[X.collaps != 0] = 1
    # joining collapsed to common variants
    X.new = cbind(X[, !rare.maf], X.collaps)
  }
  # linear regression and Wald statistic
  if (is.null(covariates)) {
    model <- lm(y ~ ., data = as.data.frame(X.new))
  } else {
    model <- lm(y ~ covariates + ., data = as.data.frame(X.new))
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
SKAT_ <- function(y, X, covariates = NULL, c) {
  require(SKAT)

  # changed
  weights <- rep(c, ncol(X))

  if (is.null(covariates)) {
    obj <- SKAT_Null_Model(y ~ 1, out_type = "C")
    skat.result <- SKAT(X, obj, method="SKATO", weights = weights)
    p <- skat.result$p.value
  } else {
    obj <- SKAT_Null_Model(y ~ covariates, out_type = "C")
    skat.result <- SKAT(X, obj, method="SKATO", weights = weights)
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

# rare variant type I error
#####
typeIerror.simulate <- function(cosi.hap, method = 6, perm = 10000, repeats = 1000, individual=1000, maf.cutoff = 0.05) {
  require(RNOmni)
  p = as.data.frame(matrix(nrow = repeats, ncol = method))
  colnames(p) <- c("CMC", "WSS-Agg", "BURDEN", "VT", "SKAT-O", "LRT-q")

  for (i in 1:repeats) {
    hap4k <- cosi.hap[sample(nrow(cosi.hap), 2*individual),]
    hap2k <- hap4k[1:individual, ] + hap4k[(individual+1):(2*individual), ]
    hap2k <- as.matrix(hap2k)
    maf.list <- colMeans(hap2k, na.rm = T) / 2
    hap2k[, maf.list > 0.5] <- 2 - hap2k[, maf.list > 0.5]
    hap2k <- hap2k[, maf.list != 0 & maf.list != 1]
    maf.list <- colMeans(hap2k, na.rm = T) / 2

    # only use rare variants
    hap2k <- hap2k[, maf.list < maf.cutoff]
    maf.list <- colMeans(hap2k, na.rm = T) / 2

    print(paste("rare variants", length(maf.list)))

    expr <- rnorm(individual) + rnorm(individual)
    expr <- rankNorm(expr)
    covariates <- NULL
    # changed

    p[i, 1] = CMC(expr, hap2k, covariates = covariates)$asym.pval
    p[i, 2] = WSS.agg(expr, hap2k, covariates = covariates)$asym.pval
    p[i, 3] = BURDEN(expr, hap2k, covariates = covariates)$asym.pval
    p[i, 4] = VT(expr, hap2k, perm = perm)$perm.pval
    p[i, 5] = SKAT_(expr, hap2k, covariates, c = 0.30)$p.value
    p[i, 6] = (lrt_perm(expr = expr, geno = hap2k, perm = perm, causal_ratio = rep(0.30, ncol(hap2k))) + 1) / (perm + 1)

    if (i %% (repeats / 10) == 0) {
      print(paste("round", i))
    }
  }  
  return(p)
}

typeIerror.check <- function(p.value) {
  fdr.table <- data.frame("alpha.level" = c(0.05, 0.01, 1e-3, 1e-4, 1e-5, 2.5e-6, 5e-6, 1e-6),
                          "fdr" = 0
  )
  for(i in 1:8) {
    fdr.table$fdr[i] <- sum(p.value < fdr.table$alpha.level[i]) / length(p.value)
  }
  return(fdr.table)
}

# main
args <- commandArgs(trailingOnly=T)
cosi.hap.file <- args[1]
perm <- as.numeric(args[2])
repeats <- as.numeric(args[3])
result.file <- args[4]

print(args)

require(Rcpp)
sourceCpp("/u/nobackup/eeskin2/k8688933/rare_var/script/lrt.cpp")

# hap is read from cosi result
require(data.table)
cosi.hap <- fread(cosi.hap.file)

typeI.p.values <- typeIerror.simulate(cosi.hap=cosi.hap, repeats=repeats, perm=perm)
write.csv(typeI.p.values, result.file, quote=F, row.names=F)
print(warnings())
