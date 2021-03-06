#include <Rcpp.h>

using namespace Rcpp;

double ss(NumericVector x) {
  return mean(pow(x-mean(x), 2));
}

NumericVector find_sigma_var(NumericVector expr, LogicalMatrix idv_rare, int K, int N) {
  NumericVector expr_x, expr_y, sigma_var(K);
  int i;
  
  for (i = 0; i < K; i++) {
    expr_x = expr[!(idv_rare(_, i))];
    expr_y = expr[(idv_rare(_, i))];
    sigma_var(i) = (sum(pow(expr_x - mean(expr_x), 2)) + sum(pow(expr_y - mean(expr_y), 2))) / N;
  }
  return sigma_var;
}

double compute_lrt(double l0, NumericVector ln_a, NumericVector ln_b_pre, NumericVector sigma_var, int K, int N) {
  NumericVector ln_b(K);
  double lr, l0l1, l1;
  ln_b = ln_b_pre - N / 2 * log(sigma_var);
  l0l1 = l0 + sum(log(1 + exp(ln_b - ln_a)));
  l1 = l0l1 + log(1 - exp(l0 - l0l1));
  lr = l1 - l0;
  return lr;
}

// [[Rcpp::export]]
double lrt_perm(NumericVector expr, IntegerMatrix geno, NumericVector causal_ratio, int perm = 1000) {
  // all variants should be rare here
  int K, N, i, t, significant_lrt_stat = 0;
  LogicalMatrix idv_rare(geno.nrow(), geno.ncol());
  NumericVector perm_sample;
  double lrt_stat, perm_pval, perm_stat, sigma, l0;
  
  K = geno.ncol();
  N = geno.nrow();
  
  for (i = 0; i < K; i++) {
    idv_rare(_, i) = geno(_, i) > 0;
  }
  
  // precompute unpermuted variables
  NumericVector ln_a(K), ln_b_pre(K), sigma_var(K);
  sigma = ss(expr);
  sigma_var = find_sigma_var(expr, idv_rare, K, N);
  ln_a = log(1 - causal_ratio) - N / 2 * log(2 * PI * sigma) - N / 2;
  ln_b_pre = log(causal_ratio) - N / 2 * log(2 * PI) - N / 2;
  l0 = sum(ln_a);

  lrt_stat = compute_lrt(l0, ln_a, ln_b_pre, sigma_var, K, N);
  
  for (t = 0; t < perm; t++) {
    perm_sample = sample(expr.size(), expr.size()) - 1;
    sigma_var = find_sigma_var(expr[perm_sample], idv_rare, K, N);
    perm_stat = compute_lrt(l0, ln_a, ln_b_pre, sigma_var, K, N);
    if (perm_stat >= lrt_stat)
      significant_lrt_stat += 1;
  }
  // std::cout << significant_lrt_stat;
  perm_pval = static_cast<double>(significant_lrt_stat + 1) / (perm + 1);
  
  // return perm_pval;
  return significant_lrt_stat;
}

// [[Rcpp::export]]
double vt_method(NumericVector expr, NumericMatrix geno, NumericVector mafs, NumericVector h_maf) {
  NumericVector z_scores(h_maf.size() - 1), expr_new = expr - mean(expr), temp;
  int i, t;
  double vt_stat, z_num, z_denom;
  LogicalVector selected;
  
  for (i = 0; i < h_maf.size() - 1; i++) {
    z_num = z_denom = 0;
    selected = mafs < h_maf(i + 1);
    for (t = 0; t < geno.nrow(); t++) {
      temp = geno(t, _);
      temp = temp[selected];
      z_num += sum(temp * expr_new(t));
      z_denom += sum(pow(temp, 2));
    }
    z_scores(i) = z_num / sqrt(z_denom);
  }
  vt_stat = max(z_scores);
  return vt_stat;
}

// [[Rcpp::export]]
double vt_perm(NumericVector expr, NumericMatrix geno, int perm = 1000) {
  NumericVector mafs, h_maf, perm_sample, temp;
  double perm_pval, vt_stat, perm_stat;
  int i, significant_vt_stat = 0;
  
  temp = (1 + colSums(geno));
  mafs = temp / (2 + 2 * geno.nrow());
  h_maf = sort_unique(mafs);
  vt_stat = vt_method(expr, geno, mafs, h_maf);
  
  for (i = 0; i < perm; i++) {
    perm_sample = sample(expr.size(), expr.size()) - 1;
    perm_stat = vt_method(expr[perm_sample], geno, mafs, h_maf);
    if (perm_stat >= vt_stat)
      significant_vt_stat += 1;
  }
  perm_pval = static_cast<double>(significant_vt_stat + 1) / (perm + 1);

  return significant_vt_stat;
}

/*** R
*/
