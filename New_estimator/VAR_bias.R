source("VAR_func.R")
source("Cov_func.R")
set.seed(101)

var_coverage <- function(N = 1e3, phi, omega, B = 1000, level = .90)
{
  
  p <- dim(phi)[1]
  
  large.coverage <- matrix(0, nrow = B, ncol = 5)
  volume <- matrix(0, nrow = B, ncol = 5)
  bias_bm <- matrix(0, nrow = B, ncol = p^2)
  bias_wbm <- matrix(0, nrow = B, ncol = p^2)
  bias_lug <- matrix(0, nrow = B, ncol = p^2)
  bias_cc <- matrix(0, nrow = B, ncol = p^2)
  
  
  truth.var <- true.sig(p = p, omega = omega, phi = phi)$final.cov
  
  for(b in 1:B)
  {
    if(b %% 10 ==0) print(b)
    chain <- var1(p = p, phi = phi, nsim = N, omega = omega)
    beta.mean <- colMeans(chain)
    
    foo <- mcse.multi(chain, r = 1)
    bn <- foo$size
    bm <- foo$cov
    
    wbm <- mcse.multi(chain, r = 2)$cov
    lug3 <- mcse.multi(chain, r = 3)$cov
    cc <- cov.sig(chain)$covariance
    
    a <- floor(N/bn)
   
    
    
    volume[b,1] <- det(bm)^(1/p)
    volume[b,2] <- det(wbm)^(1/p)
    volume[b,3] <- det(lug3)^(1/p)
    volume[b,4] <- det(cc)^(1/p)
    volume[b,5] <- det(truth.var)^(1/p)
    
    bias_bm[b,] <- as.numeric(bm - truth.var)
    bias_wbm[b,] <- as.numeric(wbm - truth.var)
    bias_lug[b,] <- as.numeric(lug3 - truth.var)
    bias_cc[b,] <- as.numeric(cc - truth.var)
    
    
    
    large.crit <- qchisq(level, df = p)/N
    
    truth <- rep(0,p)
    tester <- t(beta.mean - truth) %*% qr.solve(bm) %*% (beta.mean - truth)
    # if(tester < crit) coverage[b,1] <- 1
    if(tester < large.crit) large.coverage[b,1] <- 1
    
    tester <- t(beta.mean - truth) %*% qr.solve(wbm) %*% (beta.mean - truth)
    # if(tester < crit.wbm) coverage[b,2] <- 1
    if(tester < large.crit) large.coverage[b,2] <- 1
    
    tester <- t(beta.mean - truth) %*% qr.solve(lug3) %*% (beta.mean - truth)
    # if(tester < crit.lug3) coverage[b,3] <- 1
    if(tester < large.crit) large.coverage[b,3] <- 1
    
    tester <- t(beta.mean - truth) %*% qr.solve(cc) %*% (beta.mean - truth)
    # if(tester < crit.lug3) coverage[b,3] <- 1
    if(tester < large.crit) large.coverage[b,4] <- 1
    
    tester <- t(beta.mean - truth) %*% qr.solve(truth.var) %*% (beta.mean - truth)
    # if(tester < crit.lug3) coverage[b,3] <- 1
    if(tester < large.crit) large.coverage[b,5] <- 1
    
  }
  colnames(large.coverage) <- c("BM", "WBM", "LUG3", "Cov-Corr", "Truth")
  colnames(volume) <- c("BM", "WBM", "LUG3", "Cov_Corr", "Truth")
  return(list("large" = large.coverage, "determinant" = volume, 
              "bias_bm" = bias_bm, "bias_wbm" = bias_wbm, "bias_lug" = bias_lug,
              "bias_cc" = bias_cc))
}
