############################################################
#Package dependencies : coda, doParallel, foreach, rstan, Rcpp, RcppArmadillo

#For windows users: Rtools is also required for compliing the cpp files

#rstan version should be at least 2.17.3

#function beta.true and noise.true are required for all simulation studies (line 44 - 52)
#function poly_col is required for simulation studies related with BGL and BSGL models

#Stan requires quite long waiting time for each simulation study. Thus, it will 
#not be inclueded in the default setting. Users can change indicator variables for 
#HMC algorithms as TRUE to run HMC. 

#If one wants to include HMC algorithms, the number of iterations should be  
#increased to at least 5000 to get meaningful results since it needs a reasonable 
#amout of burn-in.

#One recommended example : n =  50, BSGL model 
# First run line 33--66, then run line 391--662

#Please set working directory to the souce file location.
#All plots will be saved in the souce file.
############################################################
install.packages("gfortran")
#install.packages('coda')
#install.packages('doParallel')
#install.packages('foreach')
#install.packages('rstan')
#install.packages('Rcpp')
#install.packages('RcppArmadillo')
#install.packages('mcmcse')

#Indicator for whether to include HMC algorithm for the bgl model
Indicator_hmc_bgl <- FALSE
#Indicator for whether to include HMC algorithm for the bsgl model
Indicator_hmc_bsgl <- FALSE
#Indicator for whether to include HMC algorithm for the bfl model
Indicator_hmc_bfl <- FALSE
#Indicator for whether to include real data analysis
Indicator_real_data <- FALSE
library(foreach)
library(doParallel)
# "True" beta and noise (for generating Y)
beta.true <- function(){
  temp <- rep(0,(p*5))
  p.nz <- ceiling(5*p*s)
  if(p.nz>0){
    temp[1:p.nz] <- rt(p.nz,2)
  }
  return(temp)
}
noise.true <- function(){rnorm(n)}
# Build design matrix by 5 order polynomial expansion 
poly_col <- function(X){
  n <- dim(X)[1]
  p <- dim(X)[2]
  matrix_out <- matrix(NA, nrow = n, ncol = (5*p))
  ind_start <- 1
  for(i in 1:p)
  {
    vec <- X[, i]
    matrix_out[, ind_start:(ind_start+4)] <- cbind(vec, vec^2, vec^3, vec^4, vec^5)
    ind_start <- ind_start + 5
  }
  return(matrix_out)
}
################################################################################
# Generate MCMC output for Figure 1(left), Figure 2(left), Figure 9(left),
# Figure 10(left) and Figure 15
################################################################################
# Random seed and dependence 
set.seed(314)
library(coda)

# Import 2BG and 3BG samplers : bgl_original performs 3BG, bgl performs 2BG
# defalut value of lambda for these two functions is 1
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload the computer
clusterEvalQ(cl, {
  Rcpp::sourceCpp('bgl.cpp') })
registerDoParallel(cl)

# n and k(denoted as p_factor) settings
n <- 50
p_factor <- c(5:9, seq(10, 50, 10))

# Sparsity of "true" coefficients for generating data
s <- .2

# Number of MCMC iterations per chain
K <- 100

# Number of chains : please set M be an integer at least 2 
M <- 5

# Generate MCMC output by 3BG for the bgl model 
lag1ac_3bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_3bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_3bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
for (j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  group_size <- rep(5, n_group)
  res_inner = foreach(k = 1:M, .noexport = c("bgl_original"), .packages = c('coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,(p*5))
    Y.raw <- drop(X%*%beta.true()+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, (p*5))
    sigma2 <- 1
    start_time <- proc.time()
    res <- bgl_original(X, Y, group_size, beta, sigma2, K = K, lambda = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1:(K * 0.1))])))
  }
  time_3bg_bgl[j, ] <- res_inner[,1]
  lag1ac_3bg_bgl[j, ] <- res_inner[,2]
  eff_3bg_bgl[j, ] <- res_inner[,3]
}

# Generate MCMC output by 2BG for the bgl model 
lag1ac_2bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_2bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_2bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
for (j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  group_size <- rep(5, n_group)
  res_inner = foreach(k = 1:M, .noexport = c("bgl"), .packages = c('coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,(p*5))
    Y.raw <- drop(X%*%beta.true()+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, (p*5))
    sigma2 <- 1
    start_time <- proc.time()
    res <- bgl(X, Y, group_size, beta, sigma2, K = K, lambda = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1:(K * 0.1))])))
  }
  time_2bg_bgl[j, ] <- res_inner[,1]
  lag1ac_2bg_bgl[j, ] <- res_inner[,2]
  eff_2bg_bgl[j, ] <- res_inner[,3]
}

# Generate HMC output by Stan for the bgl model based prior (2.5)
# Random seed and dependence 
if(Indicator_hmc_bgl){
set.seed(314)
lambda <- 1
lag1ac_hmc_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_hmc_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_hmc_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
for(j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  res_inner = foreach (k = 1:M, .packages = c('rstan', 'coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    beta_holder <- beta.true()
    Y.raw <- drop(X%*%beta_holder+noise.true())
    Y <- Y.raw-mean(Y.raw)
    stan_dat <- list(n = n,
                     p = (p*5),
                     lambda = 1,
                     X = X,
                     Y = Y,
                     group_size = rep(5, n_group),
                     n_group = n_group)
    init_f <- function () list(beta = rep(1, p*5), sigma = 1)
    fit_stan <- stan(file = 'bglasso_v2.stan', data = stan_dat, init = init_f,
                     iter = K, chains = 1, warmup = K * 0.1, refresh = 0)
    mat_sigmas <- extract(fit_stan, 'sigma' ,permuted = FALSE)
    mat_sigmas <- mat_sigmas[,,1]
    c(time = get_elapsed_time(fit_stan)[,2],
      lag1ac = autocorr(as.mcmc(mat_sigmas))[2],
      eff = effectiveSize(as.mcmc(mat_sigmas)))
  }
  time_hmc_bgl[j, ] <- res_inner[, 1]
  lag1ac_hmc_bgl[j, ] <- res_inner[, 2]
  eff_hmc_bgl[j, ] <- res_inner[, 3]
}

# Generate HMC output by Stan for the bgl model based  hierarchical prior (2.6)

# Random seed and dependence
set.seed(314)
lambda <- 1
lag1ac_hmc_hier_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_hmc_hier_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_hmc_hier_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
for(j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  res_inner = foreach (k = 1:M, .packages = c('rstan', 'coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    beta_holder <- beta.true()
    Y.raw <- drop(X%*%beta_holder+noise.true())
    Y <- Y.raw-mean(Y.raw)
    stan_dat <- list(n = n,
                     p = (p*5),
                     lambda = 1,
                     X = X,
                     Y = Y,
                     group_size = rep(5, n_group),
                     n_group = n_group)
    init_f <- function () list(beta = rep(1, p*5), sigma = 1)
    fit_stan <- stan(file = 'bglasso.stan', data = stan_dat, init = init_f,
                     iter = K, chains = 1, warmup = K * 0.1, refresh = 0)
    mat_sigmas <- extract(fit_stan, 'sigma' ,permuted = FALSE)
    mat_sigmas <- mat_sigmas[,,1]
    c(time = get_elapsed_time(fit_stan)[,2],
      lag1ac = autocorr(as.mcmc(mat_sigmas))[2],
      eff = effectiveSize(as.mcmc(mat_sigmas)))
  }
  time_hmc_hier_bgl[j, ] <- res_inner[, 1]
  lag1ac_hmc_hier_bgl[j, ] <- res_inner[, 2]
  eff_hmc_hier_bgl[j, ] <- res_inner[, 3]
}
}
stopCluster(cl)
################################################################################
# Generate plots for Figure 1(left), Figure 2(left), Figure 9(left), Figure 10(left)
# and Figure 15 using MCMC output generated by code above
################################################################################
# Left of figure 1 : bgl 50 ac plot
if(Indicator_hmc_bgl){
  npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
  stan_ac <- rowMeans(lag1ac_hmc_bgl)
  stan_hier_ac <- rowMeans(lag1ac_hmc_hier_bgl)   
  fast_ac <- rowMeans(lag1ac_2bg_bgl)
  original_ac <- rowMeans(lag1ac_3bg_bgl)
  pdf("figure1_left.pdf")
  plot(1:10,fast_ac,pch=1,type="b",ylim=c(0,1),
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
  axis(1, at=1:10, labels=npratio)
  points(1:10,original_ac,pch=2,type="b")
  points(1:10, stan_ac, pch = 3, type='b')
  points(1:10, stan_hier_ac, pch = 4, type='b')
  legend(1.032*max(2.3),(1+.75)/2,xjust=1,yjust=0.2,pch=c(1, 2, 3, 4),
         c("2BG","3BG", "Stan", "Stan_hier"),cex=.9)
  dev.off()
}
if(!Indicator_hmc_bgl){
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
fast_ac <- rowMeans(lag1ac_2bg_bgl)
original_ac <- rowMeans(lag1ac_3bg_bgl)
pdf("figure1_left.pdf")
plot(1:10,fast_ac,pch=1,type="b",ylim=c(0,1),
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
axis(1, at=1:10, labels=npratio)
points(1:10,original_ac,pch=2,type="b")
legend(1.032*max(2.3),(1+.75)/2,xjust=1,yjust=0.2,pch=c(1, 2),
       c("2BG","3BG"),cex=.9)
dev.off()
}
# Left of figure 2 : bgl 50 eff plot 
if(Indicator_hmc_bgl){
  npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
  eff_bgl50_stan <- eff_hmc_bgl / time_hmc_bgl
  eff_bgl50_stan_hier <- eff_hmc_hier_bgl / time_hmc_hier_bgl
  eff_bgl50_2bg <- eff_2bg_bgl / (time_2bg_bgl * 0.9)
  eff_bgl50_3bg <- eff_3bg_bgl / (time_3bg_bgl * 0.9)
  eff_stan <- log10(rowMeans(eff_bgl50_stan))
  eff_stan_hier <- log10(rowMeans(eff_bgl50_stan_hier))
  # 0.9 * the total time for 20000 iterations (first 2000 iterations Burin in)
  eff_o <- log10(rowMeans(eff_bgl50_3bg))
  eff_f <- log10(rowMeans(eff_bgl50_2bg))
  pdf("figure2_left.pdf")
  plot(1:10,eff_f,pch=1,type="b",ylim=c(-3,5), yaxt = "n",
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
  axis(1, at=1:10, labels=npratio)
  axis(2, at=-3:5, labels=c(-3:5))
  points(1:10,eff_o,pch=2,type="b")
  points(1:10,eff_stan,pch=3,type="b")
  points(1:10,eff_stan_hier,pch=4,type="b")
  legend(1.032*max(10),(.6 + 3.4),xjust=1,yjust=0.2,pch=c(1, 2, 3, 4),
         c("2BG","3BG", "Stan", "Stan_hier"),cex=.9) 
  dev.off()
}
if(!Indicator_hmc_bgl){
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
eff_bgl50_2bg <- eff_2bg_bgl / (time_2bg_bgl * 0.9)
eff_bgl50_3bg <- eff_3bg_bgl / (time_3bg_bgl * 0.9)
# 0.9 * the total time for 20000 iterations (first 2000 iterations Burin in)
eff_o <- log10(rowMeans(eff_bgl50_3bg))
eff_f <- log10(rowMeans(eff_bgl50_2bg))
pdf("figure2_left.pdf")
plot(1:10,eff_f,pch=1,type="b",ylim=c(-3,5), yaxt = "n",
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
axis(1, at=1:10, labels=npratio)
axis(2, at=-3:5, labels=c(-3:5))
points(1:10,eff_o,pch=2,type="b")
legend(1.032*max(10),(.6 + 3.4),xjust=1,yjust=0.2,pch=c(1, 2),
       c("2BG","3BG"),cex=.9) 
dev.off()
}
# Left of figure 9 : bgl 50 boxplot for estimates of rho_1
if(Indicator_hmc_bgl){
  pdf("figure9_left.pdf")
  boxplot( t(lag1ac_3bg_bgl)  ~ col( t(lag1ac_3bg_bgl)), xaxt = 'n', xlab='', ylab='',
           ylim = c(0, 1), col = rgb(1,0.5,0.3, alpha=0.5) )
  
  boxplot( t(lag1ac_hmc_bgl) ~ col( t(lag1ac_hmc_bgl)  ), xaxt = 'n', xlab='', ylab='',
           ylim = c(0, 1), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE )
  
  boxplot( t(lag1ac_2bg_bgl)  ~ col( t(lag1ac_2bg_bgl) ), xaxt = 'n', xlab = '', ylab='',
           ylim = c(0, 1), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE )
  
  boxplot( t(lag1ac_hmc_hier_bgl) ~ col( t(lag1ac_hmc_hier_bgl)  ), xaxt = 'n', xlab='', ylab='',
           ylim = c(0, 1), col = rgb(.6,0.1,0.6, alpha=0.5), add = TRUE )
  axis(1, at=1:10, labels=npratio)
  title(xlab = expression(frac(italic(p), italic(n))), ylab = 'Autocorrelation')
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5),
               rgb(.6,0.1,0.6, alpha=0.5))
  legend("topleft", inset=.02, 
         c("2BG","3BG","Stan", "Stan_hier"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
  # Left of figure 10 : bgl 50 boxplot for estimates of effiency
  pdf("figure10_left.pdf")
  boxplot( t(log10(eff_bgl50_3bg))  ~ col( t(log10(eff_bgl50_3bg)) ), xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  
  boxplot( t(log10(eff_bgl50_stan))  ~ col( t(log10(eff_bgl50_stan))  ), xaxt = 'n', xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE )
  
  boxplot( t(log10(eff_bgl50_2bg))  ~ col( t(log10(eff_bgl50_2bg))  ), xaxt = 'n', xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE )
  
  boxplot( t(log10(eff_bgl50_stan_hier)) ~ col( t(log10(eff_bgl50_stan_hier))), xaxt = 'n', xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(.6,0.1,0.6, alpha=0.5), add = TRUE )
  axis(1, at=1:10, labels=npratio)
  axis(2, at=-3:5, labels=c(-3:5))
  title(ylab="Efficiency (in log scale)")
  title(xlab = expression(frac(italic(p), italic(n))))
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5),
               rgb(.6,0.1,0.6, alpha=0.5))
  legend("topright", inset=.02, 
         c("2BG","3BG","Stan", "Stan_hier"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
}
if(!Indicator_hmc_bgl){
pdf("figure9_left.pdf")
boxplot( t(lag1ac_3bg_bgl)  ~ col( t(lag1ac_3bg_bgl)), xaxt = 'n', xlab = '', ylab='',
         ylim = c(0, 1), col = rgb(1,0.5,0.3, alpha=0.5) )
boxplot( t(lag1ac_2bg_bgl)  ~ col( t(lag1ac_2bg_bgl) ), xaxt = 'n', xlab = '', ylab='',
         ylim = c(0, 1), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE )
axis(1, at=1:10, labels=npratio)
title(xlab = expression(frac(italic(p), italic(n))), ylab = 'Autocorrelation')
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topleft", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()
# Left of figure 10 : bgl 50 boxplot for estimates of effiency
pdf("figure10_left.pdf")
boxplot( t(log10(eff_bgl50_3bg))  ~ col( t(log10(eff_bgl50_3bg)) ), xlab='', ylab='',
         ylim = c(-3, 5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )

boxplot( t(log10(eff_bgl50_2bg))  ~ col( t(log10(eff_bgl50_2bg))  ), xaxt = 'n', xlab = '', ylab='',
         ylim = c(-3, 5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE )
axis(1, at=1:10, labels=npratio)
axis(2, at=-3:5, labels=c(-3:5))
title(ylab="Efficiency (in log scale)")
title(xlab = expression(frac(italic(p), italic(n))))
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topright", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()
}
###############################################################################
# Generate MCMC output for Figure 3(left), Figure 4(left), Figure 11(left),
# Figure 12(left) and Figure 17
###############################################################################
# Random seed 
set.seed(314)

# Import 2BG and 3BG samplers : bsgl_original performs 3BG, bsgl performs 2BG
# defalut value of lambdas for these two functions is (1, 1)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
clusterEvalQ(cl, {
  Rcpp::sourceCpp('bsgl.cpp') })
registerDoParallel(cl)

# n and k(denoted as p_factor) settings
n <- 50
p_factor <- c(5:9, seq(10, 50, 10))

# Sparsity of "true" coefficients for generating data
s <- .2

# Number of MCMC iterations per chain
K <- 100

# Number of chains 
M <- 5

# Generate MCMC output by 3BG for the bsgl model
lag1ac_3bg_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_3bg_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_3bg_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)

for (j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  group_size <- rep(5, n_group)
  res_inner = foreach(k = 1:M, .noexport = c("bsgl_original"), .packages = c('coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,(p*5))
    Y.raw <- drop(X%*%beta.true()+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, (p*5))
    sigma2 <- 1
    start_time <- proc.time()
    res <- bsgl_original(X, Y, group_size, beta, sigma2, K = K, lambda1 = 1, lambda2 = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1: (K * 0.1))])))
  }
  time_3bg_bsgl[j, ] <- res_inner[,1]
  lag1ac_3bg_bsgl[j, ] <- res_inner[,2]
  eff_3bg_bsgl[j, ] <- res_inner[,3]
}

# Generate MCMC output by 2BG for the bsgl model
lag1ac_2bg_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_2bg_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_2bg_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)

for (j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  group_size <- rep(5, n_group)
  res_inner = foreach(k = 1:M, .noexport = c("bsgl"), .packages = c('coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,(p*5))
    Y.raw <- drop(X%*%beta.true()+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, (p*5))
    sigma2 <- 1
    start_time <- proc.time()
    res <- bsgl(X, Y, group_size, beta, sigma2, K = K, lambda1 = 1, lambda2 = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1: (K * 0.1))])))
  }
  time_2bg_bsgl[j, ] <- res_inner[,1]
  lag1ac_2bg_bsgl[j, ] <- res_inner[,2]
  eff_2bg_bsgl[j, ] <- res_inner[,3]
}

# Generate HMC output by Stan for the bsgl model based prior (2.11)

# Random seed and dependence
if(Indicator_hmc_bsgl){
set.seed(314)

lambda1 <- 1
lambda2 <- 1

lag1ac_hmc_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_hmc_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_hmc_bsgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)

for(j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  res_inner = foreach (k = 1:M, .packages = c('rstan', 'coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    beta_holder <- beta.true()
    Y.raw <- drop(X%*%beta_holder+noise.true())
    Y <- Y.raw-mean(Y.raw)
    stan_dat <- list(n = n,
                     p = p*5,
                     lambda1 = 1,
                     lambda2 = 1,
                     X = X,
                     Y = Y,
                     group_size = rep(5, n_group),
                     n_group = n_group)
    init_f <- function () list(beta = rep(1, p*5), sigma = 1)
    fit_stan <- stan(file = 'bsglasso.stan', data = stan_dat, init = init_f,
                     iter = K, chains = 1, warmup = K * 0.1, refresh = 0)
    mat_sigmas <- extract(fit_stan, 'sigma' ,permuted = FALSE)
    mat_sigmas <- mat_sigmas[,,1]
    c(time = get_elapsed_time(fit_stan)[,2],
      lag1ac = autocorr(as.mcmc(mat_sigmas))[2],
      eff = effectiveSize(as.mcmc(mat_sigmas)))
  }
  time_hmc_bsgl[j, ] <- res_inner[, 1]
  lag1ac_hmc_bsgl[j, ] <- res_inner[, 2]
  eff_hmc_bsgl[j, ] <- res_inner[, 3]
}
}
stopCluster(cl)
################################################################################
# Generate plots for Figure 3(left), Figure 4(left), Figure 11(left),
# Figure 12(left) and Figure 17
# using MCMC output generated by code above
################################################################################
if(Indicator_hmc_bsgl){
  npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
  stan_ac <- rowMeans(lag1ac_hmc_bsgl)
  fast_ac <- rowMeans(lag1ac_2bg_bsgl)
  original_ac <- rowMeans(lag1ac_3bg_bsgl)
  pdf("figure3_left.pdf")
  plot(1:10,fast_ac,pch=1,type="b",ylim=c(0,1),
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
  axis(1, at=1:10, labels=npratio)
  points(1:10,original_ac,pch=2,type="b")
  points(1:10, stan_ac, pch = 3, type='b')
  legend(1.032*max(1.8),(1+.8)/2,xjust=1,yjust=0.2,pch=c(1,2, 3),
         c("2BG","3BG", "Stan"),cex=.9)
  dev.off()
  # Left of figure 4 : bsgl 50 eff plot
  npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
  eff_bsgl50_stan <- eff_hmc_bsgl / time_hmc_bsgl
  eff_bsgl50_2bg <- eff_2bg_bsgl / (time_2bg_bsgl * 0.9)
  eff_bsgl50_3bg <- eff_3bg_bsgl / (time_3bg_bsgl * 0.9)
  eff_stan <- log10(rowMeans(eff_bsgl50_stan))
  # 0.9 * the total time for 20000 iterations (first 2000 iterations Burin in)
  eff_o <- log10(rowMeans(eff_bsgl50_3bg))
  eff_f <- log10(rowMeans(eff_bsgl50_2bg))
  pdf("figure4_left.pdf")
  plot(1:10,eff_f,pch=1,type="b",ylim=c(-3,5), yaxt = "n",
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
  axis(1, at=1:10, labels=npratio)
  axis(2, at=-3:5, labels=c(-3:5))
  points(1:10,eff_o,pch=2,type="b")
  points(1:10,eff_stan,pch=3,type="b")
  legend(1.032*max(10),(.6 + 3.4),xjust=1,yjust=0.2,pch=c(1, 2, 3),
         c("2BG","3BG", "Stan"),cex=.9) 
  dev.off()
  # Left of figure 11 : bsgl 50 boxplot for estimates of rho_1
  pdf("figure11_left.pdf")
  boxplot( t(lag1ac_3bg_bsgl)  ~ col( t(lag1ac_3bg_bsgl) ), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  
  boxplot( t(lag1ac_hmc_bsgl)  ~ col( t(lag1ac_hmc_bsgl)), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE , xaxt = 'n')
  
  boxplot( t(lag1ac_2bg_bsgl) ~ col( t(lag1ac_2bg_bsgl) ), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
  npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
  axis(1, at=1:10, labels=npratio)
  title(xlab = expression(frac(italic(p), italic(n))), ylab = 'Autocorrelation')
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5))
  legend("topleft", inset=.02, 
         c("2BG","3BG","Stan"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
  # Left of figure 12 : bsgl 50 boxplot for estimates of efficiency
  pdf("figure12_left.pdf")
  boxplot( t(log10(eff_bsgl50_3bg))  ~ col( t(log10(eff_bsgl50_3bg)) ), xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  
  boxplot( t(log10(eff_bsgl50_stan))  ~ col( t(log10(eff_bsgl50_stan))  ), xaxt = 'n', xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE )
  
  boxplot( t(log10(eff_bsgl50_2bg))  ~ col( t(log10(eff_bsgl50_2bg))  ), xaxt = 'n', xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE )
  
  npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
  axis(1, at=1:10, labels=npratio)
  axis(2, at=-3:5, labels=c(-3:5))
  title(ylab="Efficiency (in log scale)")
  title(xlab = expression(frac(italic(p), italic(n))))
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5))
  legend("topright", inset=.02, 
         c("2BG","3BG","Stan"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
}
if(!Indicator_hmc_bsgl){
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
fast_ac <- rowMeans(lag1ac_2bg_bsgl)
original_ac <- rowMeans(lag1ac_3bg_bsgl)
pdf("figure3_left.pdf")
plot(1:10,fast_ac,pch=1,type="b",ylim=c(0,1),
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
axis(1, at=1:10, labels=npratio)
points(1:10,original_ac,pch=2,type="b")
legend(1.032*max(1.8),(1+.8)/2,xjust=1,yjust=0.2,pch=c(1,2),
       c("2BG","3BG"),cex=.9)
dev.off()
# Left of figure 4 : bsgl 50 eff plot
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
eff_bsgl50_2bg <- eff_2bg_bsgl / (time_2bg_bsgl * 0.9)
eff_bsgl50_3bg <- eff_3bg_bsgl / (time_3bg_bsgl * 0.9)
# 0.9 * the total time for 20000 iterations (first 2000 iterations Burin in)
eff_o <- log10(rowMeans(eff_bsgl50_3bg))
eff_f <- log10(rowMeans(eff_bsgl50_2bg))
pdf("figure4_left.pdf")
plot(1:10,eff_f,pch=1,type="b",ylim=c(-3,5), yaxt = "n",
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
axis(1, at=1:10, labels=npratio)
axis(2, at=-3:5, labels=c(-3:5))
points(1:10,eff_o,pch=2,type="b")
legend(1.032*max(10),(.6 + 3.4),xjust=1,yjust=0.2,pch=c(1, 2),
       c("2BG","3BG"),cex=.9) 
dev.off()
# Left of figure 11 : bsgl 50 boxplot for estimates of rho_1
pdf("figure11_left.pdf")
boxplot( t(lag1ac_3bg_bsgl)  ~ col( t(lag1ac_3bg_bsgl) ), xlab='', ylab='',
         ylim = c(0, 1.1), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )

boxplot( t(lag1ac_2bg_bsgl) ~ col( t(lag1ac_2bg_bsgl) ), xlab='', ylab='',
         ylim = c(0, 1.1), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
axis(1, at=1:10, labels=npratio)
title(xlab = expression(frac(italic(p), italic(n))), ylab = 'Autocorrelation')
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topleft", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()
# Left of figure 12 : bsgl 50 boxplot for estimates of efficiency
pdf("figure12_left.pdf")
boxplot( t(log10(eff_bsgl50_3bg))  ~ col( t(log10(eff_bsgl50_3bg)) ), xlab='', ylab='',
         ylim = c(-3, 5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )

boxplot( t(log10(eff_bsgl50_2bg))  ~ col( t(log10(eff_bsgl50_2bg))  ), xaxt = 'n', xlab='', ylab='',
         ylim = c(-3, 5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE )

npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
axis(1, at=1:10, labels=npratio)
axis(2, at=-3:5, labels=c(-3:5))
title(ylab="Efficiency (in log scale)")
title(xlab = expression(frac(italic(p), italic(n))))
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topright", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()
}

################################################################################
# Generate MCMC output for Figure 5(left), Figure 6(left), Figure 13(left),
# Figure 14(left) and Figure 19
################################################################################
# Random seed 
set.seed(314)

# Import 2BG and 3BG samplers : bfl_original performs 3BG, bfl performs 2BG
# defalut value of lambdas for these two functions is (1, 1)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
clusterEvalQ(cl, {
  Rcpp::sourceCpp('bfl.cpp') })
registerDoParallel(cl)

# n and p/n settings
n <- 100
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))

# Sparsity of "true" coefficients for generating data
s <- .2

# Approximate pairwise correlation of covariates
r <- .2

# Number of MCMC iterations per chain
K <- 50

# Number of chains 
M <- 5

# Generate MCMC output by 3BG for the bfl model
lag1ac_3bg_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)
eff_3bg_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)
time_3bg_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)

for (j in 1:length(npratio))
{
  p <- n * npratio[j]
  res_inner = foreach(k = 1:M, .noexport = c("bfl_original"), .packages = c('coda'), .combine = rbind) %dopar%{
    Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,p)
    beta_holder <- c(rnorm(p/10, 1, 0.1), rep(0, p/5), rnorm(p/10, 1, 0.1), rep(0, (3*p)/5))
    Y.raw <- drop(X%*%beta_holder+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, p)
    sigma2 <- 1
    start_time <- proc.time()
    res <- bfl_original(X, Y, beta, sigma2, K = K, lambda1 = 1, lambda2 = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1:(K * 0.1))])) )
  }
  time_3bg_bfl[j, ] <- res_inner[,1]
  lag1ac_3bg_bfl[j, ] <- res_inner[,2]
  eff_3bg_bfl[j, ] <- res_inner[,3]
}

# Generate MCMC output by 2BG for the bfl model
lag1ac_2bg_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)
eff_2bg_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)
time_2bg_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)

for (j in 1:length(npratio))
{
  p <- n * npratio[j]
  res_inner = foreach(k = 1:M, .noexport = c("bfl"), .packages = c('coda'), .combine = rbind) %dopar%{
    Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,p)
    beta_holder <- c(rnorm(p/10, 1, 0.1), rep(0, p/5), rnorm(p/10, 1, 0.1), rep(0, (3*p)/5))
    Y.raw <- drop(X%*%beta_holder+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, p)
    sigma2 <- 1
    start_time <- proc.time()
    res <- bfl(X, Y, beta, sigma2, K = K, lambda1 = 1, lambda2 = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1:(K * 0.1))])) )
  }
  time_2bg_bfl[j, ] <- res_inner[,1]
  lag1ac_2bg_bfl[j, ] <- res_inner[,2]
  eff_2bg_bfl[j, ] <- res_inner[,3]
}

# Generate HMC output by Stan for the bfl model based prior (2.17)

# Random seed and dependence
if(Indicator_hmc_bfl){
set.seed(314)
lag1ac_hmc_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)
time_hmc_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)
eff_hmc_bfl <- matrix(rep(NA, (length(npratio) * M)), ncol = M)
for(j in 1:length(npratio))
{
  p <- n * npratio[j]
  res_inner = foreach (k = 1:M, .packages = c('rstan', 'coda'), .combine = rbind) %dopar%{
    Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
    X <- scale(X.raw)*sqrt(n/(n-1))
    beta_holder <- c(rnorm(p/10, 1, 0.1), rep(0, p/5), rnorm(p/10, 1, 0.1), rep(0, (3*p)/5))
    Y.raw <- drop(X%*%beta_holder+noise.true())
    Y <- Y.raw-mean(Y.raw)
    stan_dat <- list(n = n,
                     p = p,
                     lambda1 = 1,
                     lambda2 = 1,
                     X = X,
                     Y = Y)
    init_f <- function () list(beta = rep(1, p), sigma = 1)
    fit_stan <- stan(file = 'bfl.stan', data = stan_dat, init = init_f,
                     iter = K, chains = 1, warmup = K * 0.1, refresh = 0)
    mat_sigmas <- extract(fit_stan, 'sigma' ,permuted = FALSE)
    mat_sigmas <- mat_sigmas[,,1]
    c(time = get_elapsed_time(fit_stan)[,2],
      lag1ac = autocorr(as.mcmc(mat_sigmas))[2],
      eff = effectiveSize(as.mcmc(mat_sigmas)))
  }
  time_hmc_bfl[j, ] <- res_inner[, 1]
  lag1ac_hmc_bfl[j, ] <- res_inner[, 2]
  eff_hmc_bfl[j, ] <- res_inner[, 3]
}
}
stopCluster(cl)
###############################################################################
# Generate plots for Figure 5(left), Figure 6(left), Figure 13(left),
# Figure 14(left) and Figure 19
# using MCMC output generated by code above
###############################################################################
# Left of figure 5 : bfl 100 ac plot
if(Indicator_hmc_bfl){
  stan_ac <- rowMeans(lag1ac_hmc_bfl)
  fast_ac <- rowMeans(lag1ac_2bg_bfl)
  original_ac <- rowMeans(lag1ac_3bg_bfl)
  pdf("figure5_left.pdf")
  plot(1:10,fast_ac,pch=1,type="b",ylim=c(0,1),
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
  axis(1, at=1:10, labels=npratio)
  points(1:10,original_ac,pch=2,type="b")
  points(1:10, stan_ac, pch = 3, type='b')
  legend(1.032*max(1.8),(1+.8)/2,xjust=1,yjust=0.2,pch=c(1,2, 3),
         c("2BG","3BG", "Stan"),cex=.9)
  dev.off()
  # Left of figure 6 : bfl 100 eff plot
  eff_bfl100_stan <- eff_hmc_bfl / time_hmc_bfl
  eff_bfl100_2bg <- eff_2bg_bfl / (time_2bg_bfl * .9)
  eff_bfl100_3bg <- eff_3bg_bfl / (time_3bg_bfl * .9)
  eff_stan <- log10(rowMeans(eff_bfl100_stan))
  # 0.9 * the total time for 20000 iterations (first 2000 iterations Burin in)
  eff_o <- log10(rowMeans(eff_bfl100_3bg))
  eff_f <- log10(rowMeans(eff_bfl100_2bg))
  pdf("figure6_left.pdf")
  plot(1:10,eff_f,pch=1,type="b",ylim=c(-3,5), yaxt = "n",
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
  axis(1, at=1:10, labels=npratio)
  axis(2, at=-3:5, labels=c(-3:5))
  points(1:10,eff_o,pch=2,type="b")
  points(1:10,eff_stan,pch=3,type="b")
  legend(1.032*max(10),(.6 + 3.4),xjust=1,yjust=0.2,pch=c(1, 2, 3),
         c("2BG","3BG", "Stan"),cex=.9) 
  dev.off()
  # Left of figure 13 : bfl 100 boxplot for estimates of rho_1
  pdf("figure13_left.pdf")
  boxplot( t(lag1ac_3bg_bfl)  ~ col( t(lag1ac_3bg_bfl) ), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  
  boxplot( t(lag1ac_hmc_bfl)  ~ col( t(lag1ac_hmc_bfl)), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE , xaxt = 'n')
  
  boxplot( t(lag1ac_2bg_bfl) ~ col( t(lag1ac_2bg_bfl) ), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
  npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
  axis(1, at=1:10, labels=npratio)
  title(xlab = expression(frac(italic(p), italic(n))), ylab = 'Autocorrelation')
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5))
  legend("topleft", inset=.02, 
         c("2BG","3BG","Stan"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
  # Left of figure 14 : bfl 100 boxplot for estimates of efficiency
  pdf("figure14_left.pdf")
  boxplot( t(log10(eff_bfl100_3bg))  ~ col( t(log10(eff_bfl100_3bg)) ), xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  
  boxplot( t(log10(eff_bfl100_stan))  ~ col( t(log10(eff_bfl100_stan))  ), xaxt = 'n', xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE )
  
  boxplot( t(log10(eff_bfl100_2bg))  ~ col( t(log10(eff_bfl100_2bg))  ), xaxt = 'n', xlab='', ylab='',
           ylim = c(-3, 5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE )
  
  npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
  axis(1, at=1:10, labels=npratio)
  axis(2, at=-3:5, labels=c(-3:5))
  title(ylab="Efficiency (in log scale)")
  title(xlab = expression(frac(italic(p), italic(n))))
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5))
  legend("topright", inset=.02, 
         c("2BG","3BG","Stan"), fill=col_vec, horiz=TRUE, cex=0.8)
  
  dev.off()
  
}
if(!Indicator_hmc_bfl){
fast_ac <- rowMeans(lag1ac_2bg_bfl)
original_ac <- rowMeans(lag1ac_3bg_bfl)
pdf("figure5_left.pdf")
plot(1:10,fast_ac,pch=1,type="b",ylim=c(0,1),
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
axis(1, at=1:10, labels=npratio)
points(1:10,original_ac,pch=2,type="b")
legend(1.032*max(1.8),(1+.8)/2,xjust=1,yjust=0.2,pch=c(1,2),
       c("2BG","3BG"),cex=.9)
dev.off()
# Left of figure 6 : bfl 100 eff plot
eff_bfl100_2bg <- eff_2bg_bfl / (time_2bg_bfl * .9)
eff_bfl100_3bg <- eff_3bg_bfl / (time_3bg_bfl * .9)
# 0.9 * the total time for 20000 iterations (first 2000 iterations Burin in)
eff_o <- log10(rowMeans(eff_bfl100_3bg))
eff_f <- log10(rowMeans(eff_bfl100_2bg))
pdf("figure6_left.pdf")
plot(1:10,eff_f,pch=1,type="b",ylim=c(-3,5), yaxt = "n",
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
axis(1, at=1:10, labels=npratio)
axis(2, at=-3:5, labels=c(-3:5))
points(1:10,eff_o,pch=2,type="b")
legend(1.032*max(10),(.6 + 3.4),xjust=1,yjust=0.2,pch=c(1, 2),
       c("2BG","3BG"),cex=.9) 
dev.off()
# Left of figure 13 : bfl 100 boxplot for estimates of rho_1
pdf("figure13_left.pdf")
boxplot( t(lag1ac_3bg_bfl)  ~ col( t(lag1ac_3bg_bfl) ), xlab='', ylab='',
         ylim = c(0, 1.1), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )

boxplot( t(lag1ac_2bg_bfl) ~ col( t(lag1ac_2bg_bfl) ), xlab='', ylab='',
         ylim = c(0, 1.1), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
axis(1, at=1:10, labels=npratio)
title(xlab = expression(frac(italic(p), italic(n))), ylab = 'Autocorrelation')
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topleft", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()
# Left of figure 14 : bfl 100 boxplot for estimates of efficiency
pdf("figure14_left.pdf")
boxplot( t(log10(eff_bfl100_3bg))  ~ col( t(log10(eff_bfl100_3bg)) ), xlab='', ylab='',
         ylim = c(-3, 5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )

boxplot( t(log10(eff_bfl100_2bg))  ~ col( t(log10(eff_bfl100_2bg))  ), xaxt = 'n', xlab='', ylab='',
         ylim = c(-3, 5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE )

npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
axis(1, at=1:10, labels=npratio)
axis(2, at=-3:5, labels=c(-3:5))
title(ylab="Efficiency (in log scale)")
title(xlab = expression(frac(italic(p), italic(n))))
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topright", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)

dev.off()
}

################################################################################
# Generate MCMC output for Figure 30 and Figure 31 
################################################################################
# Random seed and dependence
set.seed(314)
# Import 2BG and 3BG samplers : bgl_original performs 3BG, bgl performs 2BG
# defalut value of lambda for these two functions is 1
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload the computer
clusterEvalQ(cl, {
  Rcpp::sourceCpp('bgl.cpp') })
registerDoParallel(cl)

# n and k(denoted as p_factor) settings
n <- 50
p_factor <- seq(1, 10, 1) * 100

# Sparsity of "true" coefficients for generating data
s <- .2

# Number of MCMC iterations per chain
K <- 100

# Number of chains : please set M be an integer at least 2
M <- 5

# Generate MCMC output by 3BG for the bgl model
lag1ac_3bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_3bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_3bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)

for (j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  group_size <- rep(5, n_group)
  res_inner = foreach(k = 1:M, .noexport = c("bgl_original"), .packages = c('coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,(p*5))
    Y.raw <- drop(X%*%beta.true()+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, (p*5))
    sigma2 <- 1
    start_time <- proc.time()
    res <- bgl_original(X, Y, group_size, beta, sigma2, K = K, lambda = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1:(K * 0.1))])))
  }
  time_3bg_bgl[j, ] <- res_inner[,1]
  lag1ac_3bg_bgl[j, ] <- res_inner[,2]
  eff_3bg_bgl[j, ] <- res_inner[,3]
}

# Generate MCMC output by 2BG for the bgl model
lag1ac_2bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_2bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_2bg_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)

for (j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  group_size <- rep(5, n_group)
  res_inner = foreach(k = 1:M, .noexport = c("bgl"), .packages = c('coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,(p*5))
    Y.raw <- drop(X%*%beta.true()+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, (p*5))
    sigma2 <- 1
    start_time <- proc.time()
    res <- bgl(X, Y, group_size, beta, sigma2, K = K, lambda = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1:(K * 0.1))])))
  }
  time_2bg_bgl[j, ] <- res_inner[,1]
  lag1ac_2bg_bgl[j, ] <- res_inner[,2]
  eff_2bg_bgl[j, ] <- res_inner[,3]
}

# Generate MCMC output by stan for the bgl model based on prior 2.5

# Random seed and dependence
if(Indicator_hmc_bgl){
set.seed(314)

lambda <- 1
lag1ac_hmc_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
time_hmc_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
eff_hmc_bgl <- matrix(rep(NA, (length(p_factor) * M)), ncol = M)
for(j in 1:length(p_factor))
{
  p <- p_factor[j]
  n_group <- p
  res_inner = foreach (k = 1:M, .packages = c('rstan', 'coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    beta_holder <- beta.true()
    Y.raw <- drop(X%*%beta_holder+noise.true())
    Y <- Y.raw-mean(Y.raw)
    stan_dat <- list(n = n,
                     p = (p*5),
                     lambda = 1,
                     X = X,
                     Y = Y,
                     group_size = rep(5, n_group),
                     n_group = n_group)
    init_f <- function () list(beta = rep(1, p*5), sigma = 1)
    fit_stan <- stan(file = 'bglasso_v2.stan', data = stan_dat, init = init_f,
                     iter = K, chains = 1, warmup = K * 0.1, refresh = 0)
    mat_sigmas <- extract(fit_stan, 'sigma' ,permuted = FALSE)
    mat_sigmas <- mat_sigmas[,,1]
    c(time = get_elapsed_time(fit_stan)[,2],
      lag1ac = autocorr(as.mcmc(mat_sigmas))[2],
      eff = effectiveSize(as.mcmc(mat_sigmas)))
  }
  time_hmc_bgl[j, ] <- res_inner[, 1]
  lag1ac_hmc_bgl[j, ] <- res_inner[, 2]
  eff_hmc_bgl[j, ] <- res_inner[, 3]
}
}
stopCluster(cl)

###############################################################################
# Generate plots for Figure 30 and Figure 31 
# using MCMC output generated by code above
###############################################################################
# Left of figure 30 : bgl wide lag1ac
if(Indicator_hmc_bgl){
  pdf("figure30_left.pdf")
  plot(1:10,rowMeans(lag1ac_2bg_bgl),pch=1,type="b",ylim=c(0,1),
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
  axis(1, at=1:10, labels=seq(10,100,10))
  points(1:10,rowMeans(lag1ac_3bg_bgl),pch=2,type="b")
  points(1:10, rowMeans(lag1ac_hmc_bgl), pch = 3, type='b')
  legend(1.032*max(10),(1 + .3)/2,xjust=1,yjust=0.2,pch=c(1,2, 3),
         c("2BG","3BG", "Stan"),cex=.9)
  dev.off()
  # Right of figure 30 : bgl wide boxplot for estimates of rho_1
  pdf("figure30_right.pdf")
  boxplot( t(lag1ac_3bg_bgl)  ~ col( t(lag1ac_3bg_bgl) ), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  
  boxplot( t(lag1ac_hmc_bgl)  ~ col( t(lag1ac_hmc_bgl)), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE , xaxt = 'n')
  
  boxplot( t(lag1ac_2bg_bgl) ~ col( t(lag1ac_2bg_bgl) ), xlab='', ylab='',
           ylim = c(0, 1.1), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
  axis(1, at=1:10, labels=seq(10,100,10))
  title(xlab = expression(frac(italic(p), italic(n))), ylab = 'Autocorrelation')
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5))
  legend("topleft", inset=.02, 
         c("2BG","3BG","Stan"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
  
  # Left of figure 31 : bgl wide eff
  eff_bglwide_stan <- eff_hmc_bgl / time_hmc_bgl
  eff_bglwide_2bg <- eff_2bg_bgl / (time_2bg_bgl * 0.9)
  eff_bglwide_3bg <- eff_3bg_bgl / (time_3bg_bgl * 0.9)
  pdf("figure31_left.pdf")
  plot(1:10,log10(rowMeans(eff_bglwide_2bg)),pch=1,type="b",ylim=c(-4,4), yaxt = "n", 
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
  points(1:10,log10(rowMeans(eff_bglwide_3bg)),pch=2,type="b")
  points(1:10,log10(rowMeans(eff_bglwide_stan)),pch=3,type="b")
  axis(1, at=1:10, labels=seq(10,100,10))
  axis(2, at=-4:4, labels=c(-4:4))
  legend(1.032*max(10),(1 + .3)/2,xjust=1,yjust=0.2,pch=c(1,2, 3),
         c("2BG","3BG", "Stan"),cex=.9)
  dev.off()
  
  # Right of figure 31 : bgl wide boxplot for estimates of efficiency
  pdf("figure31_right.pdf")
  boxplot( t(log10(eff_bglwide_3bg) ) ~ col( t(log10(eff_bglwide_3bg)) ), xlab='', ylab='',
           ylim = c(-4, 4), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  
  boxplot( t(log10(eff_bglwide_stan ) ) ~ col( t(log10(eff_bglwide_stan  ))  ), xlab='', ylab='',
           ylim = c(-4, 4), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE , xaxt = 'n')
  
  boxplot( t(log10(eff_bglwide_2bg) )  ~ col( t(log10(eff_bglwide_2bg ))  ), xlab='', ylab='',
           ylim = c(-4, 4), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
  axis(1, at=1:10, labels=seq(10,100,10))
  axis(2, at=-4:4, labels=c(-4:4))
  title(xlab = expression(frac(italic(p), italic(n))), ylab="Efficiency (in log scale)")
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5))
  legend("topright", inset=.02, 
         c("2BG","3BG","Stan"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
  
  
}
if(!Indicator_hmc_bgl){
pdf("figure30_left.pdf")
plot(1:10,rowMeans(lag1ac_2bg_bgl),pch=1,type="b",ylim=c(0,1),
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
axis(1, at=1:10, labels=seq(10,100,10))
points(1:10,rowMeans(lag1ac_3bg_bgl),pch=2,type="b")
legend(1.032*max(10),(1 + .3)/2,xjust=1,yjust=0.2,pch=c(1,2),
       c("2BG","3BG"),cex=.9)
dev.off()
# Right of figure 30 : bgl wide boxplot for estimates of rho_1
pdf("figure30_right.pdf")
boxplot( t(lag1ac_3bg_bgl)  ~ col( t(lag1ac_3bg_bgl) ), xlab='', ylab='',
         ylim = c(0, 1.1), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )

boxplot( t(lag1ac_2bg_bgl) ~ col( t(lag1ac_2bg_bgl) ), xlab='', ylab='',
         ylim = c(0, 1.1), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
axis(1, at=1:10, labels=seq(10,100,10))
title(xlab = expression(frac(italic(p), italic(n))), ylab = 'Autocorrelation')
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topleft", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()

# Left of figure 31 : bgl wide eff
eff_bglwide_2bg <- eff_2bg_bgl / (time_2bg_bgl * 0.9)
eff_bglwide_3bg <- eff_3bg_bgl / (time_3bg_bgl * 0.9)
pdf("figure31_left.pdf")
plot(1:10,log10(rowMeans(eff_bglwide_2bg)),pch=1,type="b",ylim=c(-4,4), yaxt = "n", 
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
points(1:10,log10(rowMeans(eff_bglwide_3bg)),pch=2,type="b")
axis(1, at=1:10, labels=seq(10,100,10))
axis(2, at=-4:4, labels=c(-4:4))
legend(1.032*max(10),(1 + .3)/2,xjust=1,yjust=0.2,pch=c(1,2),
       c("2BG","3BG"),cex=.9)
dev.off()

# Right of figure 31 : bgl wide boxplot for estimates of efficiency
pdf("figure31_right.pdf")
boxplot( t(log10(eff_bglwide_3bg) ) ~ col( t(log10(eff_bglwide_3bg)) ), xlab='', ylab='',
         ylim = c(-4, 4), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )

boxplot( t(log10(eff_bglwide_2bg) )  ~ col( t(log10(eff_bglwide_2bg ))  ), xlab='', ylab='',
         ylim = c(-4, 4), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
axis(1, at=1:10, labels=seq(10,100,10))
axis(2, at=-4:4, labels=c(-4:4))
title(xlab = expression(frac(italic(p), italic(n))), ylab="Efficiency (in log scale)")
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topright", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()
}



################################################################################
# Generate MCMC output for Figure 33 and Figure 34
################################################################################
# Random seed and dependence
set.seed(314)

# Import 2BG and 3BG samplers : bgl_original performs 3BG, bgl performs 2BG
# defalut value of lambda for these two functions is 1
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload the computer
clusterEvalQ(cl, {
  Rcpp::sourceCpp('bgl.cpp') })
registerDoParallel(cl)

# n and k(denoted as p_factor) settings
n_vec <- seq(50, 500, 50)
p_factor <- 5

# Sparsity of "true" coefficients for generating data
s <- .2

# Number of MCMC iterations per chain
K <- 100

# Number of chains : please set M be an integer at least 2
M <- 5

# Generate MCMC output by 3BG for the bgl model
lag1ac_3bg_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)
eff_3bg_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)
time_3bg_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)

for (j in 1:length(n_vec))
{
  p <- p_factor
  n_group <- p_factor
  n <- n_vec[j]
  group_size <- rep(5, n_group)
  res_inner = foreach(k = 1:M, .noexport = c("bgl_original"), .packages = c('coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,(p*5))
    Y.raw <- drop(X%*%beta.true()+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, (p*5))
    sigma2 <- 1
    start_time <- proc.time()
    res <- bgl_original(X, Y, group_size, beta, sigma2, K = K, lambda = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1:(K * 0.1))])))
  }
  time_3bg_bgl[j, ] <- res_inner[,1]
  lag1ac_3bg_bgl[j, ] <- res_inner[,2]
  eff_3bg_bgl[j, ] <- res_inner[,3]
}

# Generate MCMC output by 2BG for the bgl model
lag1ac_2bg_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)
eff_2bg_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)
time_2bg_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)

for (j in 1:length(n_vec))
{
  p <- p_factor
  n_group <- p_factor
  n <- n_vec[j]
  group_size <- rep(5, n_group)
  res_inner = foreach(k = 1:M, .noexport = c("bgl"), .packages = c('coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    X <- matrix(as.vector(X),n,(p*5))
    Y.raw <- drop(X%*%beta.true()+noise.true())
    Y <- Y.raw-mean(Y.raw)
    beta <- rep(1, (p*5))
    sigma2 <- 1
    start_time <- proc.time()
    res <- bgl(X, Y, group_size, beta, sigma2, K = K, lambda = 1)
    c(time = as.double(proc.time() - start_time)[3],
      lag1ac = autocorr(as.mcmc(res$sigma2s[-c(1:(K * 0.1))]))[2],
      eff = effectiveSize(as.mcmc(res$sigma2s[-c(1:(K * 0.1))])))
  }
  time_2bg_bgl[j, ] <- res_inner[,1]
  lag1ac_2bg_bgl[j, ] <- res_inner[,2]
  eff_2bg_bgl[j, ] <- res_inner[,3]
}

# Generate MCMC output by stan for the bgl model
if(Indicator_hmc_bgl){
lambda <- 1
lag1ac_hmc_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)
time_hmc_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)
eff_hmc_bgl <- matrix(rep(NA, (length(n_vec) * M)), ncol = M)
for(j in 1:length(n_vec))
{
  p <- p_factor
  n_group <- p
  n <- n_vec[j]
  res_inner = foreach (k = 1:M, .packages = c('rstan', 'coda'), .combine = rbind) %dopar%{
    X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X.raw <- poly_col(X.raw)
    X <- scale(X.raw)*sqrt(n/(n-1))
    beta_holder <- beta.true()
    Y.raw <- drop(X%*%beta_holder+noise.true())
    Y <- Y.raw-mean(Y.raw)
    stan_dat <- list(n = n,
                     p = (p*5),
                     lambda = 1,
                     X = X,
                     Y = Y,
                     group_size = rep(5, n_group),
                     n_group = n_group)
    init_f <- function () list(beta = rep(1, p*5), sigma = 1)
    fit_stan <- stan(file = 'bglasso_v2.stan', data = stan_dat, init = init_f,
                     iter = K, chains = 1, warmup = K * 0.1, refresh = 0)
    mat_sigmas <- extract(fit_stan, 'sigma' ,permuted = FALSE)
    mat_sigmas <- mat_sigmas[,,1]
    c(time = get_elapsed_time(fit_stan)[,2],
      lag1ac = autocorr(as.mcmc(mat_sigmas))[2],
      eff = effectiveSize(as.mcmc(mat_sigmas)))
  }
  time_hmc_bgl[j, ] <- res_inner[, 1]
  lag1ac_hmc_bgl[j, ] <- res_inner[, 2]
  eff_hmc_bgl[j, ] <- res_inner[, 3]
}
}
stopCluster(cl)
###############################################################################
# Generate plots for Figure 33 and Figure 34 
# using MCMC output generated by code above
###############################################################################
# Left of figure 33 : bgl tall lag1ac
if(Indicator_hmc_bgl){
  n_vec <- seq(50, 500, 50)
  pdf("figure33_left.pdf")
  plot(1:10,rowMeans(lag1ac_2bg_bgl),pch=1,type="b",ylim=c(-0.2,.6),
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(n), italic(p))),ylab="Autocorrelation")
  axis(1, at=1:10, labels=n_vec/25)
  points(1:10,rowMeans(lag1ac_3bg_bgl),pch=2,type="b")
  points(1:10, rowMeans(lag1ac_hmc_bgl), pch = 3, type='b')
  legend(1.032*max(10),(1 + .05)/2,xjust=1,yjust=0.2,pch=c(1,2, 3),
         c("2BG","3BG", "Stan"),cex=.9)
  dev.off()
  # Right of figure 33 : bgl tall boxplot for estimates of rho_1
  pdf("figure33_right.pdf")
  boxplot( t(lag1ac_3bg_bgl)  ~ col( t(lag1ac_3bg_bgl) ), xlab='', ylab='',
           ylim = c(-0.3, .5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  boxplot( t(lag1ac_2bg_bgl) ~ col( t(lag1ac_2bg_bgl) ), xlab='', ylab='',
           ylim = c(-.3, .5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
  boxplot( t(lag1ac_hmc_bgl)  ~ col( t(lag1ac_hmc_bgl)), xlab='', ylab='',
           ylim = c(-0.3, .5), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE , xaxt = 'n')
  axis(1, at=1:10, labels=n_vec/25)
  title(xlab = expression(frac(italic(n), italic(p))), ylab = 'Autocorrelation')
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5))
  legend("topright", inset=.02, 
         c("2BG","3BG","Stan"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
  # Left of figure 34: bgl tall eff
  eff_bgltall_2bg <- eff_2bg_bgl / (time_2bg_bgl * .9)
  eff_bgltall_3bg <- eff_3bg_bgl / (time_3bg_bgl * .9)
  eff_bgltall_stan <- eff_hmc_bgl / time_hmc_bgl
  
  pdf("figure34_left.pdf")
  plot(1:10,log10(rowMeans(eff_bgltall_2bg)),pch=1,type="b",ylim=c(-1,5), yaxt = "n", 
       cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
       xlab=expression(frac(italic(n), italic(p))),ylab="Efficiency (in log scale)")
  points(1:10,log10(rowMeans(eff_bgltall_3bg)),pch=2,type="b")
  points(1:10,log10(rowMeans(eff_bgltall_stan)),pch=3,type="b")
  axis(1, at=1:10, labels=n_vec / 25)
  axis(2, at=-1:5, labels=c(-1:5))
  legend(1.032*max(10),(1 + 7.9)/2,xjust=1,yjust=0.2,pch=c(1,2, 3),
         c("2BG","3BG", "Stan"),cex=.9)
  dev.off()
  
  # Right of figure 34: bgl tall boxplot for estimates of efficiency
  pdf("figure34_right.pdf")
  boxplot( t(log10(eff_bgltall_3bg ) ) ~ col( t(log10(eff_bgltall_3bg  )) ), xlab='', ylab='',
           ylim = c(-1, 5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
  
  boxplot( t(log10(eff_bgltall_stan) ) ~ col( t(log10(eff_bgltall_stan ))  ), xlab='', ylab='',
           ylim = c(-1, 5), col = rgb(.1,0.9,0.9, alpha=0.5), add = TRUE , xaxt = 'n')
  
  boxplot( t(log10(eff_bgltall_2bg ) )  ~ col( t(log10(eff_bgltall_2bg ))  ), xlab='', ylab='',
           ylim = c(-1, 5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
  axis(1, at=1:10, labels=n_vec / 25)
  axis(2, at=-1:4, labels=c(-1:4))
  title(xlab = expression(frac(italic(p), italic(n))), ylab="Efficiency (in log scale)")
  col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5), rgb(.1,0.9,0.9, alpha=0.5))
  legend("topright", inset=.02, 
         c("2BG","3BG","Stan"), fill=col_vec, horiz=TRUE, cex=0.8)
  dev.off()
  
  
}
if(!Indicator_hmc_bgl){
n_vec <- seq(50, 500, 50)
pdf("figure33_left.pdf")
plot(1:10,rowMeans(lag1ac_2bg_bgl),pch=1,type="b",ylim=c(-0.2,.6),
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(n), italic(p))),ylab="Autocorrelation")
axis(1, at=1:10, labels=n_vec/25)
points(1:10,rowMeans(lag1ac_3bg_bgl),pch=2,type="b")
legend(1.032*max(10),(1 + .05)/2,xjust=1,yjust=0.2,pch=c(1,2),
       c("2BG","3BG"),cex=.9)
dev.off()
# Right of figure 33 : bgl tall boxplot for estimates of rho_1
pdf("figure33_right.pdf")
boxplot( t(lag1ac_3bg_bgl)  ~ col( t(lag1ac_3bg_bgl) ), xlab='', ylab='',
         ylim = c(-0.3, .5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )
boxplot( t(lag1ac_2bg_bgl) ~ col( t(lag1ac_2bg_bgl) ),
         ylim = c(-.3, .5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
axis(1, at=1:10, labels=n_vec/25)
title(xlab = expression(frac(italic(n), italic(p))), ylab = 'Autocorrelation')
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topright", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()
# Left of figure 34: bgl tall eff
eff_bgltall_2bg <- eff_2bg_bgl / (time_2bg_bgl * .9)
eff_bgltall_3bg <- eff_3bg_bgl / (time_3bg_bgl * .9)

pdf("figure34_left.pdf")
plot(1:10,log10(rowMeans(eff_bgltall_2bg)),pch=1,type="b",ylim=c(-1,5), yaxt = "n", 
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(n), italic(p))),ylab="Efficiency (in log scale)")
points(1:10,log10(rowMeans(eff_bgltall_3bg)),pch=2,type="b")
axis(1, at=1:10, labels=n_vec / 25)
axis(2, at=-1:5, labels=c(-1:5))
legend(1.032*max(10),(1 + 7.9)/2,xjust=1,yjust=0.2,pch=c(1,2, 3),
       c("2BG","3BG", "Stan"),cex=.9)
dev.off()

# Right of figure 34: bgl tall boxplot for estimates of efficiency
pdf("figure34_right.pdf")
boxplot( t(log10(eff_bgltall_3bg ) ) ~ col( t(log10(eff_bgltall_3bg  )) ), xlab='', ylab='',
         ylim = c(-1, 5), col = rgb(1,0.5,0.3, alpha=0.5), xaxt = 'n' )

boxplot( t(log10(eff_bgltall_2bg ) )  ~ col( t(log10(eff_bgltall_2bg ))  ), xlab='', ylab='',
         ylim = c(-1, 5), col = rgb(.6,0.6,0.6, alpha=0.5), add = TRUE, xaxt = 'n' )
axis(1, at=1:10, labels=n_vec / 25)
axis(2, at=-1:4, labels=c(-1:4))
title(xlab = expression(frac(italic(p), italic(n))), ylab="Efficiency (in log scale)")
col_vec <- c(rgb(.6,0.6,0.6, alpha=0.5), rgb(1,0.5,0.3, alpha=0.5))
legend("topright", inset=.02, 
       c("2BG","3BG"), fill=col_vec, horiz=TRUE, cex=0.8)
dev.off()
}



############################################################
# Generate MCMC output for Subsection 5.1
# (Gene Expression Data)
############################################################
if(Indicator_real_data){
Rcpp::sourceCpp('bgl.cpp') 
library(gglasso)
library(coda)
library(rstan)
data(bardet)
group_size <- rep(5, 20)
X.raw <- bardet$x
Y.raw <- bardet$y
n <- dim(X.raw)[1]
p <- dim(X.raw)[2]
X <- scale(X.raw)*sqrt(n/(n-1))
X <- matrix(as.vector(X),n,p)
Y <- Y.raw-mean(Y.raw)
start_time <- proc.time()
res_2bg_bgl <- bgl(X, Y,group_size, rep(1,p), 1, lambda = 0.0601 ,K = 5000)
time_2bg_gene <- as.double(proc.time() - start_time)[3]
start_time <- proc.time()
res_3bg_bgl <- bgl_original(X, Y,group_size, rep(1,p), 1, lambda = 0.0601 ,K = 5000)
time_3bg_gene <- as.double(proc.time() - start_time)[3]
stan_dat <- list(n = n,
                 p = p,
                 lambda = 0.0601,
                 X = X,
                 Y = Y,
                 group_size = group_size,
                 n_group = 20)
fit_hmc_bgl <- stan(file = 'bglasso_v2.stan', data = stan_dat, control = list(max_treedepth = 15),
                    iter = 5000, chains = 1, warmup = 500)
time_hmc_gene <- get_elapsed_time(fit_hmc_bgl)[,2]
t <- stan_ac(fit_hmc_bgl, pars = 'sigma')
t1 <- t$data

# Autocorrelation
t1$ac
autocorr(as.mcmc(res_2bg_bgl$sigma2s[-c(1:500)]))
autocorr(as.mcmc(res_3bg_bgl$sigma2s[-c(1:500)]))

# Time
time_2bg_gene * .9
time_3bg_gene * .9
time_hmc_gene

# effective sample size
effectiveSize(as.mcmc(res_2bg_bgl$sigma2s[-c(1:500)]))
effectiveSize(as.mcmc(res_3bg_bgl$sigma2s[-c(1:500)]))
summary(fit_hmc_bgl)$summary['sigma','n_eff']

############################################################
# Generate MCMC output for Subsection 5.2
# (Economic Data)
############################################################
Rcpp::sourceCpp('bsgl.cpp') 
require(BSGS)
data(Crisis2008BalancedData)
var.names <- colnames(Crisis2008BalancedData)[-1]
country.all <- rownames(Crisis2008BalancedData)
cov.of.interest <- colnames(Crisis2008BalancedData)[-1]
Y <- Crisis2008BalancedData[, 1]
Y <- Y - mean(Y)
X <- Crisis2008BalancedData[, -1]
dummy.variable <- cov.of.interest[lapply(apply(X, 2, unique), length) == 2]
non.dummy.X <- X[, !(colnames(X) %in% dummy.variable)]
X.normalized <- apply(non.dummy.X, 2, function(XX) (XX - mean(XX))/sd(XX))
X[, !(colnames(X) %in% dummy.variable)] <- X.normalized
n <- dim(X)[1]
p <- dim(X)[2]
group_size <- c(10, 3, 4, 2, 4, 11, 4, 1, 12)

start_time <- proc.time()
res_2bg_bsgl <- bsgl(X, Y, group_size, rep(1, dim(X)[2]), 1, lambda1 = .104, lambda2 = .082, K = 5000)
time_2bg_econ <- as.double(proc.time() - start_time)[3]
start_time <- proc.time()
res_3bg_bsgl <- bsgl_original(X, Y, group_size, rep(1, dim(X)[2]), 1, lambda1 = .104, lambda2 = .082, K = 5000)
time_3bg_econ <- as.double(proc.time() - start_time)[3]
stan_dat <- list(n = n,
                 p = p,
                 lambda1 = .104,
                 lambda2 = .082,
                 X = X,
                 Y = Y,
                 group_size = group_size,
                 n_group = 9)

fit_hmc_bsgl <- stan(file = 'bsglasso.stan', data = stan_dat, control = list(max_treedepth = 15),
                     iter = 5000, chains = 1, warmup = 500)
time_hmc_econ <- get_elapsed_time(fit_hmc_bsgl)[,2]
t <- stan_ac(fit_hmc_bsgl, pars = 'sigma')
t1 <- t$data

# Autocorrelation
t1$ac
autocorr(as.mcmc(res_2bg_bsgl$sigma2s[-c(1:500)]))
autocorr(as.mcmc(res_3bg_bsgl$sigma2s[-c(1:500)]))

# Time
time_2bg_econ * .9
time_3bg_econ * .9
time_hmc_econ

# effective sample size
effectiveSize(as.mcmc(res_2bg_bsgl$sigma2s[-c(1:500)]))
effectiveSize(as.mcmc(res_3bg_bsgl$sigma2s[-c(1:500)]))
summary(fit_hmc_bsgl)$summary['sigma','n_eff']

############################################################
# Generate MCMC output for Subsection 5.3
# (CGH Data)
############################################################
Rcpp::sourceCpp('bfl.cpp') 
library(cghFLasso)
library(plotrix)
data(CGH)
Y <- CGH$GBM.y[1:200]
n <- length(Y)
p <- n
X <- diag(p)
start_time <- proc.time()
res_2bg_bfl <- bfl(X, Y, rep(1,p)+rnorm(p), 1, lambda1 = .129, lambda2 = .962, K = 5000)
time_2bg_cgh <- as.double(proc.time() - start_time)[3]
start_time <- proc.time()
res_3bg_bfl <- bfl_original(X, Y, rep(1,p)+rnorm(p), 1, lambda1 = .129, lambda2 = .962, K = 5000)
time_3bg_cgh <- as.double(proc.time() - start_time)[3]
stan_dat <- list(n = n,
                 p = p,
                 lambda1 = .129,
                 lambda2 = .962,
                 X = X,
                 Y = Y)
fit_hmc_bfl <- stan(file = 'bfl.stan', data = stan_dat, control = list(max_treedepth = 15),
                    iter = 5000, chains = 1, warmup = 500)
time_hmc_cgh <- get_elapsed_time(fit_hmc_bfl)[,2]
t <- stan_ac(fit_hmc_bfl, pars = 'sigma')
t1 <- t$data

# Autocorrelation
t1$ac
autocorr(as.mcmc(res_2bg_bfl$sigma2s[-c(1:500)]))
autocorr(as.mcmc(res_3bg_bfl$sigma2s[-c(1:500)]))

# Time
time_2bg_cgh * .9
time_3bg_cgh * .9
time_hmc_cgh

# effective sample size
effectiveSize(as.mcmc(res_2bg_bfl$sigma2s[-c(1:500)]))
effectiveSize(as.mcmc(res_3bg_bfl$sigma2s[-c(1:500)]))
summary(fit_hmc_bfl)$summary['sigma','n_eff']
}