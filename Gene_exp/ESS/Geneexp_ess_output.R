rm(list=ls())
Rcpp::sourceCpp('../bgl.cpp') 
library(gglasso)
library("coda")
library(bpr)
library(foreach)
library(doParallel)
library(Rcpp)
library(RcppDist)
source("Covariance_gene.R")
################################Initialization
data(bardet)
group_size <- rep(5, 20)
X.raw <- bardet$x
Y.raw <- bardet$y
n <- dim(X.raw)[1]
p <- dim(X.raw)[2]
X <- scale(X.raw)*sqrt(n/(n-1))
X <- matrix(as.vector(X),n,p)
Y <- Y.raw-mean(Y.raw)

N = 5e5
nloops <- 10
B <- 10
subsize <- floor(seq(1e4, N, length = nloops))

#----------------------------------------------------------------------------------#
Table = list()

parallel::detectCores()
n.cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = 12)



# For loop for getting ESS and norm values corresponding to different chain sizes
Table = foreach(b=1:B, .packages = c("mcmcse"))%dopar%{
  #if(b %% 1 ==0) print(b)
  chain <- bgl(X, Y,group_size, rep(1,p), 1, lambda = 0.0601 ,K = N)$beta
  combine = list()
  for(j in 1:length(subsize)){
    print(j)
    minichain <- chain[1:subsize[j], ]
    ess_track <- var_track(minichain)

    ess_track_bm <- multiESS(minichain, covmat = ess_track$bm_est)/subsize[j]
    ess_track_ise <- multiESS(minichain, covmat = ess_track$ise_est)/subsize[j]
    ess_track_cc <- multiESS(minichain, covmat = ess_track$cc_est)/subsize[j]
    ess_track_sve <- multiESS(minichain, covmat = ess_track$sve_est)/subsize[j]

    ess_list = list(ess_track_bm, ess_track_ise, ess_track_cc, ess_track_sve)

    combine = append(combine, list(ess_list))
  }
  combine
}
  
save(Table, subsize, nloops, B, file = "ESS_data.Rdata")
