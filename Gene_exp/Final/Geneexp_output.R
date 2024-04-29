rm(list=ls())
Rcpp::sourceCpp('../bfl.cpp') 
library(cghFLasso)
library(plotrix)
library(foreach)
library(doParallel)
library(Rcpp)
library(RcppDist)
source("Covariance_gene.R")
################################Initialization
data(CGH)
Y <- CGH$GBM.y[1:200]
n <- length(Y)
q <- n
X <- diag(q)

N = 5e5
nloops <- 10
B <- 100
subsize <- floor(seq(1e4, N, length = nloops))

#----------------------------------------------------------------------------------#
Table = list()

parallel::detectCores()
n.cores <- 50
doParallel::registerDoParallel(cores = n.cores)


Table = list()
# For loop for getting ESS and norm values corresponding to different chain sizes
for(b in 1:B){
#Table = foreach(b=1:B, .packages = c("mcmcse"))%dopar%{
  #if(b %% 1 ==0) print(b)
  dat <-  bfl(X, Y, rep(1,q)+rnorm(q), 1, lambda1 = .129, lambda2 = .962, K = N)
  chain <- cbind(dat$beta, dat$sigma2)
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

    time_list = list(ess_track$bm_time, ess_track$ise_time, ess_track$cc_time, ess_track$sve_time)

    combine = append(combine, list(ess_list, time_list))
  }
  Table[[b]] = combine
}
  
save(Table, subsize, nloops, B, file = "final_data.Rdata")
