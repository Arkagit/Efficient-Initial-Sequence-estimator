rm(list=ls())
library("coda")
library(bpr)
library(foreach)
library(doParallel)
library(Rcpp)
library(RcppDist)
load("..//data_calcium_imaging_for_poisson.RData")
str(data)
#install.packages("bpr")
#install.packages("RcppDist")

library(mcmcse)
source("Covariance_calc.R")
#### quadratic on depth
data$depth2 = data$depth^2
X = model.matrix(~ ., data = data[,c(2:4,6,9)])
str(X)
p = ncol(X)
n = nrow(X)
y = data$n_spikes
str(y)

mle.start <- c(summary(glm(y ~ X - 1, family = "poisson"(link = "log")))$coef[,1])
str(mle.start)


N = 5e5
burnin = 1:100
nloops <- 100
B <- 10
subsize <- floor(seq(1e4, N, length = nloops))
data = data.frame(y=y, X)

#----------------------------------------------------------------------------------#
Table = list()

parallel::detectCores()
n.cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n.cores)



# For loop for getting ESS and norm values corresponding to different chain sizes
Table = foreach(b=1:B, .packages = c("mcmcse"))%dopar%{
  #if(b %% 1 ==0) print(b)
  chain <- bpr::sample_bpr(y ~ . - 1, data = data,
                      iter = N, burnin = max(burnin),
                      prior = list(type="gaussian", b = rep(0,p), B = diag(p)*2), 
                      pars = list(max_dist = 1e+6),
                      state = mle.start)$sim$beta
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
  
save(Table, subsize, nloops, B, file = "ESS_Frob_data.Rdata")
