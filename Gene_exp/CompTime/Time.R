set.seed(100)
library(mcmcse)
library(foreach)
library(doParallel)
library(gglasso)
library(Rcpp)
library(RcppDist)

Rcpp::sourceCpp('../bgl.cpp') 
N = seq(1e5, 1e6, length.out = 10)
repet = 50
source("../Cov_func.R")
#library(matrixcalc)
data(bardet)
group_size <- rep(5, 20)
X.raw <- bardet$x
Y.raw <- bardet$y
n <- dim(X.raw)[1]
p <- dim(X.raw)[2]
X <- scale(X.raw)*sqrt(n/(n-1))
X <- matrix(as.vector(X),n,p)
Y <- Y.raw-mean(Y.raw)

#----------------------------------------------------------------------------------#



bm_time = rep(0, length(N))
ise_time = rep(0, length(N))
sve_time = rep(0, length(N))
cc_time = rep(0, length(N))

Table = list()

parallel::detectCores()
n.cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = 12)



Table = foreach(i = 1:repet, .packages = c("mcmcse"))%dopar%{
	print(i)
  	chain <- bgl(X, Y,group_size, rep(1,p), 1, lambda = 0.0601 ,K = max(N))$beta
	combine = list()

	for(j in 1:length(N)){
		minichain = chain[1:N[j],]
			
		bm_time[j] = system.time(mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov)[3]
		
		ise_time[j] = system.time(mcse.initseq(minichain)$cov)[3]
		
		sve_time[j] = system.time(mcse.multi(minichain, r = 1,method = "tukey", adjust = FALSE)$cov)[3]
		
		cc_time[j] = system.time(cov.sig(minichain, type = "geyer")$covariance)[3]
		
	}

	comb = list(bm_time, ise_time, sve_time, cc_time)

	comb
}

save(Table, N, repet, file = "Time_d.Rdata")