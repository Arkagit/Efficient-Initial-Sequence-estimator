set.seed(111)
library(mcmcse)
library(foreach)
library(doParallel)

source("../Asymp_var.R")
source("../VAR_func.R")
source("../Cov_func.R")
#library(matrixcalc)
N = seq(1e5, 1e6, length.out = 10)
repet = 1e2

bm_time = rep(0, length(N))
ise_time = rep(0, length(N))
sve_time = rep(0, length(N))
cc_time = rep(0, length(N))
mls_time = rep(0, length(N))

Table = list()

parallel::detectCores()
n.cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n.cores)



Table = foreach(i = 1:repet, .packages = c("mcmcse"))%dopar%{
	print(i)
  	chain <- var1(p = p, phi = phi, nsim = max(N), omega = omega)
	combine = list()

	for(j in 1:length(N)){
		minichain = chain[1:N[j],]
			
		bm_time[j] = system.time(mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov)[3]
		
		ise_time[j] = system.time(mcse.initseq(minichain)$cov)[3]
		
		sve_time[j] = system.time(mcse.multi(minichain, r = 1,method = "tukey", adjust = FALSE)$cov)[3]
		
		cc_time[j] = system.time(cov.sig(minichain, type = "geyer")$covariance)[3]

		mls_time[j] = system.time(cov.sig(minichain, type = "MomentLS")$covariance)[3]
		
	}

	comb = list(bm_time, ise_time, sve_time, cc_time, mls_time)

	comb
}

save(Table, N, repet, file = "Time_d.Rdata")