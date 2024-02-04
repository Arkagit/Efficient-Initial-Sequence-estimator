set.seed(111)
library(mcmcse)
library("coda")
library(bpr)
library(foreach)
library(doParallel)
library(Rcpp)
library(RcppDist)
load("..//data_calcium_imaging_for_poisson.RData")
str(data)

source("../Cov_func.R")
#library(matrixcalc)
data$depth2 = data$depth^2
X = model.matrix(~ ., data = data[,c(2:4,6,9)])
str(X)
p = ncol(X)
n = nrow(X)
y = data$n_spikes
str(y)

mle.start <- c(summary(glm(y ~ X - 1, family = "poisson"(link = "log")))$coef[,1])
str(mle.start)

N = seq(1e5, 1e6, length.out = 10)
repet = 1e2
burnin = 1:100
nloops <- 100
#subsize <- floor(seq(1e4, N, length = nloops))
data = data.frame(y=y, X)

#----------------------------------------------------------------------------------#



bm_time = rep(0, length(N))
lug_time = rep(0, length(N))
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
	 chain <- bpr::sample_bpr(y ~ . - 1, data = data,
                      iter = max(N), burnin = max(burnin),
                      prior = list(type="gaussian", b = rep(0,p), B = diag(p)*2), 
                      pars = list(max_dist = 1e+6),
                      state = mle.start)$sim$beta
	for(j in 1:length(N)){
		minichain = chain[1:N[j],]

				
		bm_time[j] = system.time(mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov)[3]
		
		lug_time[j] = system.time(mcse.multi(minichain, r = 3, method = "bm", adjust = FALSE)$cov)[3]
		
		ise_time[j] = system.time(mcse.initseq(minichain)$cov)[3]
		
		sve_time[j] = system.time(mcse.multi(minichain, r = 1,method = "tukey", adjust = FALSE)$cov)[3]
		
		cc_time[j] = system.time(cov.sig(minichain, type = "geyer")$covariance)[3]
		
		mls_time[j] = system.time(cov.sig(minichain, type = "MomentLS")$covariance)[3]
	}

	comb = list(bm_time, lug_time, ise_time, sve_time, cc_time, mls_time)

	comb
}

save(Table, N, repet, file = "Time_d.Rdata")