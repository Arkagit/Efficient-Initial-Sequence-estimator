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

N = c(5e3, 8e3, 1e4, 3e4, 5e4, 8e4, 1e5, 3e5, 5e5)
repet = 1e2
burnin = 1:100
nloops <- 100
#subsize <- floor(seq(1e4, N, length = nloops))
data = data.frame(y=y, X)

#----------------------------------------------------------------------------------#



bm_time = rep(0, length(N))
ise_time = rep(0, length(N))
sve_time = rep(0, length(N))
cc_time = rep(0, length(N))
mls_time = rep(0, length(N))




Table = list()

parallel::detectCores()
n.cores <- 50
doParallel::registerDoParallel(cores = n.cores)



Table = foreach(i = 1:repet, .packages = c("mcmcse"))%dopar%{
	
	 chain <- bpr::sample_bpr(y ~ . - 1, data = data,
                      iter = max(N), burnin = max(burnin),
                      prior = list(type="gaussian", b = rep(0,p), B = diag(p)*2), 
                      pars = list(max_dist = 1e+6),
                      state = mle.start)$sim$beta
	 mat = list()

	for(j in 1:length(N)){
		minichain = chain[1:N[j],]

		gam <- var(minichain)

		bm_time[j] = system.time(bm <- mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov)[3]

		ise_time[j] = system.time(ise <- mcse.initseq(minichain)$cov)[3]
		
		sve_time[j] = system.time(sve <- mcse.multi(minichain, r = 1,method = "tukey", adjust = FALSE)$cov)[3]
		
		cc_time[j] = system.time(cc <- cov.sig(minichain, type = "geyer")$covariance)[3]
		
		mls_time[j] = system.time(mls <- cov.sig(minichain, type = "MomentLS")$covariance)[3]

		mat = append(mat, list(gam, bm, ise, sve, cc, mls))
	}

	print(i)

	comb_time = list(bm_time, ise_time, sve_time, cc_time, mls_time)

	list(mat, comb_time)
}

save(Table, N, repet, file = "Time_d.Rdata")