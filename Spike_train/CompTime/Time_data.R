set.seed(100)
library(mcmcse)
library(foreach)
library(doParallel)

source("../Asymp_var.R")
source("../VAR_func.R")
source("../Cov_func.R")
#library(matrixcalc)
N = 1e6
repet = 1e4
#Ess.cov.sig <- function(data){
#	chain <- as.matrix(data)
 #   if(det(cov.sig(chain)$covariance) > 0){
  #  	ess =  nrow(chain)*exp((log(det(var(chain))) - log(det(cov.sig(chain)$covariance)))/ncol(chain))
   # }else{
    #	stop("covariance matrix is not positive definite")
    #}
    #return("ESS" = ess)
#}

# function for estimates of any given chain
bm_time = 0
lug_time = 0
ise_time = 0
sve_time = 0
cc_time = 0
mls_time = 0

Table = list()

parallel::detectCores()
n.cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n.cores)


Table = foreach(i = 1:repet, .packages = c("mcmcse"))%dopar%{
	print(i)
	minichain = var1(p = p, phi = phi, nsim = N, omega = omega)
	t = Sys.time()
	bm_est = mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov
	bm_time = Sys.time() - t
	t = Sys.time()
	lug_est = mcse.multi(minichain, r = 3, method = "bm", adjust = FALSE)$cov
	lug_time = Sys.time() - t
	t = Sys.time()
	ise_est = mcse.initseq(minichain)$cov
	ise_time = Sys.time() - t
	t = Sys.time()
	sve_est = mcse.multi(minichain, r = 1,method = "tukey", adjust = FALSE)$cov
	sve_time = Sys.time() - t
	t = Sys.time()
	cc_est = cov.sig(minichain, type = "geyer")$covariance
	cc_time = Sys.time() - t
	t = Sys.time()
	mls_est = cov.sig(minichain, type = "MomentLS")$covariance
	mls_time = Sys.time() - t

	comb = list(bm_time, lug_time, ise_time, sve_time, cc_time, mls_time)

	comb
}




save(Table, file = "Time.Rdata")
