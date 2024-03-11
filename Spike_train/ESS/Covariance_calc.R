source("../Cov_func.R")
#library(matrixcalc)


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
var_track <- function(minichain){

	t = Sys.time()
	bm_est = mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov
	bm_time = Sys.time() - t
	#lug_est = mcse.multi(minichain, r = 3, method = "bm", adjust = FALSE)$cov
	t = Sys.time()
	ise_est = mcse.initseq(minichain)$cov
	ise_time = Sys.time() - t

	t = Sys.time()
	sve_est = mcse.multi(minichain, r = 1,method = "tukey", adjust = FALSE)$cov
	sve_time = Sys.time() - t

	t = Sys.time()
	cc_est = cov.sig(minichain, type = "geyer")$covariance
	cc_time = Sys.time() - t
	#mls_est = cov.sig(minichain, type = "MomentLS")$covariance
	return(list("bm_est" = bm_est, "ise_est" = ise_est, "sve_est" = sve_est, "cc_est" = cc_est,
		"bm_time" = bm_time, "ise_time" = ise_time, "sve_time" = sve_time,
		"cc_time" = cc_time))
}


