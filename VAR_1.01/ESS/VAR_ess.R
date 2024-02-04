######### Sourcing new covariance files
source("../Cov_func.R")

# function for estimates of any given chain
var_track <- function(minichain){
	bm_est = mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov
	lug_est = mcse.multi(minichain, r = 3, method = "bm", adjust = FALSE)$cov
	ise_est = mcse.initseq(minichain)$cov
	sve_est = mcse.multi(minichain, r = 1,method = "tukey", adjust = FALSE)$cov
	cc_est = cov.sig(minichain, type = "geyer")$covariance
	mls_est = cov.sig(minichain, type = "MomentLS")$covariance
	return(list("bm_est" = bm_est, "lug_est" = lug_est, "ise_est" = ise_est,
		"sve_est" = sve_est, "cc_est" = cc_est, "mls_est" = mls_est))
}


