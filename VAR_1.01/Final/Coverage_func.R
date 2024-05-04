source("../VAR_func.R")
source("../Cov_func.R")
source("../Asymp_var.R")
#library(matrixcalc)
sourceCpp("inseq.cpp")


coverage <- function(subsize, phi, omega, level = .90){
	p <- dim(phi)[1]
	nloops <- 50
	qchisq_qnt <- qchisq(level, p)
	count_mat <- matrix(0, nrow = length(subsize), ncol = 5)
	colnames(count_mat) <- c("BM", "ISE", "New ISE (Geyer)", "SVE", "New ISE (MLS)")
	rownames(count_mat) <- subsize
	true_mean = rep(0, p)
	time = matrix(0, nrow = length(subsize), ncol = 5)
	est = list()
	ess = list()
	trunc_ise = matrix(0, nrow = length(subsize), ncol = 2)
	chain = var1(p = p, phi = phi, nsim = max(subsize), omega = omega)

	for (i in 1:length(subsize)) {

    	minichain <- matrix(chain[1:subsize[i],], nrow = subsize[i])
    	
		time[i,1] = system.time(bm_est <-  mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE))[3]
		
		time[i,2] = system.time(ise_est <- inseq(minichain))[3]
		
		time[i,3] = system.time(cc_est <- cov.sig(minichain, type = "geyer"))[3]
		
		time[i,4] = system.time(sve_est <- mcse.multi(minichain, r = 1, method = "tukey", adjust = FALSE))[3]
		
		time[i,5] = system.time(mls_est <- cov.sig(minichain, type = "MomentLS"))[3]

		trunc_ise[i,] = c(ise_est$trunc, max(cc_est$stopping_lag))

		bm_est = bm_est$cov
		ise_est = ise_est$Sig
		cc_est = cc_est$covariance
		sve_est = sve_est$cov
		mls_est = mls_est$covariance

		chain_mean = matrix(colMeans(minichain), nrow = p)


		if(subsize[i]*t(chain_mean - true_mean)%*%solve(bm_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,1] = count_mat[i,1] + 1

		if(subsize[i]*t(chain_mean - true_mean)%*%solve(ise_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,2] = count_mat[i,2] + 1

		if(subsize[i]*t(chain_mean - true_mean)%*%solve(cc_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,3] = count_mat[i,3] + 1

		if(subsize[i]*t(chain_mean - true_mean)%*%solve(sve_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,4] = count_mat[i,4] + 1

		if(subsize[i]*t(chain_mean - true_mean)%*%solve(mls_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,5] = count_mat[i,5] + 1

		bm_ess = multiESS(minichain, covmat = bm_est)

		ise_ess = multiESS(minichain, covmat = ise_est)

		cc_ess = multiESS(minichain, covmat = cc_est)

		sve_ess = multiESS(minichain, covmat = sve_est)

		mls_ess = multiESS(minichain, covmat = mls_est)

		est = append(est, list(bm_est, ise_est, cc_est, sve_est, mls_est))

		ess = append(ess, list(bm_ess, ise_ess, cc_ess, sve_ess, mls_ess))
	}

	return(list("count_mat" = count_mat, "Time" = time, "estimates" = est, "ESS" = ess, "Truncation" = trunc_ise))
}


