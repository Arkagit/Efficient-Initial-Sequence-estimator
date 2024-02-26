source("../VAR_func.R")
source("../Cov_func.R")
source("../Asymp_var.R")
#library(matrixcalc)


coverage <- function(subsize, phi, omega, level = .90){
	p <- dim(phi)[1]
	nloops <- 50
	qchisq_qnt <- qchisq(level, p)
	count_mat <- matrix(0, nrow = length(subsize), ncol = 5)
	colnames(count_mat) <- c("BM", "ISE", "New ISE (Geyer)", "SVE", "New ISE (MLS)")
	rownames(count_mat) <- subsize
	chain = var1(p = p, phi = phi, nsim = max(subsize), omega = omega)
	true_mean = rep(0, p)
	time = matrix(0, nrow = length(subsize), ncol = 5)
	est = list()
	ess = list()

	for (i in 1:length(subsize)) {

    	minichain <- data.frame(chain[1:subsize[i],])
    	
    	t = Sys.time()
		bm_est = mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov
		time[i,1] = Sys.time() - t
		
		t = Sys.time()
		ise_est = mcse.initseq(minichain)$cov
		time[i,2] = Sys.time() - t
		
		t = Sys.time()
		cc_est = cov.sig(minichain, type = "geyer")$covariance
		time[i,3] = Sys.time() - t
		
		t = Sys.time()
		sve_est = mcse.multi(minichain, r = 1, method = "tukey", adjust = FALSE)$cov
		time[i,4] = Sys.time() - t
		
		t = Sys.time()
		mls_est = cov.sig(minichain, type = "MomentLS")$covariance
		time[i,5] = Sys.time() - t

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

	return(list("count_mat" = count_mat, "Time" = time, "estimates" = est, "ESS" = ess))
}


