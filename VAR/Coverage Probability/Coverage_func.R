source("../VAR_func.R")
source("../Cov_func.R")
library(matrixcalc)


coverage <- function(subsize, phi, omega, level = .90){
	p <- dim(phi)[1]
	nloops <- 50
	qchisq_qnt <- qchisq(level, p)
	count_mat <- matrix(0, nrow = length(subsize), ncol = 4)
	colnames(count_mat) <- c("BM", "Lugsail", "ISE", "New ISE")
	rownames(count_mat) <- subsize
	chain = var1(p = p, phi = phi, nsim = max(subsize), omega = omega)
	true_mean = rep(0, p)

	for (i in 1:length(subsize)) {
    	minichain <- chain[1:subsize[i],]
		bm_est = mcse.multi(minichain, r = 1, method = "bm", adjust = FALSE)$cov
		lug_est = mcse.multi(minichain, r = 3, method = "bm", adjust = FALSE)$cov
		ise_est = mcse.initseq(data.frame(minichain))$cov
		cc_est = cov.sig(data.frame(minichain))$covariance
		chain_mean = matrix(colMeans(minichain), nrow = p)


		if(subsize[i]*t(chain_mean - true_mean)%*%solve(bm_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,1] = count_mat[i,1] + 1

		if(subsize[i]*t(chain_mean - true_mean)%*%solve(lug_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,2] = count_mat[i,2] + 1

		if(subsize[i]*t(chain_mean - true_mean)%*%solve(ise_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,3] = count_mat[i,3] + 1

		if(subsize[i]*t(chain_mean - true_mean)%*%solve(cc_est)%*%(chain_mean - true_mean) < qchisq_qnt)
			count_mat[i,4] = count_mat[i,4] + 1
	}

	return(count_mat)
}