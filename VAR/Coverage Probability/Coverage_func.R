source(paste(dirname(getwd()),"/VAR_func.R", sep = ""))
source(paste(dirname(getwd()),"/Cov_func.R", sep = ""))
library(matrixcalc)
set.seed(101)

coverage <- function(N = 1e5, phi, omega, level = .90){
	p <- dim(phi)[1]
	truth <- true.sig(rep(1,p), p, omega, phi)$final.cov
	nloops <- 50
	subsize <- c(5e3, 1e4, 5e4)
	H <- qchisq(level, p)
	M <- matrix(0, nrow = length(subsize), ncol = 4)
	colnames(M) <- c("BM", "Lugsail", "ISE", "New ISE")
	rownames(M) <- subsize
	for(i in 1:length(subsize)){
		print(i)
		for (j in 1:1000) {
			minichain = var1(p = p, phi = phi, nsim = N, omega = omega)[1:subsize[i],]
			bm_est = mcse.multi(minichain, r = 1, method = "bm", adjust = "FALSE")$cov
			lug_est = mcse.multi(minichain, r = 3, method = "bm", adjust = "FALSE")$cov
			ise_est = mcse.initseq(data.frame(minichain))$cov
			cc_est = cov.sig(data.frame(minichain))$covariance

			if(subsize[i]*t(matrix(colMeans(minichain), nrow = 10))%*%solve(bm_est)%*%matrix(colMeans(minichain), nrow = 10) < H)
				M[i,1] = M[i,1] + 1/1000

			if(subsize[i]*t(matrix(colMeans(minichain), nrow = 10))%*%solve(lug_est)%*%matrix(colMeans(minichain), nrow = 10) < H)
				M[i,2] = M[i,2] + 1/1000

			if(subsize[i]*t(matrix(colMeans(minichain), nrow = 10))%*%solve(ise_est)%*%matrix(colMeans(minichain), nrow = 10) < H)
				M[i,3] = M[i,3] + 1/1000

			if(subsize[i]*t(matrix(colMeans(minichain), nrow = 10))%*%solve(cc_est)%*%matrix(colMeans(minichain), nrow = 10) < H)
				M[i,4] = M[i,4] + 1/1000
		}
	}

	return(M)
}