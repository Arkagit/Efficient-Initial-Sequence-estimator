source(paste(dirname(getwd()),"/VAR_func.R", sep = ""))
source(paste(dirname(getwd()),"/Cov_func.R", sep = ""))
library(matrixcalc)
set.seed(101)


#Ess.cov.sig <- function(data){
#	chain <- as.matrix(data)
 #   if(det(cov.sig(chain)$covariance) > 0){
  #  	ess =  nrow(chain)*exp((log(det(var(chain))) - log(det(cov.sig(chain)$covariance)))/ncol(chain))
   # }else{
    #	stop("covariance matrix is not positive definite")
    #}
    #return("ESS" = ess)
#}


var_track <- function(N = 5e4, phi, omega, B = 10, level = .90)
{
	p <- dim(phi)[1]
	nloops <- 50
	subsize <- floor(seq(5e3, N, length = nloops))
	ess_track_lug <- matrix(0, nrow = B, ncol = length(subsize))
	ess_track_bm <- matrix(0, nrow = B, ncol = length(subsize))
	ess_track_cc <- matrix(0, nrow = B, ncol = length(subsize))
	ess_track_ise <- matrix(0, nrow = B, ncol = length(subsize))
	

    fro_track_lug <- matrix(0, nrow = B, ncol = length(subsize))
	fro_track_bm <- matrix(0, nrow = B, ncol = length(subsize))
	fro_track_cc <- matrix(0, nrow = B, ncol = length(subsize))
	fro_track_ise <- matrix(0, nrow = B, ncol = length(subsize))

	
	colnames(ess_track_bm) <- subsize
	colnames(ess_track_lug) <- subsize
	colnames(ess_track_cc) <- subsize
	colnames(ess_track_ise) <- subsize

	
	colnames(fro_track_bm) <- subsize
	colnames(fro_track_lug) <- subsize
	colnames(fro_track_cc) <- subsize
	colnames(fro_track_ise) <- subsize

	for(b in 1:B)
	{
		if(b %% 1 ==0) print(b)
		chain <- var1(p = p, phi = phi, nsim = N, omega = omega)

		for(j in 1:length(subsize))
		{
			print(j)
			minichain <- chain[1:subsize[j], ]
			ess_track_bm[b, j] <- multiESS(minichain, r = 1, method = "bm", adjust = FALSE)/subsize[j]
			ess_track_lug[b,j] <- multiESS(minichain, r = 3, method = "bm", adjust = FALSE)/subsize[j]
			ess_track_ise[b,j] <- multiESS(minichain, covmat = mcse.initseq(minichain)$cov)/subsize[j]
			ess_track_cc[b,j] <- multiESS(minichain, covmat = cov.sig(minichain)$covariance)/subsize[j]

			fro_track_bm[b, j] <- norm(mcse.multi(minichain, method = "bm", r = 1)$cov, type = "F")
			fro_track_lug[b,j] <- norm(mcse.multi(minichain, method = "bm", r = 3)$cov, type = "F")
			fro_track_ise[b,j] <- norm(mcse.initseq(minichain)$cov, type = "F")
			fro_track_cc[b,j] <- norm(cov.sig(minichain)$covariance, type = "F")

		}
	}
	return(list("BM" = ess_track_bm, "Frob_bm" = fro_track_bm, "LUG" = ess_track_lug, "Frob_lug" = fro_track_lug, "CC" = ess_track_cc, "Frob_cc" = fro_track_cc, "ISE" = ess_track_ise, "Frob_ise" = fro_track_ise))
}



