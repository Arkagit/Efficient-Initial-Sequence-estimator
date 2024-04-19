library(phonTools)
library(mcmcse)
library(multichainACF)

parallel_gacf = function(par_chain, M, ch_l){
	
	sigma2 = numeric(length = p)
	Mat  = matrix(0, nrow = ch_l, ncol = M)

	for(i1 in 1:p){

		GACF = globalACF(par_chain, component = i1, type = "covariance", lag.max = ch_l - 1, plot = FALSE)

		for(i2 in 1:M){
			Mat[,i2] = GACF$chainsACF[[i2]][[1]]
		}
		covar = rowMeans(Mat)

		for(k in 1:floor(ch_l/2 - 1)){

			dummy1 = covar[2*k - 1]
			dummy2 = covar[2*k]
			if(dummy1 + dummy2 > 0){
				sigma2[i1] = sigma2[i1] + 2*(dummy1 + dummy2)
			}else{
				break
			}
		}
	} 
	return(sigma2)
}


parallel_sig = function(par_chain, M, ch_l){
	
	sigma = sqrt(parallel_gacf(par_chain, M, ch_l))
	bm = matrix(0, nrow = p, ncol = p)
	for(i3 in 1:M){
		dat = par_chain[[i3]]
		bm = bm + mcse.multi(dat, method = "bm", r = 1, adjust = FALSE)$cov/M
	}
	est = diag(sigma)%*%cov2cor(bm)%*%diag(sigma)
	return(est)
}

Gamma0 = function(par_chain, M){
	var0 = matrix(0, nrow = p, ncol = p)
	for(i4 in 1:M){
		var0 = var0 + var(par_chain[[i4]])/M
	}
	return(var0)
}


#M = 10

#N = 1e4

#p = 2

#chain_val = 1

#par_chain = list()

#for(i in 1:M){
#	par_chain = append(par_chain, list(matrix(c(rnorm(N), rgamma(N, 2)), nrow = N, byrow = FALSE)))
#}

#parallel_gacf(par_chain, M, N, 1)

#parallel_sig(par_chain, M, N, 1)
