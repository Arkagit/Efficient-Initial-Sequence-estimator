library(phonTools)
library(mcmcse)
library(multichainACF)

parallel_gacf = function(par_chain, M, ch_l, chain_val){
	
	sigma2 = numeric(length = p)

	for(i in 1:p){
		GACF = globalACF(par_chain, chains = chain_val, component = i, type = "covariance", lag.max = ch_l - 1, plot = FALSE)
		sigma2[i] = - GACF$chainsACF[[chain_val]][[1]][1]

		for(k in 1:floor(ch_l/2 - 1)){
			dummy1 = GACF$chainsACF[[chain_val]][[1]][2*k - 1]
			dummy2 = GACF$chainsACF[[chain_val]][[1]][2*k]
			if(dummy1 + dummy2 > 0){
				sigma2[i] = sigma2[i] + 2*(dummy1 + dummy2)
			}else{
				break
			}
		}
	} 
	return(sigma2)
}


parallel_sig = function(par_chain, M, ch_l, chain_val){
	
	sigma = sqrt(parallel_gacf(par_chain, M, ch_l, chain_val))
	dat = par_chain[[chain_val]]
	bm = mcse.multi(dat, method = "bm", r = 1, adjust = FALSE)$cov
	est = diag(sigma)%*%cov2cor(bm)%*%diag(sigma)
	return(est)
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
