library(phonTools)
library(mcmcse)

parallel_gacf = function(par_chain, M, N, chain_val){
	chain_mean = matrix(0, nrow = M, ncol = p)
	for(i in 1:M){
		chain_mean[i,] = colMeans(par_chain[[i]])
	}
	gl_mean = colMeans(chain_mean)
	m_diff = gl_mean - chain_mean[chain_val,]
	sigma2 = numeric(length = p)

	for(i in 1:p){

		acf_c = var(par_chain[[chain_val]][,i])*fastacf(par_chain[[chain_val]][,i], show = "FALSE", window = "rectangular", lag.max = floor(N/2))$acf
		sigma2[i] = - acf_c[1]

		for(k in 1:floor(N/2)){
			dummy1 = acf_c[2*k - 1] + (N - 2*k + 1)*(gl_mean[i]^(2) - chain_mean[chain_val, i]^(2))/N - 
			         m_diff[i]*(sum(par_chain[[chain_val]][1:(N-2*k + 1),i]) + sum(par_chain[[chain_val]][(2*k):N,i]))/N

	        dummy2 = acf_c[2*k] + (N - 2*k)*(gl_mean[i]^(2) - chain_mean[chain_val, i]^(2))/N - 
			         m_diff[i]*(sum(par_chain[[chain_val]][1:(N-2*k),i]) + sum(par_chain[[chain_val]][(2*k+1):N,i]))/N
			
			if(dummy1 + dummy2 > 0){
				sigma2[i] = sigma2[i] + 2*(dummy1 + dummy2)
			}else{
				break
			}
		}
	} 
	return(sigma2)
}


parallel_sig = function(par_chain, M, N, chain_val){
	
	sigma = sqrt(parallel_gacf(par_chain, M, N, chain_val))
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
