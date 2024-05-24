parallel_stan = function(par_chain, M, N){
	sigma2 = numeric(length = p)

	for(k in 1:p){
		X = matrix(0, nrow = N, ncol = M)
		rho = matrix(0, nrow = N, ncol = M)

		final_rho = numeric(length = N)
		for(l in 1:M){
			X[,l] = par_chain[[l]][,k]
		}

		# Mean of chains
		mean_x = colMeans(X)

		# Combined Mean
		comb_mean = mean(mean_x)

		# Within sample variances
		sm2 = numeric(length = M)

		for(i in 1:M){
			sm2[i] = sum((X[,i] - mean_x[i])^(2))/(N-1)
		}

		# Between sample variance
		B = sum((mean_x - comb_mean)^(2))*N/(M-1)

		W = mean(sm2)

		var_plus = (N-1)*W/N + B/N

		for(i in 1:M){
			rho[,i] = acf(X[,i], lag.max = (N-1), plot = FALSE)$acf
			#print(i)
		}

		# Final auto-corr estimates
		final_rho = 1 - (W - rho%*%sm2/M)/var_plus

		# Initial positive sequence estimator from final_rho
		Rho = numeric(length = floor(N/2))

		for(i in 1:floor(N/2)){
			Rho[i] = final_rho[2*i-1] + final_rho[2*i]
		}

		tau = - final_rho[1]

		for(i in 1:floor(N/2)){
			Trunc_time = i
			if(Rho[i] > 0){
				tau = tau + 2*Rho[i]
			}else{
				break
			}
		}
		sigma2[k] = tau*var_plus
		}

	return(sigma2)
}



parallel_stan_sig = function(par_chain, M, N){
	
	sigma = sqrt(parallel_stan(par_chain, M, N))
	bm = matrix(0, nrow = p, ncol = p)
	for(i3 in 1:M){
		dat = par_chain[[i3]]
		bm = bm + mcse.multi(dat, method = "bm", r = 1, adjust = FALSE)$cov/M
	}
	est = diag(sigma)%*%cov2cor(bm)%*%diag(sigma)
	return(est)
}
















