set.seed(112)

source("../Asymp_var.R")
source("../VAR_func.R")
source("../global_var.R")
source("../Stan_code.R")

library(foreach)
library(doParallel)


B = 100

M = c(2, 4, 8, 16)

subsize <- c(1e3, 5e3, 1e4, 5e4, 1e5)

N = max(subsize)

#chain_val = 1

true_var = true.sig.gen(p, omega, phi)$final.cov

parallel::detectCores()
n.cores <- 50
doParallel::registerDoParallel(cores = n.cores)

Estimate = list()

Estimate = foreach(b = 1:B, .packages = c("mcmcse"))%dopar%{
#for(b in 1:B){
	data = list()
	for(i in 1:max(M)){
		data = append(data, list(var1(p = p, phi = phi, nsim = N, omega = omega)))
	}

	result = list()

	print(4*b+1)

	for(s in 1:length(subsize)){
		result[[s]] = list()
		Sigma_global = list()
		Sigma_stan = list()
		Gamma = list()

		for(u in 1:length(M)){
			minichain = list()
			for(v in 1:M[u]){
				minichain[[v]] = data[[v]][1:subsize[s],]
			}
			print(4*b+2)

			Sigma_global[[u]] = parallel_sig(par_chain = minichain, M = M[u], ch_l = subsize[s])
			#Sigma_global[[u]] = norm(true_var - parallel_sig(par_chain = minichain, M = M[u], ch_l = subsize[s], chain_val = chain_val), type = "F")
			print(4*b + 3)
			Sigma_stan[[u]] = parallel_stan_sig(minichain, M[u], subsize[s])
			#Sigma_stan[[u]] = norm(true_var - parallel_stan_sig(minichain, M[u], subsize[s]), type = "F")
			Gamma[[u]] = Gamma0(minichain, M[u])
		}
		result[[s]] = list(Sigma_global, Sigma_stan, Gamma)
	}
	result
}

save(Estimate, B, M, N, subsize, file = "parallel_est.Rdata")