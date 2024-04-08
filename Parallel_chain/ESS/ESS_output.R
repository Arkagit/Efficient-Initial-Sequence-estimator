set.seed(111)

source("../Asymp_var.R")
source("../VAR_func.R")
source("../global_var.R")
source("../Stan_code.R")

library(foreach)
library(doParallel)


B = 100

M = 10

subsize <- c(5e3, 8e3, 1e4, 3e4, 5e4, 8e4, 1e5)

N = max(subsize)

chain_val = 1

#parallel::detectCores()
#n.cores <- 50
#doParallel::registerDoParallel(cores = n.cores)

Estimate = foreach(b = 1:B, .packages = c("mcmcse"))%dopar%{
	data = list()
	for(i in 1:M){
		data = append(data, list(var1(p = p, phi = phi, nsim = N, omega = omega)))
	}

	print(4*b+1)

	Sigma_global = list()
	Sigma_stan = list()
	for(u in 1:length(subsize)){
		minichain = list()
		for(v in 1:M){
			minichain[[v]] = data[[v]][1:subsize[u],]
		}
		print(4*b+2)

		Sigma_global[[u]] = multiESS(minichain[[1]], parallel_sig(minichain, M, subsize[u], chain_val))

		print(4*b + 3)
		Sigma_stan[[u]] = multiESS(minichain[[1]], parallel_stan_sig(minichain, M, subsize[u]))
	}
	
	list(Sigma_global, Sigma_stan)
}

save(Estimate, B, M, N, subsize, chain_val, file = "parallel_est.Rdata")