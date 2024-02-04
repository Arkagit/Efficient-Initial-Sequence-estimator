set.seed(100)

#### Calling libraries
library(mcmcse)
library(foreach)
library(doParallel)

#### Sourcing function and data files
source("VAR_ess.R")
source("../VAR_func.R")
source("../Asymp_var.R")

# Initializing variable values
N <- 5e5
nloops <- 50
subsize <- floor(seq(1e4, N, length = nloops))
B <- 10
truth <- true.sig.gen(p = p, omega = omega, phi = phi)
ess_true <- (det(truth$tar.var)/det(truth$final.cov))^(1/p)


########### Initializing Output list
Table = list()

########### Setting up parallel programming config
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n.cores)



# For loop for getting ESS and norm values corresponding to different chain sizes
Table = foreach(b=1:B, .packages = c("mcmcse"))%dopar%{
	#if(b %% 1 ==0) print(b)
	chain <- var1(p = p, phi = phi, nsim = max(subsize), omega = omega)
	combine = list()
	for(j in 1:length(subsize)){
		print(j)
		minichain <- chain[1:subsize[j], ]
		ess_track <- var_track(minichain)

		ess_track_bm <- multiESS(minichain, covmat = ess_track$bm_est)/subsize[j]
		#ess_track_lug <- multiESS(minichain, covmat = ess_track$lug_est)/subsize[j]
		ess_track_ise <- multiESS(minichain, covmat = ess_track$ise_est)/subsize[j]
		ess_track_cc <- multiESS(minichain, covmat = ess_track$cc_est)/subsize[j]
		ess_track_sve <- multiESS(minichain, covmat = ess_track$sve_est)/subsize[j]
		ess_track_mls <- multiESS(minichain, covmat = ess_track$mls_est)/subsize[j]

		ess_list = list(ess_track_bm, ess_track_ise, ess_track_cc,
			ess_track_sve, ess_track_mls)


		fro_track_bm <- norm(ess_track$bm_est, type = "F")/norm(truth$final.cov, type = "F")
		#fro_track_lug <- norm(ess_track$lug_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_ise <- norm(ess_track$ise_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_cc <- norm(ess_track$cc_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_sve <- norm(ess_track$sve_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_mls <- norm(ess_track$mls_est, type = "F")/norm(truth$final.cov, type = "F")

		frob_list = list(fro_track_bm, fro_track_ise, fro_track_cc,
			fro_track_sve, fro_track_mls)

		m_list = list(ess_list, frob_list)

		combine = append(combine, list(m_list))
	}	

	combine
	#ess_list
}

############ Saving The List
save(Table, subsize, ess_true, truth, nloops, B, file = "ESS_data_1.01.Rdata")
