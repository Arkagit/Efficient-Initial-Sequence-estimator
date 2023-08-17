set.seed(100)
source("VAR_ess.R")

# Initializing variable values
N <- 5e4
nloops <- 50
subsize <- floor(seq(5e3, N, length = nloops))
B <- 10
p <- 10
rho <- 0.99
phi <- diag(rep(rho,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))
truth <- true.sig(p = p, omega = omega, phi = phi)
ess_true <- (det(truth$tar.var)/det(truth$final.cov))^(1/p)


# Allocating memory for ess and f_norm values
ess_track_lug <- matrix(0, nrow = B, ncol = length(subsize))
ess_track_bm <- matrix(0, nrow = B, ncol = length(subsize))
ess_track_cc <- matrix(0, nrow = B, ncol = length(subsize))
ess_track_ise <- matrix(0, nrow = B, ncol = length(subsize))
ess_track_sve <- matrix(0, nrow = B, ncol = length(subsize))
ess_track_mls <- matrix(0, nrow = B, ncol = length(subsize))


fro_track_lug <- matrix(0, nrow = B, ncol = length(subsize))
fro_track_bm <- matrix(0, nrow = B, ncol = length(subsize))
fro_track_cc <- matrix(0, nrow = B, ncol = length(subsize))
fro_track_ise <- matrix(0, nrow = B, ncol = length(subsize))
fro_track_sve <- matrix(0, nrow = B, ncol = length(subsize))
fro_track_mls <- matrix(0, nrow = B, ncol = length(subsize))

# Naming columns of final table
colnames(ess_track_bm) <- subsize
colnames(ess_track_lug) <- subsize
colnames(ess_track_cc) <- subsize
colnames(ess_track_ise) <- subsize
colnames(ess_track_sve) <- subsize
colnames(ess_track_mls) <- subsize

colnames(fro_track_bm) <- subsize
colnames(fro_track_lug) <- subsize
colnames(fro_track_cc) <- subsize
colnames(fro_track_ise) <- subsize
colnames(fro_track_sve) <- subsize
colnames(fro_track_mls) <- subsize


# For loop for getting ESS and norm values corresponding to different chain sizes
for(b in 1:B){
	if(b %% 1 ==0) print(b)
	chain <- var1(p = p, phi = phi, nsim = max(subsize), omega = omega)
	for(j in 1:length(subsize)){
		print(j)
		minichain <- chain[1:subsize[j], ]
		ess_track <- var_track(minichain)

		ess_track_bm[b, j] <- multiESS(minichain, covmat = ess_track$bm_est)/subsize[j]
		ess_track_lug[b,j] <- multiESS(minichain, covmat = ess_track$lug_est)/subsize[j]
		ess_track_ise[b,j] <- multiESS(minichain, covmat = ess_track$ise_est)/subsize[j]
		ess_track_cc[b,j] <- multiESS(minichain, covmat = ess_track$cc_est)/subsize[j]
		ess_track_sve[b,j] <- multiESS(minichain, covmat = ess_track$sve_est)/subsize[j]
		ess_track_mls[b,j] <- multiESS(minichain, covmat = ess_track$mls_est)/subsize[j]


		fro_track_bm[b, j] <- norm(ess_track$bm_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_lug[b,j] <- norm(ess_track$lug_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_ise[b,j] <- norm(ess_track$ise_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_cc[b,j] <- norm(ess_track$cc_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_sve[b,j] <- norm(ess_track$sve_est, type = "F")/norm(truth$final.cov, type = "F")
		fro_track_mls[b,j] <- norm(ess_track$mls_est, type = "F")/norm(truth$final.cov, type = "F")
	}	
}

# plot(subsize, colMeans(ess_track$BM), type= 'l', ylim = range(ess_track))
# lines(subsize, colMeans(ess_track$WBM), col = "red")
# lines(subsize, colMeans(ess_track$LUG), col = "blue")
# abline(h = ess_true)
save(ess_track_bm, fro_track_bm, ess_track_lug, fro_track_lug,
     ess_track_cc, fro_track_cc, ess_track_ise, fro_track_ise, ess_track_sve, fro_track_sve,
       ess_track_mls, fro_track_mls, subsize, ess_true, nloops, file = "ESS_data_0.99.Rdata")
