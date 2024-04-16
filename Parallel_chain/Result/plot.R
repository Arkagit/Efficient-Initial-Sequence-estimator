set.seed(100)

source("../Asymp_var.R")
source("../VAR_func.R")
source("../global_var.R")
source("../Stan_code.R")


load("parallel_est.Rdata")

############################################
################## ESS #####################


ess_true <- (det(true.sig.gen(p = p, omega = omega, phi = phi)$tar.var)/det(true.sig.gen(p = p, omega = omega, phi = phi)$final.cov))^(1/p)
ess_true


ess_track_glob = array(0, dim = c(B, length(subsize), length(M)))
ess_track_stan = array(0, dim = c(B, length(subsize), length(M)))

for(i in 1:B){
  for(j in 1:length(subsize)){
    for(k in 1:length(M)){
      ess_track_glob[i,j,k] = (det(Estimate[[i]][[j]][[3]][[k]])/det(Estimate[[i]][[j]][[1]][[k]]))^(1/p)
      ess_track_stan[i,j,k] = (det(Estimate[[i]][[j]][[3]][[k]])/det(Estimate[[i]][[j]][[2]][[k]]))^(1/p)
    }
  }
}



mean_ess_glob = numeric(0)
mean_ess_stan = numeric(0)

sd_ess_glob = numeric(0)
sd_ess_stan = numeric(0)

for(i in 1:length(M)){
	mean_ess_glob[i] = mean(ess_track_glob[,7,i])
	sd_ess_glob[i] = sd(ess_track_glob[,7,i])
	mean_ess_stan[i] = mean(ess_track_stan[,7,i])
	sd_ess_stan[i] = sd(ess_track_stan[,7,i])

}


final_mat_ess = matrix(c(mean_ess_glob, mean_ess_stan), byrow = TRUE, nrow = 2)
rownames(final_mat_ess) = c("Globally Centered", "STAN")
colnames(final_mat_ess) = c("M = 2", "M = 4", "M = 8", "M = 16")
final_mat_ess

sd_mat_ess = matrix(c(sd_ess_glob, sd_ess_stan), byrow = TRUE, nrow = 2)
rownames(sd_mat_ess) = c("Globally Centered", "STAN")
colnames(sd_mat_ess) = c("M = 2", "M = 4", "M = 8", "M = 16")
sd_mat_ess


############################################
################## Frob #####################

Truth = true.sig.gen(p = p, omega = omega, phi = phi)$final.cov

frob_track_glob = array(0, dim = c(B, length(subsize), length(M)))
frob_track_stan = array(0, dim = c(B, length(subsize), length(M)))

for(i in 1:B){
  for(j in 1:length(subsize)){
    for(k in 1:length(M)){
      frob_track_glob[i,j,k] = norm(Truth - Estimate[[i]][[j]][[1]][[k]], type = "F")
      frob_track_stan[i,j,k] = norm(Truth - Estimate[[i]][[j]][[2]][[k]], type = "F")
    }
  }
}



mean_frob_glob = numeric(0)
mean_frob_stan = numeric(0)

sd_frob_glob = numeric(0)
sd_frob_stan = numeric(0)

for(i in 1:length(M)){
	mean_frob_glob[i] = mean(frob_track_glob[,7,i])
	sd_frob_glob[i] = sd(frob_track_glob[,7,i])
	mean_frob_stan[i] = mean(frob_track_stan[,7,i])
	sd_frob_stan[i] = sd(frob_track_stan[,7,i])

}


final_mat_frob = matrix(c(mean_frob_glob, mean_frob_stan), byrow = TRUE, nrow = 2)
rownames(final_mat_frob) = c("Globally Centered", "STAN")
colnames(final_mat_frob) = c("M = 2", "M = 4", "M = 8", "M = 16")
final_mat_frob

sd_mat_frob = matrix(c(sd_frob_glob, sd_frob_stan), byrow = TRUE, nrow = 2)
rownames(sd_mat_frob) = c("Globally Centered", "STAN")
colnames(sd_mat_frob) = c("M = 2", "M = 4", "M = 8", "M = 16")
sd_mat_frob
