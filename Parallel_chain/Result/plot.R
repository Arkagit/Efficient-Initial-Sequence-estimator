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
	mean_ess_glob[i] = mean(ess_track_glob[,5,i])
	sd_ess_glob[i] = sd(ess_track_glob[,5,i])
	mean_ess_stan[i] = mean(ess_track_stan[,5,i])
	sd_ess_stan[i] = sd(ess_track_stan[,5,i])

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
      frob_track_glob[i,j,k] = norm(Truth - Estimate[[i]][[j]][[1]][[k]], type = "F")/norm(Truth, type = "F")
      frob_track_stan[i,j,k] = norm(Truth - Estimate[[i]][[j]][[2]][[k]], type = "F")/norm(Truth, type = "F")
    }
  }
}



mean_frob_glob = numeric(0)
mean_frob_stan = numeric(0)

sd_frob_glob = numeric(0)
sd_frob_stan = numeric(0)

for(i in 1:length(M)){
	mean_frob_glob[i] = mean(frob_track_glob[,5,i])
	sd_frob_glob[i] = sd(frob_track_glob[,5,i])
	mean_frob_stan[i] = mean(frob_track_stan[,5,i])
	sd_frob_stan[i] = sd(frob_track_stan[,5,i])

}


final_mat_frob = matrix(c(mean_frob_glob, mean_frob_stan), byrow = TRUE, nrow = 2)
rownames(final_mat_frob) = c("Globally Centered", "STAN")
colnames(final_mat_frob) = c("M = 2", "M = 4", "M = 8", "M = 16")
final_mat_frob

sd_mat_frob = matrix(c(sd_frob_glob, sd_frob_stan), byrow = TRUE, nrow = 2)
rownames(sd_mat_frob) = c("Globally Centered", "STAN")
colnames(sd_mat_frob) = c("M = 2", "M = 4", "M = 8", "M = 16")
sd_mat_frob

#################### boxplot ##############

#boxplot(frob_track_glob[,7,i])

#plot(frob_track_glob[,7,i])

df = data.frame(
  id = c(rep("global", B), rep("stan", B)),
  "2" = c(frob_track_glob[,5,1], frob_track_stan[,5,1]),
  "4" = c(frob_track_glob[,5,2], frob_track_stan[,5,2]),
  "8" = c(frob_track_glob[,5,3], frob_track_stan[,5,3]),
  "16" = c(frob_track_glob[,5,4], frob_track_stan[,5,4]))


pdf("Boxplot.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

boxplot(df[,-1], boxfill = NA, border = NA, names = M, xlab = "Number of parallel chains", 
  ylab = "Relative Frobenius Norm", ylim = c(0, 0.2)) #invisible boxes - only axes and plot area
boxplot(df[df$id=="global", -1], xaxt = "n", add = TRUE, boxfill="cadetblue1", 
  boxwex=0.25, at = 1:ncol(df[,-1]) - 0.15) #shift these left by -0.15
boxplot(df[df$id=="stan", -1], xaxt = "n", add = TRUE, boxfill="firebrick1", 
  boxwex=0.25, at = 1:ncol(df[,-1]) + 0.15) #shift to the right by +0.15
legend("topright", legend = c("Globally Centered", "STAN"), 
       fill = c("cadetblue1", "firebrick1"), border = "black")
dev.off()

######################### ESS boxplot ######################

df1 = data.frame(
  id = c(rep("global", B), rep("stan", B)),
  "2" = c(ess_track_glob[,5,1], ess_track_stan[,5,1]),
  "4" = c(ess_track_glob[,5,2], ess_track_stan[,5,2]),
  "8" = c(ess_track_glob[,5,3], ess_track_stan[,5,3]),
  "16" = c(ess_track_glob[,5,4], ess_track_stan[,5,4]))


pdf("ESS_Boxplot.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

boxplot(df1[,-1], boxfill = NA, border = NA, names = M, xlab = "Number of parallel chains", 
  ylab = "Effective Sample Size", ylim = c(0.02, 0.03)) #invisible boxes - only axes and plot area
boxplot(df1[df1$id=="global", -1], xaxt = "n", add = TRUE, boxfill="cadetblue1", 
  boxwex=0.25, at = 1:ncol(df1[,-1]) - 0.15) #shift these left by -0.15
boxplot(df1[df1$id=="stan", -1], xaxt = "n", add = TRUE, boxfill="firebrick1", 
  boxwex=0.25, at = 1:ncol(df1[,-1]) + 0.15) #shift to the right by +0.15
legend("topright", legend = c("Globally Centered", "STAN"), 
       fill = c("cadetblue1", "firebrick1"), border = "black")
dev.off()


########################### Comparison Frob ##################

df2 = data.frame(
  id = c(rep("global", B), rep("stan", B)),
  "1e3" = c(frob_track_glob[,1,2], frob_track_stan[,1,2]),
  "5e3" = c(frob_track_glob[,2,2], frob_track_stan[,2,2]),
  "1e4" = c(frob_track_glob[,3,2], frob_track_stan[,3,2]),
  "5e4" = c(frob_track_glob[,4,2], frob_track_stan[,4,2]),
  "1e5" = c(frob_track_glob[,5,2], frob_track_stan[,5,2]))


pdf("Comp_Boxplot_frob.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

boxplot(df2[,-1], boxfill = NA, border = NA, names = subsize, xlab = "Length of parallel chains (m = 4)", 
  ylab = "Relative Frobenious Norm", ylim = c(0, 2)) #invisible boxes - only axes and plot area
boxplot(df2[df2$id=="global", -1], xaxt = "n", add = TRUE, boxfill="cadetblue1", 
  boxwex=0.25, at = 1:ncol(df2[,-1]) - 0.15) #shift these left by -0.15
boxplot(df2[df2$id=="stan", -1], xaxt = "n", add = TRUE, boxfill="firebrick1", 
  boxwex=0.25, at = 1:ncol(df2[,-1]) + 0.15) #shift to the right by +0.15
legend("topright", legend = c("Globally Centered", "STAN"), 
       fill = c("cadetblue1", "firebrick1"), border = "black")
dev.off()



########################### Comparison ESS ##################

df3 = data.frame(
  id = c(rep("global", B), rep("stan", B)),
  "1e3" = c(ess_track_glob[,1,2], ess_track_stan[,1,2]),
  "5e3" = c(ess_track_glob[,2,2], ess_track_stan[,2,2]),
  "1e4" = c(ess_track_glob[,3,2], ess_track_stan[,3,2]),
  "5e4" = c(ess_track_glob[,4,2], ess_track_stan[,4,2]),
  "1e5" = c(ess_track_glob[,5,2], ess_track_stan[,5,2]))


pdf("Comp_Boxplot_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

boxplot(df3[,-1], boxfill = NA, border = NA, names = subsize, xlab = "Length of parallel chains (m = 4)", 
  ylab = "Effective Sample Size", ylim = c(0, 0.03)) #invisible boxes - only axes and plot area
boxplot(df3[df3$id=="global", -1], xaxt = "n", add = TRUE, boxfill="cadetblue1", 
  boxwex=0.25, at = 1:ncol(df3[,-1]) - 0.15) #shift these left by -0.15
boxplot(df3[df3$id=="stan", -1], xaxt = "n", add = TRUE, boxfill="firebrick1", 
  boxwex=0.25, at = 1:ncol(df3[,-1]) + 0.15) #shift to the right by +0.15
legend("topright", legend = c("Globally Centered", "STAN"), 
       fill = c("cadetblue1", "firebrick1"), border = "black")
dev.off()
