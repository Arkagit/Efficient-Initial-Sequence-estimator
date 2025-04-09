############################################################
####Coverage vs Comp Time################################### 
source("Examples/VAR_1.01/VAR_func.R")
source("Examples/VAR_1.01/Asymp_var.R")

load("Examples/VAR_1.01/Results/dat_matices.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

truth <- true.sig.gen(p = p, omega = omega, phi = phi)

nloops = length(subsize)


chart = matrix(0, nrow = length(subsize), ncol = 5)
time = matrix(0, nrow = length(subsize), ncol = 5)

time_track_bm = matrix(0, nrow = B, ncol = nloops)
time_track_ise = matrix(0, nrow = B, ncol = nloops)
time_track_cc = matrix(0, nrow = B, ncol = nloops)
time_track_sve = matrix(0, nrow = B, ncol = nloops)
time_track_mls = matrix(0, nrow = B, ncol = nloops)

se_time_bm <- numeric(length(subsize))
se_time_ise <- numeric(length(subsize))
se_time_cc <- numeric(length(subsize))
se_time_sve <- numeric(length(subsize))
se_time_mls <- numeric(length(subsize))

for(i in 1:B){
  chart = chart + (cover[[i]]$count_mat)/B
  time = time + log(cover[[i]]$Time)/(log(10)*B)
  for(j in 1:length(subsize)){
    time_track_bm[i,j] = log(cover[[i]]$Time[j,1])/log(10)
    time_track_ise[i,j] = log(cover[[i]]$Time[j,2])/log(10)
    time_track_cc[i,j] = log(cover[[i]]$Time[j,3])/log(10)
    time_track_sve[i,j] = log(cover[[i]]$Time[j,4])/log(10)
    time_track_mls[i,j] = log(cover[[i]]$Time[j,5])/log(10)
  }
}

se_time_bm <- apply(time_track_bm, 2, sd)/sqrt(B)
se_time_ise <- apply(time_track_ise, 2, sd)/sqrt(B)
se_time_cc <- apply(time_track_cc, 2, sd)/sqrt(B)
se_time_sve <- apply(time_track_sve, 2, sd)/sqrt(B)
se_time_mls <- apply(time_track_mls, 2, sd)/sqrt(B)
chart;time




add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#time = log(time)/ log(10)

pdf("Plots/coverage_time1.pdf", height = 6, width = 8)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(chart[,1], time[,1], xlim = c(0.43, 1), ylim = c(-2, 2),type = 'b', xlab = "Coverage",
	ylab = "log(Time)", pch = 1:dim(time)[1])
lines(chart[,2], time[,2],type = 'b', col = "red", pch = 1:dim(time)[1])
lines(chart[,3], time[,3],type = 'b', col = "purple", pch = 1:dim(time)[1])
lines(chart[,4], time[,4],type = 'b', col = "skyblue", pch = 1:dim(time)[1])
lines(chart[,5], time[,5],type = 'b', col = "brown", pch = 1:dim(time)[1])
legend("topleft", bty = "n", legend = subsize,
       lwd = 1, cex = 0.75, col = "black", pch = 1:(dim(time)[1]))

add_legend("topright", bty = "n", legend = c("BM", "CC - ISE", "mISE", "SVE", "CC - MLS"), 
	col = c("black", "purple", "red", "skyblue", "brown"), lty = 1, cex=0.65)
dev.off()

############################################################
###ESS Plot############################################ 
source("Examples/VAR_1.01/VAR_func.R")
source("Examples/VAR_1.01/Asymp_var.R")

load("Examples/VAR_1.01/Results/dat_matices.Rdata")


add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

ess_true <- (det(true.sig.gen(p = p, omega = omega, phi = phi)$tar.var)/det(true.sig.gen(p = p, omega = omega, phi = phi)$final.cov))^(1/p)
ess_true
nloops = length(subsize)
names = subsize

ess_track_bm = matrix(0, nrow = B, ncol = nloops)
ess_track_ise = matrix(0, nrow = B, ncol = nloops)
ess_track_cc = matrix(0, nrow = B, ncol = nloops)
ess_track_sve = matrix(0, nrow = B, ncol = nloops)
ess_track_mls = matrix(0, nrow = B, ncol = nloops)

se_ess_bm <- numeric(length(subsize))
se_ess_ise <- numeric(length(subsize))
se_ess_cc <- numeric(length(subsize))
se_ess_sve <- numeric(length(subsize))
se_ess_mls <- numeric(length(subsize))

for(i in 1:B){
  for(j in 1:nloops){
    ess_track_bm[i,j] = cover[[i]]$ESS[[5*(j-1) + 1]]/subsize[j]
    ess_track_ise[i,j] = cover[[i]]$ESS[[5*(j-1) + 2]]/subsize[j]
    ess_track_cc[i,j] = cover[[i]]$ESS[[5*(j-1) + 3]]/subsize[j]
    ess_track_sve[i,j] = cover[[i]]$ESS[[5*(j-1) + 4]]/subsize[j]
    ess_track_mls[i,j] = cover[[i]]$ESS[[5*(j-1) + 5]]/subsize[j]
  }
}

se_ess_bm <- apply(ess_track_bm, 2, sd)/sqrt(B)
se_ess_ise <- apply(ess_track_ise, 2, sd)/sqrt(B)
se_ess_cc <- apply(ess_track_cc, 2, sd)/sqrt(B)
se_ess_sve <- apply(ess_track_sve, 2, sd)/sqrt(B)
se_ess_mls <- apply(ess_track_mls, 2, sd)/sqrt(B)

#subsize = log(subsize)/log(10)


pdf("Plots/VAR_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track_bm), type = "l", xlab = "Chain length",
   ylim = c(0.020, 0.040), log = 'x',  ylab = "ESS/n")
segments(x0 = subsize, y0 = colMeans(ess_track_bm) - 1.96*se_ess_bm, 
  y1 = colMeans(ess_track_bm) + 1.96*se_ess_bm)

lines(subsize, colMeans(ess_track_ise), col = "red")
segments(x0 = subsize, y0 = colMeans(ess_track_ise) - 1.96*se_ess_ise, 
  y1 = colMeans(ess_track_ise) + 1.96*se_ess_bm, col = "red")

lines(subsize, colMeans(ess_track_cc), col = "purple")
segments(x0 = subsize, y0 = colMeans(ess_track_cc) - 1.96*se_ess_cc, 
  y1 = colMeans(ess_track_cc) + 1.96*se_ess_cc, col = "purple")

lines(subsize, colMeans(ess_track_sve), col = "skyblue")
segments(x0 = subsize, y0 = colMeans(ess_track_sve) - 1.96*se_ess_sve, 
  y1 = colMeans(ess_track_sve) + 1.96*se_ess_sve, col = "skyblue")

lines(subsize, colMeans(ess_track_mls), col = "brown")
segments(x0 = subsize, y0 = colMeans(ess_track_mls) - 1.96*se_ess_mls, 
  y1 = colMeans(ess_track_mls) + 1.96*se_ess_mls, col = "brown")
abline(h = ess_true, lty = 2)
legend("topright", bty = "n",legend = c("BM", "mISE", "SVE", "CC - ISE", "CC - MLS", "True"), 
  col = c("black",  "red", "skyblue", "purple", "brown", "black"), lty = c(1,1,1,1,1,2), cex=0.65)

dev.off()

############################################################
####Frobenius Norm Plot#################################### 
source("Examples/VAR_1.01/VAR_func.R")
source("Examples/VAR_1.01/Asymp_var.R")

load("Examples/VAR_1.01/Results/dat_matices.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

true_norm <- norm(true.sig.gen(p = p, omega = omega, phi = phi)$final.cov, type = "F")
Truth <- true.sig.gen(p = p, omega = omega, phi = phi)$final.cov

nloops = length(subsize)

norm_track_bm = matrix(0, nrow = B, ncol = nloops)
norm_track_ise = matrix(0, nrow = B, ncol = nloops)
norm_track_cc = matrix(0, nrow = B, ncol = nloops)
norm_track_sve = matrix(0, nrow = B, ncol = nloops)
norm_track_mls = matrix(0, nrow = B, ncol = nloops)

se_norm_bm <- numeric(length(subsize))
se_norm_ise <- numeric(length(subsize))
se_norm_cc <- numeric(length(subsize))
se_norm_sve <- numeric(length(subsize))
se_norm_mls <- numeric(length(subsize))

for(i in 1:B){
  for(j in 1:nloops){
    norm_track_bm[i,j] = norm(Truth - cover[[i]]$estimates[[5*(j-1) + 1]], type = "F")/norm(Truth, type = "F")
    norm_track_ise[i,j] = norm(Truth - cover[[i]]$estimates[[5*(j-1) + 2]], type = "F")/norm(Truth, type = "F")
    norm_track_cc[i,j] = norm(Truth - cover[[i]]$estimates[[5*(j-1) + 3]], type = "F")/norm(Truth, type = "F")
    norm_track_sve[i,j] = norm(Truth - cover[[i]]$estimates[[5*(j-1) + 4]], type = "F")/norm(Truth, type = "F")
    norm_track_mls[i,j] = norm(Truth - cover[[i]]$estimates[[5*(j-1) + 5]], type = "F")/norm(Truth, type = "F")
  }
}

se_norm_bm <- apply(norm_track_bm, 2, sd)/sqrt(B)
se_norm_ise <- apply(norm_track_ise, 2, sd)/sqrt(B)
se_norm_cc <- apply(norm_track_cc, 2, sd)/sqrt(B)
se_norm_sve <- apply(norm_track_sve, 2, sd)/sqrt(B)
se_norm_mls <- apply(norm_track_mls, 2, sd)/sqrt(B)

#subsize = log(subsize)/log(10)


pdf("Plots/VAR_Frob.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(norm_track_bm), type = "l", xlab = "Chain length",
  ylim = c(0, 0.6), log = 'x', ylab = "Relative Frobenius norm")
segments(x0 = subsize, y0 = colMeans(norm_track_bm) - 1.96*se_norm_bm, 
  y1 = colMeans(norm_track_bm) + 1.96*se_norm_bm)

lines(subsize, colMeans(norm_track_ise), col = "red")
segments(x0 = subsize, y0 = colMeans(norm_track_ise) - 1.96*se_norm_ise, 
  y1 = colMeans(norm_track_ise) + 1.96*se_norm_bm, col = "red")

lines(subsize, colMeans(norm_track_cc), col = "purple")
segments(x0 = subsize, y0 = colMeans(norm_track_cc) - 1.96*se_norm_cc, 
  y1 = colMeans(norm_track_cc) + 1.96*se_norm_cc, col = "purple")

lines(subsize, colMeans(norm_track_sve), col = "skyblue")
segments(x0 = subsize, y0 = colMeans(norm_track_sve) - 1.96*se_norm_sve, 
  y1 = colMeans(norm_track_sve) + 1.96*se_norm_sve, col = "skyblue")

lines(subsize, colMeans(norm_track_mls), col = "brown")
segments(x0 = subsize, y0 = colMeans(norm_track_mls) - 1.96*se_norm_mls, 
  y1 = colMeans(norm_track_mls) + 1.96*se_norm_mls, col = "brown")
#abline(h = true_norm, lty = 2)
legend("topright", bty = "n",legend = c("BM", "mISE", "SVE", "CC - ISE", "CC - MLS"), 
  col = c("black",  "red", "skyblue", "purple", "brown"), lty = c(1,1,1,1,1), cex=0.75)

dev.off()

############################################################
####Eigen values sup dif Plot################################ 
source("Examples/VAR_1.01/VAR_func.R")
source("Examples/VAR_1.01/Asymp_var.R")

load("Examples/VAR_1.01/Results/dat_matices.Rdata")


add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

truth <- true.sig.gen(p = p, omega = omega, phi = phi)$final.cov

nloops = length(subsize)

bias_track_bm = matrix(0, nrow = B, ncol = nloops)
bias_track_ise = matrix(0, nrow = B, ncol = nloops)
bias_track_cc = matrix(0, nrow = B, ncol = nloops)
bias_track_sve = matrix(0, nrow = B, ncol = nloops)
bias_track_mls = matrix(0, nrow = B, ncol = nloops)

se_bias_bm <- numeric(length(subsize))
se_bias_ise <- numeric(length(subsize))
se_bias_cc <- numeric(length(subsize))
se_bias_sve <- numeric(length(subsize))
se_bias_mls <- numeric(length(subsize))

for(i in 1:B){
  for(j in 1:length(subsize)){
    bias_track_bm[i,j] = max(abs(eigen(cover[[i]]$estimates[[5*(j-1) + 1]])$values - eigen(truth)$value))
    bias_track_ise[i,j] = max(abs(eigen(cover[[i]]$estimates[[5*(j-1) + 2]])$values - eigen(truth)$values))
    bias_track_cc[i,j] = max(abs(eigen(cover[[i]]$estimates[[5*(j-1) + 3]])$values - eigen(truth)$values))
    bias_track_sve[i,j] = max(abs(eigen(cover[[i]]$estimates[[5*(j-1) + 4]])$values - eigen(truth)$values))
    bias_track_mls[i,j] = max(abs(eigen(cover[[i]]$estimates[[5*(j-1) + 5]])$values - eigen(truth)$values))
  }
}

se_bias_bm <- apply(bias_track_bm, 2, sd)/sqrt(B)
se_bias_ise <- apply(bias_track_ise, 2, sd)/sqrt(B)
se_bias_cc <- apply(bias_track_cc, 2, sd)/sqrt(B)
se_bias_sve <- apply(bias_track_sve, 2, sd)/sqrt(B)
se_bias_mls <- apply(bias_track_mls, 2, sd)/sqrt(B)


#subsize = log(subsize)/log(10)


pdf("Plots/VAR_eigen_bias.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(bias_track_bm), type = "l", xlab = "Chain length",
  ylim = c(0, 6000), log = 'x', ylab = "Absolute max error of Eigen values")
segments(x0 = subsize, y0 = colMeans(bias_track_bm) - 1.96*se_bias_bm, 
  y1 = colMeans(bias_track_bm) + 1.96*se_bias_bm)

lines(subsize, colMeans(bias_track_ise), col = "red")
segments(x0 = subsize, y0 = colMeans(bias_track_ise) - 1.96*se_bias_ise, 
  y1 = colMeans(bias_track_ise) + 1.96*se_bias_bm, col = "red")

lines(subsize, colMeans(bias_track_cc), col = "purple")
segments(x0 = subsize, y0 = colMeans(bias_track_cc) - 1.96*se_bias_cc, 
  y1 = colMeans(bias_track_cc) + 1.96*se_bias_cc, col = "purple")

lines(subsize, colMeans(bias_track_sve), col = "skyblue")
segments(x0 = subsize, y0 = colMeans(bias_track_sve) - 1.96*se_bias_sve, 
  y1 = colMeans(bias_track_sve) + 1.96*se_bias_sve, col = "skyblue")

lines(subsize, colMeans(bias_track_mls), col = "brown")
segments(x0 = subsize, y0 = colMeans(bias_track_mls) - 1.96*se_bias_mls, 
  y1 = colMeans(bias_track_mls) + 1.96*se_bias_mls, col = "brown")

legend("topright", bty = "n",legend = c("BM", "mISE", "SVE", "CC - ISE", "CC - MLS"), 
  col = c("black",  "red", "brown", "purple","skyblue"), lty = 1, cex=0.75)

dev.off()


############################################################
##Computational Time Plot#################################### 

source("Examples/VAR_1.01/VAR_func.R")
source("Examples/VAR_1.01/Asymp_var.R")

load("Examples/VAR_1.01/Results/dat_matices.Rdata")

truth <- true.sig.gen(p = p, omega = omega, phi = phi)

nloops = length(subsize)


chart = matrix(0, nrow = length(subsize), ncol = 5)
time = matrix(0, nrow = length(subsize), ncol = 5)

time_track_bm = matrix(0, nrow = B, ncol = nloops)
time_track_ise = matrix(0, nrow = B, ncol = nloops)
time_track_cc = matrix(0, nrow = B, ncol = nloops)
time_track_sve = matrix(0, nrow = B, ncol = nloops)
time_track_mls = matrix(0, nrow = B, ncol = nloops)

se_time_bm <- numeric(length(subsize))
se_time_ise <- numeric(length(subsize))
se_time_cc <- numeric(length(subsize))
se_time_sve <- numeric(length(subsize))
se_time_mls <- numeric(length(subsize))

for(i in 1:B){
  chart = chart + (cover[[i]]$count_mat)/B
  time = time + log(cover[[i]]$Time)/(log(10)*B)
  for(j in 1:length(subsize)){
    time_track_bm[i,j] = log(cover[[i]]$Time[j,1])
    time_track_ise[i,j] = log(cover[[i]]$Time[j,2])
    time_track_cc[i,j] = log(cover[[i]]$Time[j,3])
    time_track_sve[i,j] = log(cover[[i]]$Time[j,4])
    time_track_mls[i,j] = log(cover[[i]]$Time[j,5])
  }
}

se_time_bm <- apply(time_track_bm, 2, sd)/sqrt(B)
se_time_ise <- apply(time_track_ise, 2, sd)/sqrt(B)
se_time_cc <- apply(time_track_cc, 2, sd)/sqrt(B)
se_time_sve <- apply(time_track_sve, 2, sd)/sqrt(B)
se_time_mls <- apply(time_track_mls, 2, sd)/sqrt(B)
chart;time




add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#subsize = log(subsize)/log(10)

pdf("Plots/VAR_comptime.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(time_track_bm), type = "l", xlab = "Chain length",
  ylim = c(-6, 7), log = 'x', ylab = "Log Computational Time (sec)")
segments(x0 = subsize, y0 = colMeans(time_track_bm) - 1.96*se_time_bm, 
  y1 = colMeans(time_track_bm) + 1.96*se_time_bm)

lines(subsize, colMeans(time_track_ise), col = "red")
segments(x0 = subsize, y0 = colMeans(time_track_ise) - 1.96*se_time_ise, 
  y1 = colMeans(time_track_ise) + 1.96*se_time_bm, col = "red")

lines(subsize, colMeans(time_track_cc), col = "purple")
segments(x0 = subsize, y0 = colMeans(time_track_cc) - 1.96*se_time_cc, 
  y1 = colMeans(time_track_cc) + 1.96*se_time_cc, col = "purple")

lines(subsize, colMeans(time_track_sve), col = "skyblue")
segments(x0 = subsize, y0 = colMeans(time_track_sve) - 1.96*se_time_sve, 
  y1 = colMeans(time_track_sve) + 1.96*se_time_sve, col = "skyblue")

lines(subsize, colMeans(time_track_mls), col = "brown")
segments(x0 = subsize, y0 = colMeans(time_track_mls) - 1.96*se_time_mls, 
  y1 = colMeans(time_track_mls) + 1.96*se_time_mls, col = "brown")

legend("topleft", bty = "n",legend = c("BM", "mISE",  "SVE", "CC - ISE", "CC - MLS"), 
  col = c("black",  "red",  "skyblue","purple","brown"), lty = 1, cex=0.75)

dev.off()




############################################################
## Truncation Time Plot #################################### 

source("Examples/VAR_1.01/VAR_func.R")
source("Examples/VAR_1.01/Asymp_var.R")

load("Examples/VAR_1.01/Results/dat_matices.Rdata")

library(latex2exp)

nloops = length(subsize)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

Trunc_ise = matrix(0, nrow = B, ncol = nloops)

for(i in 1:B){
  for(j in 1:nloops){
    Trunc_ise[i,j] = (cover[[i]]$Truncation[j,1])
  }
}

Trunc_ise = cbind(log(Trunc_ise))
l.subsize = log(subsize)

se_ise <- apply(Trunc_ise, 2, sd)/sqrt(B)


pdf("Plots/VAR_theoretical_complexity.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

plot(l.subsize, colMeans(Trunc_ise), type = "b", col = "red",
 pch = 16, xlab = TeX(r'($\log (n)$)'), ylab = TeX(r'(Log of $t_{n}$)'))
segments(x0 = l.subsize, y0 = colMeans(Trunc_ise) - 1.96*se_ise, 
  y1 = colMeans(Trunc_ise) + 1.96*se_ise, col = "red")


dev.off()



#################################################
rm()
################################################# Calcium Example

load("Examples/Spike_train/Results/Output_d.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

bm_time = matrix(0, nrow = repet, ncol = length(N))
ise_time = matrix(0, nrow = repet, ncol = length(N))
sve_time = matrix(0, nrow = repet, ncol = length(N))
cc_time = matrix(0, nrow = repet, ncol = length(N))
mls_time = matrix(0, nrow = repet, ncol = length(N))


for(i in 1:repet){
  bm_time[i,] = log(as.numeric(Table[[i]][[2]][[1]]))
  ise_time[i,] = log(as.numeric(Table[[i]][[2]][[2]]))
  sve_time[i,] = log(as.numeric(Table[[i]][[2]][[3]]))
  cc_time[i,] = log(as.numeric(Table[[i]][[2]][[4]]))
  mls_time[i,] = log(as.numeric(Table[[i]][[2]][[5]]))
}

se_time_bm <- apply(bm_time, 2, sd)/sqrt(repet)
se_time_ise <- apply(ise_time, 2, sd)/sqrt(repet)
se_time_cc <- apply(cc_time, 2, sd)/sqrt(repet)
se_time_sve <- apply(sve_time, 2, sd)/sqrt(repet)
se_time_mls <- apply(mls_time, 2, sd)/sqrt(repet)

#N = log(N)/log(10)


pdf("Plots/Calcium_spike_comptime.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(N, colMeans(bm_time),col = "black", xlab = "Chain Length", ylab = "Log Computational time (sec)", 
  ylim = c(-10, 10), log = 'x', type = "l")
segments(x0 = N, y0 = colMeans(bm_time) - 1.96*se_time_bm, 
  y1 = colMeans(bm_time) + 1.96*se_time_bm)

lines(N, colMeans(sve_time), col = "skyblue", lty = 1)
segments(x0 = N, y0 = colMeans(sve_time) - 1.96*se_time_sve, 
  y1 = colMeans(sve_time) + 1.96*se_time_sve, col = "skyblue")

lines(N, colMeans(cc_time), col = "purple", lty = 1)
segments(x0 = N, y0 = colMeans(cc_time) - 1.96*se_time_cc, 
  y1 = colMeans(cc_time) + 1.96*se_time_cc, col = "purple")

lines(N, colMeans(ise_time), col = "red", lty = 1)
segments(x0 = N, y0 = colMeans(ise_time) - 1.96*se_time_ise, 
  y1 = colMeans(ise_time) + 1.96*se_time_ise, col = "red")

lines(N, colMeans(mls_time), col = "brown", lty = 1)
segments(x0 = N, y0 = colMeans(mls_time) - 1.96*se_time_mls, 
  y1 = colMeans(mls_time) + 1.96*se_time_mls, col = "brown")

legend("topleft", bty = "n",legend = c("BM", "SVE",  "mISE", "CC - ISE","CC - MLS"),
 col = c("black", "skyblue",  "red", "purple","brown"), cex = 0.8,lty = 1)
dev.off()

#####################################################################
###### ESS Plot #####################################################

load("Examples/Spike_train/Results/Output_d.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

p = dim(Table[[1]][[1]][[1]])[1]

bm_ess = matrix(0, nrow = repet, ncol = length(N))
ise_ess = matrix(0, nrow = repet, ncol = length(N))
sve_ess = matrix(0, nrow = repet, ncol = length(N))
cc_ess = matrix(0, nrow = repet, ncol = length(N))
mls_ess = matrix(0, nrow = repet, ncol = length(N))

for(i in 1:repet){
  for(j in 1:length(N)){
    bm_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 2]]))^(1/p)
    ise_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 3]]))^(1/p)
    sve_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 4]]))^(1/p)
    cc_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 5]]))^(1/p)
    mls_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 6]]))^(1/p)
  }
}


se_ess_bm <- apply(bm_ess, 2, sd)/sqrt(repet)
se_ess_sve <- apply(ise_ess, 2, sd)/sqrt(repet)
se_ess_cc <- apply(cc_ess, 2, sd)/sqrt(repet)
se_ess_ise <- apply(sve_ess, 2, sd)/sqrt(repet)
se_ess_mls <- apply(mls_ess, 2, sd)/sqrt(repet)

#N = log(N)/log(10)

pdf("Plots/Calcium_spike_ess.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(N, colMeans(bm_ess),col = "black", xlab = "Chain Length", ylab = "ESS/n", 
  ylim = c(0.035, 0.075), log = 'x', type = "l")
segments(x0 = N, y0 = colMeans(bm_ess) - 1.96*se_ess_bm, 
  y1 = colMeans(bm_ess) + 1.96*se_ess_bm)

lines(N, colMeans(sve_ess), col = "skyblue", lty = 1)
segments(x0 = N, y0 = colMeans(sve_ess) - 1.96*se_ess_sve, 
  y1 = colMeans(sve_ess) + 1.96*se_ess_sve, col = "skyblue")

lines(N, colMeans(cc_ess), col = "purple", lty = 1)
segments(x0 = N, y0 = colMeans(cc_ess) - 1.96*se_ess_cc, 
  y1 = colMeans(cc_ess) + 1.96*se_ess_cc, col = "purple")

lines(N, colMeans(ise_ess), col = "red", lty = 1)
segments(x0 = N, y0 = colMeans(ise_ess) - 1.96*se_ess_ise, 
  y1 = colMeans(ise_ess) + 1.96*se_ess_ise, col = "red")

lines(N, colMeans(mls_ess), col = "brown", lty = 1)
segments(x0 = N, y0 = colMeans(mls_ess) - 1.96*se_ess_mls, 
  y1 = colMeans(mls_ess) + 1.96*se_ess_mls, col = "brown")

#abline(h = 0, lty = 2)
legend("topright", bty = "n",legend = c("BM", "SVE", "mISE", "CC - ISE", "CC - MLS"),
 col = c("black", "skyblue", "red", "purple", "brown"), cex = 0.8,lty = 1)
dev.off()


#################################################
rm()
################################################# Parallel Chain Example

set.seed(100)

source("Examples/Parallel_chain/Asymp_var.R")
source("Examples/Parallel_chain/VAR_func.R")
source("Examples/Parallel_chain/global_var.R")
source("Examples/Parallel_chain/Stan_code.R")


load("Examples/Parallel_chain/Results/parallel_est.Rdata")

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
####### n = 1e5

df = data.frame(
  id = c(rep("global", B), rep("stan", B)),
  "2" = c(frob_track_glob[,5,1], frob_track_stan[,5,1]),
  "4" = c(frob_track_glob[,5,2], frob_track_stan[,5,2]),
  "8" = c(frob_track_glob[,5,3], frob_track_stan[,5,3]),
  "16" = c(frob_track_glob[,5,4], frob_track_stan[,5,4]))


pdf("Plots/Boxplot_frob_1e5.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

boxplot(df[,-1], boxfill = NA, border = NA, names = M, xlab = "Number of parallel chains (n = 1e5)", 
  ylab = "Relative Frobenius Norm", ylim = c(0, 0.3)) #invisible boxes - only axes and plot area
boxplot(df[df$id=="global", -1], xaxt = "n", add = TRUE, boxfill="cadetblue1", 
  boxwex=0.25, at = 1:ncol(df[,-1]) - 0.15) #shift these left by -0.15
boxplot(df[df$id=="stan", -1], xaxt = "n", add = TRUE, boxfill="firebrick1", 
  boxwex=0.25, at = 1:ncol(df[,-1]) + 0.15) #shift to the right by +0.15
legend("topright", bty = "n", legend = c("Globally Centered", "STAN"), 
       fill = c("cadetblue1", "firebrick1"), border = "black")
dev.off()


#################### boxplot ##############
####### n = 1e4

df = data.frame(
  id = c(rep("global", B), rep("stan", B)),
  "2" = c(frob_track_glob[,1,1], frob_track_stan[,1,1]),
  "4" = c(frob_track_glob[,1,2], frob_track_stan[,1,2]),
  "8" = c(frob_track_glob[,1,3], frob_track_stan[,1,3]),
  "16" = c(frob_track_glob[,1,4], frob_track_stan[,1,4]))


pdf("Plots/Boxplot_frob_1e3.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

boxplot(df[,-1], boxfill = NA, border = NA, names = M, xlab = "Number of parallel chains (n = 1e3)", 
  ylab = "Relative Frobenius Norm", ylim = c(0, 4.1)) #invisible boxes - only axes and plot area
boxplot(df[df$id=="global", -1], xaxt = "n", add = TRUE, boxfill="cadetblue1", 
  boxwex=0.25, at = 1:ncol(df[,-1]) - 0.15) #shift these left by -0.15
boxplot(df[df$id=="stan", -1], xaxt = "n", add = TRUE, boxfill="firebrick1", 
  boxwex=0.25, at = 1:ncol(df[,-1]) + 0.15) #shift to the right by +0.15
legend("topright", bty = "n", legend = c("Globally Centered", "STAN"), 
       fill = c("cadetblue1", "firebrick1"), border = "black")
dev.off()

