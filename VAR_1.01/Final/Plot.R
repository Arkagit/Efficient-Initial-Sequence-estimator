############################################################
####Coverage vs Comp Time################################### 
source("../VAR_func.R")
source("../Asymp_var.R")

load("dat_matices.Rdata")

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

pdf("coverage_time1.pdf", height = 6, width = 8)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(chart[,1], time[,1], xlim = c(0.43, 1), ylim = c(-2, 2),type = 'b', xlab = "Coverage",
	ylab = "log(Time)", pch = 1:dim(time)[1])
lines(chart[,2], time[,2],type = 'b', col = "red", pch = 1:dim(time)[1])
lines(chart[,3], time[,3],type = 'b', col = "purple", pch = 1:dim(time)[1])
lines(chart[,4], time[,4],type = 'b', col = "skyblue", pch = 1:dim(time)[1])
lines(chart[,5], time[,5],type = 'b', col = "brown", pch = 1:dim(time)[1])
legend("topleft", legend = subsize,
       lwd = 1, cex = 0.75, col = "black", pch = 1:(dim(time)[1]))

add_legend("topright", bty = "n", legend = c("BM", "CC - ISE", "ISE", "SVE", "CC - MLS"), 
	col = c("black", "purple", "red", "skyblue", "brown"), lty = 1, cex=0.65)
dev.off()

############################################################
###ESS Plot############################################ 
source("../VAR_func.R")
source("../Asymp_var.R")

load("dat_matices.Rdata")

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


pdf("VAR_ess.pdf", height = 6, width = 6)
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
legend("topright", bty = "n",legend = c("BM", "ISE", "SVE", "CC - ISE", "CC - MLS", "True"), 
  col = c("black",  "red", "skyblue", "purple", "brown", "black"), lty = c(1,1,1,1,1,2), cex=0.65)

dev.off()

############################################################
####Frobenius Norm Plot#################################### 
source("../VAR_func.R")
source("../Asymp_var.R")

load("dat_matices.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

true_norm <- norm(true.sig.gen(p = p, omega = omega, phi = phi)$final.cov, type = "F")
true_norm
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
    norm_track_bm[i,j] = norm(cover[[i]]$estimates[[5*(j-1) + 1]], type = "F")
    norm_track_ise[i,j] = norm(cover[[i]]$estimates[[5*(j-1) + 2]], type = "F")
    norm_track_cc[i,j] = norm(cover[[i]]$estimates[[5*(j-1) + 3]], type = "F")
    norm_track_sve[i,j] = norm(cover[[i]]$estimates[[5*(j-1) + 4]], type = "F")
    norm_track_mls[i,j] = norm(cover[[i]]$estimates[[5*(j-1) + 5]], type = "F")
  }
}

se_norm_bm <- apply(norm_track_bm, 2, sd)/sqrt(B)
se_norm_ise <- apply(norm_track_ise, 2, sd)/sqrt(B)
se_norm_cc <- apply(norm_track_cc, 2, sd)/sqrt(B)
se_norm_sve <- apply(norm_track_sve, 2, sd)/sqrt(B)
se_norm_mls <- apply(norm_track_mls, 2, sd)/sqrt(B)

#subsize = log(subsize)/log(10)


pdf("VAR_Frob.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(norm_track_bm), type = "l", xlab = "Chain length",
  ylim = c(4000, 11000), log = 'x', ylab = "Frobenius norm")
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
abline(h = true_norm, lty = 2)
legend("bottomright", bty = "n",legend = c("BM", "ISE", "SVE", "CC - ISE", "CC - MLS", "True"), 
  col = c("black",  "red", "skyblue", "purple", "brown", "black"), lty = c(1,1,1,1,1,2), cex=0.75)

dev.off()

############################################################
####Eigen values sup dif Plot################################ 
source("../VAR_func.R")
source("../Asymp_var.R")

load("dat_matices.Rdata")

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


pdf("VAR_eigen_bias.pdf", height = 6, width = 6)
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

legend("topright", bty = "n",legend = c("BM", "ISE", "SVE", "CC - ISE", "CC - MLS"), 
  col = c("black",  "red", "brown", "purple","skyblue"), lty = 1, cex=0.75)

dev.off()


############################################################
##Computational Time Plot#################################### 

source("../VAR_func.R")
source("../Asymp_var.R")

load("dat_matices.Rdata")

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
    time_track_bm[i,j] = cover[[i]]$Time[j,1]
    time_track_ise[i,j] = cover[[i]]$Time[j,2]
    time_track_cc[i,j] = cover[[i]]$Time[j,3]
    time_track_sve[i,j] = cover[[i]]$Time[j,4]
    time_track_mls[i,j] = cover[[i]]$Time[j,5]
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

pdf("VAR_comptime.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(time_track_bm), type = "l", xlab = "Chain length",
  ylim = c(0, 80), log = 'x', ylab = "Computational Time (sec)")
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

legend("topleft", bty = "n",legend = c("BM", "ISE",  "SVE", "CC - ISE", "CC - MLS"), 
  col = c("black",  "red",  "skyblue","purple","brown"), lty = 1, cex=0.75)

dev.off()



############################################################
##Computational Time Plot#################################### 

source("../VAR_func.R")
source("../Asymp_var.R")

load("dat_matices.Rdata")

nloops = length(subsize)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

Trunc_ise = matrix(0, nrow = B, ncol = nloops)
Trunc_cc = matrix(0, nrow = B, ncol = nloops)

for(i in 1:B){
  for(j in 1:nloops){
    Trunc_ise[i,j] = log(subsize[j]*(cover[[i]]$Truncation[j,1])*(p*p))
    Trunc_cc[i,j] = log(subsize[j]*log(subsize[j])*p)
  }
}

se_ise <- apply(Trunc_ise, 2, sd)/sqrt(B)
se_cc <- apply(Trunc_cc, 2, sd)/sqrt(B)


pdf("VAR_theoretical_complexity.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

plot(subsize, colMeans(Trunc_ise), type = "l", xlab = "Chain length",
   ylim = c(10, 25), ylab = "Log Theoretical Complexity", col = "red")
segments(x0 = subsize, y0 = colMeans(Trunc_ise) - 1.96*se_ise, 
  y1 = colMeans(Trunc_ise) + 1.96*se_ise, col = "red")

lines(subsize, colMeans(Trunc_cc), type = "l", xlab = "Chain length",
   ylab = "Theoretical Time/(n d)", col = "purple")
segments(x0 = subsize, y0 = colMeans(Trunc_cc) - 1.96*se_cc, 
  y1 = colMeans(Trunc_cc) + 1.96*se_cc, col = "purple")

legend("bottomright", bty = "n",legend = c("ISE complexity", "CC-Geyer complexity"), 
  col = c("red", "purple"), lty = 1, cex=0.75)


dev.off()


