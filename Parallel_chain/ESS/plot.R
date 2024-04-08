set.seed(100)

source("../Asymp_var.R")
source("../VAR_func.R")
source("../global_var.R")
source("../Stan_code.R")


load("parallel_est.Rdata")

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

ess_track_glob = matrix(0, nrow = B, ncol = length(subsize))
ess_track_stan = matrix(0, nrow = B, ncol = length(subsize))

se_ess_glob <- numeric(length(subsize))
se_ess_stan <- numeric(length(subsize))

for(i in 1:B){
  for(j in 1:length(subsize)){
    ess_track_glob[i,j] = Estimate[[i]][[1]][[j]]/subsize[j]
    ess_track_stan[i,j] = Estimate[[i]][[2]][[j]]/subsize[j]
  }
}

se_ess_glob <- apply(ess_track_glob, 2, sd)/sqrt(B)
se_ess_stan <- apply(ess_track_stan, 2, sd)/sqrt(B)

subsize = log(subsize)/log(10)


pdf("parallel_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track_glob), type = "l", xlab = "Log chain length",
   ylim = c(0.020, 0.040), ylab = "ESS/n")
segments(x0 = subsize, y0 = colMeans(ess_track_glob) - 1.96*se_ess_glob, 
  y1 = colMeans(ess_track_glob) + 1.96*se_ess_glob)

lines(subsize, colMeans(ess_track_stan), col = "blue")
segments(x0 = subsize, y0 = colMeans(ess_track_stan) - 1.96*se_ess_stan, 
  y1 = colMeans(ess_track_stan) + 1.96*se_ess_glob, col = "blue")

abline(h = ess_true, lty = 2, col = "red")
legend("topright", bty = "n",legend = c("CC - Globally Centered", "CC - Stan", "True"), 
  col = c("black", "blue", "red"), lty = c(1,1,2), cex=0.65)

dev.off()