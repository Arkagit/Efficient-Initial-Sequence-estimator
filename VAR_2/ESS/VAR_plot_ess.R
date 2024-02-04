################ Sourcing Data files
source("VAR_ess.R")
load("ESS_data_2.Rdata")
source("../Asymp_var.R")


# Initialization
#p <- 10
#rho <- 0.99
#phi <- diag(rep(rho,p))
#omega <- matrix(.9, nrow = p, ncol = p)
#diag(omega) <- 1
#omega <- omega^(abs(col(omega)-row(omega)))

#truth <- true.sig.gen(p = p, omega = omega, phi = phi)



############ ESS and Frobenius norms of different data lengths in matrices
ess_track_bm = matrix(0, nrow = B, ncol = nloops)
#ess_track_lug = matrix(0, nrow = B, ncol = nloops)
ess_track_ise = matrix(0, nrow = B, ncol = nloops)
ess_track_cc = matrix(0, nrow = B, ncol = nloops)
ess_track_sve = matrix(0, nrow = B, ncol = nloops)
ess_track_mls = matrix(0, nrow = B, ncol = nloops)

fro_track_bm = matrix(0, nrow = B, ncol = nloops)
#fro_track_lug = matrix(0, nrow = B, ncol = nloops)
fro_track_ise = matrix(0, nrow = B, ncol = nloops)
fro_track_cc = matrix(0, nrow = B, ncol = nloops)
fro_track_sve = matrix(0, nrow = B, ncol = nloops)
fro_track_mls = matrix(0, nrow = B, ncol = nloops)


for(i in 1:B){
	for(j in 1:nloops){
		ess_track_bm[i,j] = as.numeric(Table[[i]][[j]][[1]][[1]])
		#ess_track_lug[i,j] = as.numeric(Table[[i]][[j]][[1]][[2]])
		ess_track_ise[i,j] = as.numeric(Table[[i]][[j]][[1]][[2]])
		ess_track_cc[i,j] = as.numeric(Table[[i]][[j]][[1]][[3]])
		ess_track_sve[i,j] = as.numeric(Table[[i]][[j]][[1]][[4]])
		ess_track_mls[i,j] = as.numeric(Table[[i]][[j]][[1]][[5]])

		fro_track_bm[i,j] = as.numeric(Table[[i]][[j]][[2]][[1]])
		#fro_track_lug[i,j] = as.numeric(Table[[i]][[j]][[2]][[2]])
		fro_track_ise[i,j] = as.numeric(Table[[i]][[j]][[2]][[2]])
		fro_track_cc[i,j] = as.numeric(Table[[i]][[j]][[2]][[3]])
		fro_track_sve[i,j] = as.numeric(Table[[i]][[j]][[2]][[4]])
		fro_track_mls[i,j] = as.numeric(Table[[i]][[j]][[2]][[5]])
	}
}

############# Calculating Standard Errors
se_ess_bm <- apply(ess_track_bm, 2, sd)/sqrt(nloops)
#se_ess_lug <- apply(ess_track_lug, 2, sd)/sqrt(nloops)
se_ess_ise <- apply(ess_track_ise, 2, sd)/sqrt(nloops)
se_ess_cc <- apply(ess_track_cc, 2, sd)/sqrt(nloops)
se_ess_sve <- apply(ess_track_sve, 2, sd)/sqrt(nloops)
se_ess_mls <- apply(ess_track_mls, 2, sd)/sqrt(nloops)

se_fro_bm <- apply(fro_track_bm, 2, sd)/sqrt(nloops)
#se_fro_lug <- apply(fro_track_lug, 2, sd)/sqrt(nloops)
se_fro_ise <- apply(fro_track_ise, 2, sd)/sqrt(nloops)
se_fro_cc <- apply(fro_track_cc, 2, sd)/sqrt(nloops)
se_fro_sve <- apply(fro_track_sve, 2, sd)/sqrt(nloops)
se_fro_mls <- apply(fro_track_mls, 2, sd)/sqrt(nloops)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


# ESS plot
pdf("plot_ess_2.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track_bm), ylim = c(0.8, 0.9),
	type = 'l', xlab = "Time", ylab = expression(paste(hat(ESS)/n)))
segments(x0 = subsize, y0 = colMeans(ess_track_bm) - 1.96*se_ess_bm, 
	y1 = colMeans(ess_track_bm) + 1.96*se_ess_bm)

#lines(subsize, colMeans(ess_track_lug), col = "red" )
#segments(x0 = subsize, y0 = colMeans(ess_track_lug) - 1.96*se_ess_lug, 
#	y1 = colMeans(ess_track_lug) + 1.96*se_ess_lug, col = "red")

lines(subsize, colMeans(ess_track_cc), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track_cc) - 1.96*se_ess_cc, 
	y1 = colMeans(ess_track_cc) + 1.96*se_ess_cc, col = "purple")

lines(subsize, colMeans(ess_track_ise), col = "green" )
segments(x0 = subsize, y0 = colMeans(ess_track_ise) - 1.96*se_ess_ise, 
	y1 = colMeans(ess_track_ise) + 1.96*se_ess_ise, col = "green")

lines(subsize,colMeans(ess_track_sve), col = "skyblue" )
segments(x0 = subsize, y0 = colMeans(ess_track_sve) - 1.96*se_ess_sve, 
	y1 = colMeans(ess_track_sve) + 1.96*se_ess_sve, col = "skyblue")

lines(subsize,colMeans(ess_track_mls), col = "brown" )
segments(x0 = subsize, y0 = colMeans(ess_track_mls) - 1.96*se_ess_mls, 
	y1 = colMeans(ess_track_mls) + 1.96*se_ess_mls, col = "brown")

abline(h = ess_true, lty = 2)
add_legend("topright", bty = "n",legend = c("BM", "New ISE (Geyer)", "ISE", "SVE", "New ISE (MLS)"), 
	col = c("black", "purple", "green", "skyblue", "brown"), lty = 1, cex=0.5)
dev.off()


# Frobenius norm plot
pdf("plot_frob_2.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(fro_track_bm),type = 'l', xlab = "Time", ylim = c(0.9, 1),
	ylab = "Frobenius Norm")
segments(x0 = subsize, y0 = colMeans(fro_track_bm) - 1.96*se_fro_bm, 
	y1 = colMeans(fro_track_bm) + 1.96*se_fro_bm)

#lines(subsize,colMeans(fro_track_lug), col = "red" )
#segments(x0 = subsize, y0 = colMeans(fro_track_lug) - 1.96*se_fro_lug, 
#	y1 = colMeans(fro_track_lug) + 1.96*se_fro_lug, col = "red")

lines(subsize,colMeans(fro_track_cc), col = "purple" )
segments(x0 = subsize, y0 = colMeans(fro_track_cc) - 1.96*se_fro_cc, 
	y1 = colMeans(fro_track_cc) + 1.96*se_fro_cc, col = "purple")

lines(subsize,colMeans(fro_track_ise), col = "green" )
segments(x0 = subsize, y0 = colMeans(fro_track_ise) - 1.96*se_fro_ise, 
	y1 = colMeans(fro_track_ise) + 1.96*se_fro_ise, col = "green")

lines(subsize,colMeans(fro_track_sve), col = "skyblue" )
segments(x0 = subsize, y0 = colMeans(fro_track_sve) - 1.96*se_fro_sve, 
	y1 = colMeans(fro_track_sve) + 1.96*se_fro_sve, col = "skyblue")

lines(subsize,colMeans(fro_track_mls), col = "brown" )
segments(x0 = subsize, y0 = colMeans(fro_track_mls) - 1.96*se_fro_mls, 
	y1 = colMeans(fro_track_mls) + 1.96*se_fro_mls, col = "brown")

abline(h = 1, lty = 2)
add_legend("topright", bty = "n",legend = c("BM", "New ISE (Geyer)", "ISE", "SVE", "New ISE (MLS)"), 
	col = c("black", "purple", "green", "skyblue", "brown"), lty = 1, cex=0.5)
dev.off()









