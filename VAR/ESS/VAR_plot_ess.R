source("VAR_ess.R")
load("ESS_data_0.99.Rdata")


# Initialization
p <- 10
rho <- 0.99
phi <- diag(rep(rho,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))

truth <- true.sig(p = p, omega = omega, phi = phi)

se_ess_bm <- apply(ess_track_bm, 2, sd)/sqrt(nloops)
se_ess_lug <- apply(ess_track_lug, 2, sd)/sqrt(nloops)
se_ess_ise <- apply(ess_track_ise, 2, sd)/sqrt(nloops)
se_ess_cc <- apply(ess_track_cc, 2, sd)/sqrt(nloops)
se_ess_sve <- apply(ess_track_sve, 2, sd)/sqrt(nloops)
se_ess_mls <- apply(ess_track_mls, 2, sd)/sqrt(nloops)

se_fro_bm <- apply(fro_track_bm, 2, sd)/sqrt(nloops)
se_fro_lug <- apply(fro_track_lug, 2, sd)/sqrt(nloops)
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
pdf("plot_ess_0.99.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track_bm), ylim = c(0, 0.012),
	type = 'l', xlab = "Time", ylab = expression(paste(hat(ESS)/n," (", 
		rho, " = 0.99)")))
segments(x0 = subsize, y0 = colMeans(ess_track_bm) - 1.96*se_ess_bm, 
	y1 = colMeans(ess_track_bm) + 1.96*se_ess_bm)

lines(subsize,colMeans(ess_track_lug), col = "red" )
segments(x0 = subsize, y0 = colMeans(ess_track_lug) - 1.96*se_ess_lug, 
	y1 = colMeans(ess_track_lug) + 1.96*se_ess_lug, col = "red")

lines(subsize,colMeans(ess_track_cc), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track_cc) - 1.96*se_ess_cc, 
	y1 = colMeans(ess_track_cc) + 1.96*se_ess_cc, col = "purple")

lines(subsize,colMeans(ess_track_ise), col = "green" )
segments(x0 = subsize, y0 = colMeans(ess_track_ise) - 1.96*se_ess_ise, 
	y1 = colMeans(ess_track_ise) + 1.96*se_ess_ise, col = "green")

lines(subsize,colMeans(ess_track_sve), col = "skyblue" )
segments(x0 = subsize, y0 = colMeans(ess_track_sve) - 1.96*se_ess_sve, 
	y1 = colMeans(ess_track_sve) + 1.96*se_ess_sve, col = "skyblue")

lines(subsize,colMeans(ess_track_mls), col = "brown" )
segments(x0 = subsize, y0 = colMeans(ess_track_mls) - 1.96*se_ess_mls, 
	y1 = colMeans(ess_track_mls) + 1.96*se_ess_mls, col = "brown")

abline(h = ess_true, lty = 2)
add_legend("topright", bty = "n",legend = c("BM", "New ISE (Geyer)", "Over Lugsail", "ISE", "SVE", "New ISE (MLS)"), 
	col = c("black", "purple", "red","green", "skyblue", "brown"), lty = 1, cex=0.5)
dev.off()


# Frobenius norm plot
pdf("plot_frob_0.99.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(fro_track_bm),type = 'l', xlab = "Time", 
	ylim = c(0.7, 1.5), ylab = "Frobenius Norm")
segments(x0 = subsize, y0 = colMeans(fro_track_bm) - 1.96*se_fro_bm, 
	y1 = colMeans(fro_track_bm) + 1.96*se_fro_bm)

lines(subsize,colMeans(fro_track_lug), col = "red" )
segments(x0 = subsize, y0 = colMeans(fro_track_lug) - 1.96*se_fro_lug, 
	y1 = colMeans(fro_track_lug) + 1.96*se_fro_lug, col = "red")

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
add_legend("topright", bty = "n",legend = c("BM", "New ISE (Geyer)", "Over Lugsail", "ISE", "SVE", "New ISE (MLS)"), 
	col = c("black", "purple", "red","green", "skyblue", "brown"), lty = 1, cex=0.5)
dev.off()









