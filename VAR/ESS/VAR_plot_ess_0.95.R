source("VAR_ess.R")
load("ESS_data_0.95.Rdata")

p <- 10
rho <- 0.95
phi <- diag(rep(rho,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))

truth <- true.sig(p = p, omega = omega, phi = phi)

se_ess_bm <- apply(ess_track$BM, 2, sd)/sqrt(100)
se_ess_lug <- apply(ess_track$LUG, 2, sd)/sqrt(100)
se_ess_ise <- apply(ess_track$ISE, 2, sd)/sqrt(100)
se_ess_cc <- apply(ess_track$CC, 2, sd)/sqrt(100)

se_fro_bm <- apply(ess_track$Frob_bm, 2, sd)/sqrt(100)
se_fro_lug <- apply(ess_track$Frob_lug, 2, sd)/sqrt(100)
se_fro_ise <- apply(ess_track$Frob_ise, 2, sd)/sqrt(100)
se_fro_cc <- apply(ess_track$Frob_cc, 2, sd)/sqrt(100)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


pdf("plot_ess_0.95.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track$BM), ylim = c(0.024, 0.040),
	type = 'l', xlab = "Time", ylab = expression(paste("hat(ESS)/n (", rho, " = 0.05)")))
segments(x0 = subsize, y0 = colMeans(ess_track$BM) - 1.96*se_ess_bm, y1 = colMeans(ess_track$BM) + 1.96*se_ess_bm)

lines(subsize,colMeans(ess_track$LUG), col = "red" )
segments(x0 = subsize, y0 = colMeans(ess_track$LUG) - 1.96*se_ess_lug, y1 = colMeans(ess_track$LUG) + 1.96*se_ess_lug, col = "red")

lines(subsize,colMeans(ess_track$CC), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track$CC) - 1.96*se_ess_cc, y1 = colMeans(ess_track$CC) + 1.96*se_ess_cc, col = "purple")

lines(subsize,colMeans(ess_track$ISE), col = "green" )
segments(x0 = subsize, y0 = colMeans(ess_track$ISE) - 1.96*se_ess_ise, y1 = colMeans(ess_track$ISE) + 1.96*se_ess_ise, col = "green")

abline(h = ess_true, lty = 2)
add_legend("topright", bty = "n",legend = c("BM", "CC", "Over Lugsail", "ISE"), col = c("black", "purple", "red","green"), lty = 1, horiz=TRUE, cex=0.7)
dev.off()



pdf("plot_frob_0.95.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track$Frob_bm),type = 'l', xlab = "Time", ylim = c(2300, 4300),
	ylab = expression(paste("Frobenius Norm (", rho, " = 0.05)")))
segments(x0 = subsize, y0 = colMeans(ess_track$Frob_bm) - 1.96*se_fro_bm, y1 = colMeans(ess_track$Frob_bm) + 1.96*se_fro_bm)

lines(subsize,colMeans(ess_track$Frob_lug), col = "red" )
segments(x0 = subsize, y0 = colMeans(ess_track$Frob_lug) - 1.96*se_fro_lug, y1 = colMeans(ess_track$Frob_lug) + 1.96*se_fro_lug, col = "red")

lines(subsize,colMeans(ess_track$Frob_cc), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track$Frob_cc) - 1.96*se_fro_cc, y1 = colMeans(ess_track$Frob_cc) + 1.96*se_fro_cc, col = "purple")

lines(subsize,colMeans(ess_track$Frob_ise), col = "green" )
segments(x0 = subsize, y0 = colMeans(ess_track$Frob_ise) - 1.96*se_fro_ise, y1 = colMeans(ess_track$Frob_ise) + 1.96*se_fro_ise, col = "green")

abline(h = norm(truth$final.cov, type = "F"), lty = 2)
add_legend("topright", bty = "n",legend = c("BM", "CC", "Over Lugsail", "ISE"), col = c("black", "purple", "red","green"), lty = 1, horiz=TRUE, cex=0.7)
dev.off()