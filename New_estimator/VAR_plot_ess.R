load("ESS_data.Rdata")

se_ess_bm <- apply(ess_track$BM, 2, sd)/sqrt(100)
se_ess_wbm <- apply(ess_track$WBM, 2, sd)/sqrt(100)
se_ess_lug <- apply(ess_track$LUG, 2, sd)/sqrt(100)
se_ess_ise <- apply(ess_track$ISE, 2, sd)/sqrt(100)
se_ess_cc <- apply(ess_track$CC, 2, sd)/sqrt(100)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


pdf("plot_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track$BM), ylim = c(.024, .040), 
	type = 'l', xlab = "Time", ylab = expression(hat(ESS)/n))
segments(x0 = subsize, y0 = colMeans(ess_track$BM) - 1.96*se_ess_bm, y1 = colMeans(ess_track$BM) + 1.96*se_ess_bm)

lines(subsize,colMeans(ess_track$WBM), col = "blue" )
segments(x0 = subsize, y0 = colMeans(ess_track$WBM) - 1.96*se_ess_wbm, y1 = colMeans(ess_track$WBM) + 1.96*se_ess_wbm, col = "blue")

lines(subsize,colMeans(ess_track$LUG), col = "red" )
segments(x0 = subsize, y0 = colMeans(ess_track$LUG) - 1.96*se_ess_lug, y1 = colMeans(ess_track$LUG) + 1.96*se_ess_lug, col = "red")

lines(subsize,colMeans(ess_track$CC), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track$CC) - 1.96*se_ess_cc, y1 = colMeans(ess_track$CC) + 1.96*se_ess_cc, col = "purple")

lines(subsize,colMeans(ess_track$ISE), col = "green" )
segments(x0 = subsize, y0 = colMeans(ess_track$ISE) - 1.96*se_ess_ise, y1 = colMeans(ess_track$ISE) + 1.96*se_ess_ise, col = "green")

abline(h = ess_true, lty = 2)
add_legend("topright", bty = "n",legend = c("BM", "Zero Lugsail", "Adapt Lugsail", "Over Lugsail", "ISE"), col = c("black", "blue", "purple", "red","green"), lty = 1, horiz=TRUE, cex=0.7)
dev.off()