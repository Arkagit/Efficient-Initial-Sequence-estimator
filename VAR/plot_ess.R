load("ESS_data.Rdata")

se_ess_bm <- apply(ess_track$BM, 2, sd)/sqrt(nloops)
se_ess_wbm <- apply(ess_track$WBM, 2, sd)/sqrt(nloops)
se_ess_lug <- apply(ess_track$LUG, 2, sd)/sqrt(nloops)
se_ess_cc <- apply(ess_track$CC, 2, sd)/sqrt(nloops)
se_ess_ise <- apply(ess_track$ISE, 2, sd)/sqrt(nloops)

pdf("plot_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

plot(subsize,colMeans(ess_track$ISE), col = "green", type = 'l', ylim = c(0.023, 0.04), xlab = "Time", ylab = expression(hat(ESS)/n))
segments(x0 = subsize, y0 = colMeans(ess_track$ISE) - 1.96*se_ess_ise, y1 = colMeans(ess_track$ISE) + 1.96*se_ess_ise, col = "green")

lines(subsize, colMeans(ess_track$BM))
segments(x0 = subsize, y0 = colMeans(ess_track$BM) - 1.96*se_ess_bm, y1 = colMeans(ess_track$BM) + 1.96*se_ess_bm)

lines(subsize,colMeans(ess_track$WBM), col = "blue" )
segments(x0 = subsize, y0 = colMeans(ess_track$WBM) - 1.96*se_ess_wbm, y1 = colMeans(ess_track$WBM) + 1.96*se_ess_wbm, col = "blue")

lines(subsize,colMeans(ess_track$LUG), col = "red" )
segments(x0 = subsize, y0 = colMeans(ess_track$LUG) - 1.96*se_ess_lug, y1 = colMeans(ess_track$LUG) + 1.96*se_ess_lug, col = "red")

lines(subsize,colMeans(ess_track$CC), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track$CC) - 1.96*se_ess_cc, y1 = colMeans(ess_track$CC) + 1.96*se_ess_cc, col = "purple")


abline(h = ess_true, lty = 2)
legend("topright", bty = "n",legend = c("BM", "Zero Lugsail", "CC", "Over Lugsail", "ISE"), col = c("black", "blue", "purple", "red", "green"), lty = 1)

dev.off()

