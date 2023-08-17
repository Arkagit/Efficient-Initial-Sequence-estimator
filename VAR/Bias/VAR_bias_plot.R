load("bias_data.Rdata")

sizes <- c(1e3, 5e3, 1e4, 5e4)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

# Saving bias plots for comparing estimation methods
pdf("plot_bias.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(sizes, b_bm, type = "b", col = "black", ylab = "Average bias on diagonals", xlab = "Chain length",
 ylim = c(-0.3, 0.2))
lines(sizes, b_sve, col = "skyblue", type = "b", lty = 1)
lines(sizes, b_cc, col = "purple", type = "b", lty = 1)
lines(sizes, b_lug, col = "red", type = "b", lty = 1)
lines(sizes, b_ise, col = "green", type = "b", lty = 1)
lines(sizes, b_mls, col = "brown", type = "b", lty = 1)
abline(h = 0, lty = 2)
add_legend("top", inset=c(-0.3, 0), bty = "n",legend = c("BM", "SVE", "New ISE (Geyer)", "Over Lugsail", "ISE", "New ISE (MLS)"),
 col = c("black", "skyblue", "purple", "red", "green", "brown"), cex=0.50,lty = 1)
dev.off()