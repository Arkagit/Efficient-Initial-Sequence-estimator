load("bias_data.Rdata")

sizes <- c(1e3, 5e3, 1e4, 5e4)

# Saving bias plots for comparing estimation methods
pdf("plot_bias.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(sizes, b_bm, type = "b", col = "black", ylab = "Average bias on diagonals", xlab = "Chain length", ylim = c(-0.3, 0.2))
lines(sizes, b_wbm, col = "blue", type = "b", lty = 1)
lines(sizes, b_cc, col = "purple", type = "b", lty = 1)
lines(sizes, b_lug, col = "red", type = "b", lty = 1)
abline(h = 0, lty = 2)
legend("topright", bty = "n",legend = c("BM", "Zero Lugsail", "Cov-Corr", "Over Lugsail"), col = c("black", "blue", "purple", "red"), lty = 1)
dev.off()