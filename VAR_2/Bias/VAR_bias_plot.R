load("bias_data.Rdata")

sizes <- subsize

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
b_bm = c(BS[[1]][[1]], BS[[2]][[1]], BS[[3]][[1]],BS[[4]][[1]],BS[[5]][[1]])
b_ise = c(BS[[1]][[2]], BS[[2]][[2]], BS[[3]][[2]], BS[[4]][[2]], BS[[5]][[2]])
b_cc = c(BS[[1]][[3]], BS[[2]][[3]], BS[[3]][[3]], BS[[4]][[3]], BS[[5]][[3]])
b_sve = c(BS[[1]][[4]], BS[[2]][[4]], BS[[3]][[4]], BS[[4]][[4]], BS[[5]][[4]])
b_mls = c(BS[[1]][[5]], BS[[2]][[5]], BS[[3]][[5]], BS[[4]][[5]], BS[[5]][[5]])
# Saving bias plots for comparing estimation methods
pdf("plot_bias_2.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(sizes, b_bm, type = "b", col = "black", ylab = "Average bias of determinants", xlab = "Chain length",
 ylim = c(-0.05, 0.05))
lines(sizes, b_sve, col = "skyblue", type = "b", lty = 1)
lines(sizes, b_cc, col = "purple", type = "b", lty = 1)
lines(sizes, b_ise, col = "green", type = "b", lty = 1)
lines(sizes, b_mls, col = "brown", type = "b", lty = 1)
abline(h = 0, lty = 2)
add_legend("top", inset=c(-0.3, 0), bty = "n",legend = c("BM", "SVE", "New ISE (Geyer)", "New ISE (MLS)", "ISE"),
 col = c("black", "skyblue", "purple", "brown", "green"), cex=0.50,lty = 1)
dev.off()