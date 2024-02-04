load("Time_d.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

bm_time = log(as.numeric(Table[[1]][[1]]))
ise_time = log(as.numeric(Table[[1]][[2]]))
sve_time = log(as.numeric(Table[[1]][[3]]))
cc_time = log(as.numeric(Table[[1]][[4]]))
mls_time = log(as.numeric(Table[[1]][[5]]))

for(i in 2:repet){
	bm_time = rbind(bm_time, log(as.numeric(Table[[i]][[1]])))
	ise_time = rbind(ise_time, log(as.numeric(Table[[i]][[2]])))
	sve_time = rbind(sve_time, log(as.numeric(Table[[i]][[3]])))
	cc_time = rbind(cc_time, log(as.numeric(Table[[i]][[4]])))
	mls_time = rbind(cc_time, log(as.numeric(Table[[i]][[5]])))
}

se_time_bm <- apply(bm_time, 2, sd)/sqrt(repet)
se_time_ise <- apply(ise_time, 2, sd)/sqrt(repet)
se_time_sve <- apply(sve_time, 2, sd)/sqrt(repet)
se_time_cc <- apply(cc_time, 2, sd)/sqrt(repet)
se_time_mls <- apply(mls_time, 2, sd)/sqrt(repet)

# Saving bias plots for comparing estimation methods
pdf("plot_comptime_1.01.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(N, colMeans(bm_time),col = "black", xlab = "Chain Length", ylab = "Log computational time", 
	ylim = c(-5, 7), type = "l")
segments(x0 = N, y0 = colMeans(bm_time) - 1.96*se_time_bm, 
	y1 = colMeans(bm_time) + 1.96*se_time_bm)

lines(N, colMeans(sve_time), col = "skyblue", lty = 1)
segments(x0 = N, y0 = colMeans(sve_time) - 1.96*se_time_sve, 
	y1 = colMeans(sve_time) + 1.96*se_time_sve, col = "skyblue")

lines(N, colMeans(cc_time), col = "purple", lty = 1, lwd = 2)
segments(x0 = N, y0 = colMeans(cc_time) - 1.96*se_time_cc, 
	y1 = colMeans(cc_time) + 1.96*se_time_cc, col = "purple", lwd = 2)

lines(N, colMeans(ise_time), col = "green", lty = 1)
segments(x0 = N, y0 = colMeans(ise_time) - 1.96*se_time_ise, 
	y1 = colMeans(ise_time) + 1.96*se_time_ise, col = "green")

lines(N, colMeans(mls_time), col = "brown", lty = 3)
segments(x0 = N, y0 = colMeans(mls_time) - 1.96*se_time_mls, 
	y1 = colMeans(mls_time) + 1.96*se_time_mls, col = "brown")

abline(h = 0, lty = 2)
legend("bottomright",legend = c("BM", "SVE", "New ISE (Geyer)", "ISE", "MLS"),
 col = c("black", "skyblue", "purple", "green", "brown"), cex = 0.8, lty = 1, lwd = c(1,1,2,1,3))
dev.off()