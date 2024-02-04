load("Time.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

bm_time = log(as.numeric(Table[[1]][[1]]))
lug_time = log(as.numeric(Table[[1]][[2]]))
ise_time = log(as.numeric(Table[[1]][[3]]))
sve_time = log(as.numeric(Table[[1]][[4]]))
cc_time = log(as.numeric(Table[[1]][[5]]))
mls_time = log(as.numeric(Table[[1]][[6]]))
for(i in 2:1000){
	bm_time = append(bm_time, log(as.numeric(Table[[i]][[1]])))
	lug_time = append(lug_time, log(as.numeric(Table[[i]][[2]])))
	ise_time = append(ise_time, log(as.numeric(Table[[i]][[3]])))
	sve_time = append(sve_time, log(as.numeric(Table[[i]][[4]])))
	cc_time = append(cc_time, log(as.numeric(Table[[i]][[5]])))
	mls_time = append(mls_time, log(as.numeric(Table[[i]][[6]])))
}



# Saving bias plots for comparing estimation methods
pdf("plot_time.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(density(bm_time),col = "black", ylab = "density", xlab = "logarithmic computational time", main = "Density plot of computational time", 
	xlim = c(-4, 4), ylim = c(0, 4.5))
lines(density(sve_time), col = "skyblue", lty = 1)
lines(density(cc_time), col = "purple", lty = 1)
lines(density(lug_time), col = "red", lty = 1)
lines(density(ise_time), col = "green", lty = 1)
lines(density(mls_time), col = "brown", lty = 1)
abline(h = 0, lty = 2)
legend("topright",legend = c("BM", "SVE", "New ISE (Geyer)", "Over Lugsail", "ISE", "MLS"),
 col = c("black", "skyblue", "purple", "red", "green", "brown"), cex=0.50,lty = 1)
dev.off()