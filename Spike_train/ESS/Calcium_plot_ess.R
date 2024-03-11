source("VAR_ess.R")
load("ESS_Frob_data.Rdata")
Table


ess_track_bm = matrix(0, nrow = B, ncol = nloops)
ess_track_ise = matrix(0, nrow = B, ncol = nloops)
ess_track_cc = matrix(0, nrow = B, ncol = nloops)
ess_track_sve = matrix(0, nrow = B, ncol = nloops)

time_track_bm = matrix(0, nrow = B, ncol = nloops)
time_track_ise = matrix(0, nrow = B, ncol = nloops)
time_track_cc = matrix(0, nrow = B, ncol = nloops)
time_track_sve = matrix(0, nrow = B, ncol = nloops)


for(i in 1:B){
	for(j in 1:nloops){
		ess_track_bm[i,j] = as.numeric(Table[[i]][[j]][[(j-1)+1]][[1]])
		ess_track_ise[i,j] = as.numeric(Table[[i]][[j]][[(j-1)+1]][[2]])
		ess_track_cc[i,j] = as.numeric(Table[[i]][[j]][[(j-1)+1]][[3]])
		ess_track_sve[i,j] = as.numeric(Table[[i]][[j]][[(j-1)+1]][[4]])

		time_track_bm[i,j] = as.numeric(Table[[i]][[j]][[(j-1)+2]][[1]])
		time_track_ise[i,j] = as.numeric(Table[[i]][[j]][[(j-1)+2]][[2]])
		time_track_cc[i,j] = as.numeric(Table[[i]][[j]][[(j-1)+2]][[3]])
		time_track_sve[i,j] = as.numeric(Table[[i]][[j]][[(j-1)+2]][[4]])
	}
}


se_ess_bm <- apply(ess_track_bm, 2, sd)/sqrt(B)
se_ess_ise <- apply(ess_track_ise, 2, sd)/sqrt(B)
se_ess_cc <- apply(ess_track_cc, 2, sd)/sqrt(B)
se_ess_sve <- apply(ess_track_sve, 2, sd)/sqrt(B)

se_time_bm <- apply(time_track_bm, 2, sd)/sqrt(B)
se_time_ise <- apply(time_track_ise, 2, sd)/sqrt(B)
se_time_cc <- apply(time_track_cc, 2, sd)/sqrt(B)
se_time_sve <- apply(time_track_sve, 2, sd)/sqrt(B)


add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


# ESS plot
pdf("plot_ess_calc.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track_bm), ylim = c(0.05, 0.1),
	type = 'l', xlab = "Time", ylab = expression(paste(hat(ESS)/n)))
segments(x0 = subsize, y0 = colMeans(ess_track_bm) - 1.96*se_ess_bm, 
	y1 = colMeans(ess_track_bm) + 1.96*se_ess_bm)

lines(subsize, colMeans(ess_track_cc), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track_cc) - 1.96*se_ess_cc, 
	y1 = colMeans(ess_track_cc) + 1.96*se_ess_cc, col = "purple")

lines(subsize, colMeans(ess_track_ise), col = "green" )
segments(x0 = subsize, y0 = colMeans(ess_track_ise) - 1.96*se_ess_ise, 
	y1 = colMeans(ess_track_ise) + 1.96*se_ess_ise, col = "green")

lines(subsize,colMeans(ess_track_sve), col = "skyblue" )
segments(x0 = subsize, y0 = colMeans(ess_track_sve) - 1.96*se_ess_sve, 
	y1 = colMeans(ess_track_sve) + 1.96*se_ess_sve, col = "skyblue")

#abline(h = ess_true, lty = 2)
legend("topright", bty = "n",legend = c("BM", "New ISE (Geyer)", "ISE", "SVE"), 
	col = c("black", "purple","green", "skyblue"), lty = 1, cex=0.8)
dev.off()












