load("dat_matices.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

time = log(time)

pdf("coverage_time1.pdf", height = 6, width = 8)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(chart[,1], time[,1], xlim = c(0.1,1), ylim = c(0,9.3), type = 'b', xlab = "Coverage",
	ylab = "log(Time)", main = "Log Computational Time vs Coverage Probability", pch = 1:dim(time)[1])
#lines(chart[,2], avg_time[,2],type = 'b', col = "red", pch = 1:length(subsize))
lines(chart[,3], time[,3],type = 'b', col = "green", pch = 1:dim(time)[1])
lines(chart[,4], time[,4],type = 'b', col = "purple", pch = 1:dim(time)[1])
lines(chart[,5], time[,5],type = 'b', col = "skyblue", pch = 1:dim(time)[1])
lines(chart[,6], time[,6],type = 'b', col = "brown", pch = 1:dim(time)[1])
add_legend("topright", bty = "n",legend = c("BM", "New ISE (Geyer)", "ISE", "SVE", "New ISE (MLS)"), 
	col = c("black", "purple", "green", "skyblue", "brown"), lty = 1, cex=0.5)

dev.off()