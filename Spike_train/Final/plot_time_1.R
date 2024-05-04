load("Time_d.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

bm_time = matrix(0, nrow = repet, ncol = length(N))
ise_time = matrix(0, nrow = repet, ncol = length(N))
sve_time = matrix(0, nrow = repet, ncol = length(N))
cc_time = matrix(0, nrow = repet, ncol = length(N))
mls_time = matrix(0, nrow = repet, ncol = length(N))


for(i in 1:repet){
	bm_time[i,] = as.numeric(Table[[i]][[2]][[1]])
	ise_time[i,] = as.numeric(Table[[i]][[2]][[2]])
	sve_time[i,] = as.numeric(Table[[i]][[2]][[3]])
	cc_time[i,] = as.numeric(Table[[i]][[2]][[4]])
	mls_time[i,] = as.numeric(Table[[i]][[2]][[5]])
}

se_time_bm <- apply(bm_time, 2, sd)/sqrt(repet)
se_time_ise <- apply(ise_time, 2, sd)/sqrt(repet)
se_time_cc <- apply(cc_time, 2, sd)/sqrt(repet)
se_time_sve <- apply(sve_time, 2, sd)/sqrt(repet)
se_time_mls <- apply(mls_time, 2, sd)/sqrt(repet)

#N = log(N)/log(10)


pdf("Calcium_spike_comptime.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(N, colMeans(bm_time),col = "black", xlab = "Chain Length", ylab = "Computational time (sec)", 
	ylim = c(0, 250), log = 'x', type = "l")
segments(x0 = N, y0 = colMeans(bm_time) - 1.96*se_time_bm, 
	y1 = colMeans(bm_time) + 1.96*se_time_bm)

lines(N, colMeans(sve_time), col = "skyblue", lty = 1)
segments(x0 = N, y0 = colMeans(sve_time) - 1.96*se_time_sve, 
	y1 = colMeans(sve_time) + 1.96*se_time_sve, col = "skyblue")

lines(N, colMeans(cc_time), col = "purple", lty = 1)
segments(x0 = N, y0 = colMeans(cc_time) - 1.96*se_time_cc, 
	y1 = colMeans(cc_time) + 1.96*se_time_cc, col = "purple")

lines(N, colMeans(ise_time), col = "red", lty = 1)
segments(x0 = N, y0 = colMeans(ise_time) - 1.96*se_time_ise, 
	y1 = colMeans(ise_time) + 1.96*se_time_ise, col = "red")

lines(N, colMeans(mls_time), col = "brown", lty = 1)
segments(x0 = N, y0 = colMeans(mls_time) - 1.96*se_time_mls, 
	y1 = colMeans(mls_time) + 1.96*se_time_mls, col = "brown")

legend("topleft",legend = c("BM", "SVE",  "ISE", "CC - ISE","CC - MLS"),
 col = c("black", "skyblue",  "red", "purple","brown"), cex = 0.8,lty = 1)
dev.off()

#####################################################################
###### ESS Plot #####################################################

load("Time_d.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

p = dim(Table[[1]][[1]][[1]])[1]

bm_ess = matrix(0, nrow = repet, ncol = length(N))
ise_ess = matrix(0, nrow = repet, ncol = length(N))
sve_ess = matrix(0, nrow = repet, ncol = length(N))
cc_ess = matrix(0, nrow = repet, ncol = length(N))
mls_ess = matrix(0, nrow = repet, ncol = length(N))

for(i in 1:repet){
	for(j in 1:length(N)){
		bm_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 2]]))^(1/p)
		ise_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 3]]))^(1/p)
		sve_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 4]]))^(1/p)
		cc_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 5]]))^(1/p)
		mls_ess[i,j] = (det(Table[[i]][[1]][[6*(j-1) + 1]])/det(Table[[i]][[1]][[6*(j-1) + 6]]))^(1/p)
	}
}


se_ess_bm <- apply(bm_ess, 2, sd)/sqrt(repet)
se_ess_sve <- apply(ise_ess, 2, sd)/sqrt(repet)
se_ess_cc <- apply(cc_ess, 2, sd)/sqrt(repet)
se_ess_ise <- apply(sve_ess, 2, sd)/sqrt(repet)
se_ess_mls <- apply(mls_ess, 2, sd)/sqrt(repet)

#N = log(N)/log(10)

pdf("Calcium_spike_ess.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(N, colMeans(bm_ess),col = "black", xlab = "Chain Length", ylab = "ESS/n", 
	ylim = c(0.05, 0.1), log = 'x', type = "l")
segments(x0 = N, y0 = colMeans(bm_ess) - 1.96*se_ess_bm, 
	y1 = colMeans(bm_ess) + 1.96*se_ess_bm)

lines(N, colMeans(sve_ess), col = "skyblue", lty = 1)
segments(x0 = N, y0 = colMeans(sve_ess) - 1.96*se_ess_sve, 
	y1 = colMeans(sve_ess) + 1.96*se_ess_sve, col = "skyblue")

lines(N, colMeans(cc_ess), col = "purple", lty = 1)
segments(x0 = N, y0 = colMeans(cc_ess) - 1.96*se_ess_cc, 
	y1 = colMeans(cc_ess) + 1.96*se_ess_cc, col = "purple")

lines(N, colMeans(ise_ess), col = "red", lty = 1)
segments(x0 = N, y0 = colMeans(ise_ess) - 1.96*se_ess_ise, 
	y1 = colMeans(ise_ess) + 1.96*se_ess_ise, col = "red")

lines(N, colMeans(mls_ess), col = "brown", lty = 1)
segments(x0 = N, y0 = colMeans(mls_ess) - 1.96*se_ess_mls, 
	y1 = colMeans(mls_ess) + 1.96*se_ess_mls, col = "brown")

#abline(h = 0, lty = 2)
legend("topright",legend = c("BM", "SVE", "ISE", "CC - ISE", "CC - MLS"),
 col = c("black", "skyblue", "red", "purple", "brown"), cex = 0.8,lty = 1)
dev.off()
