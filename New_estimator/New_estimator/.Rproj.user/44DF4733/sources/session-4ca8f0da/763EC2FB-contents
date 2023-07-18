set.seed(100)


########### Lugsail lag window ############
b <- 100
r <- 2
x <- (0:(2*b))/b
b <- 1

bart <- function(x, b)
{
  (1 - x/b)*(x <= b)
}

th <- function(x, b)
{
 ( (1 + cos(pi*x/b) )/2 )*(x <= b)
}

qs <- function(x, b)
{
  x <- x/b
  25/(12*pi^2*x^2) * (sin(6*pi*x/5)/(6*pi*x/5) - cos(6*pi*x/5)   )
}

lugsail <- function(x, n, r = 2, c = 1/2, b = 1, change = TRUE, window = "bart")
{
	if(window == "bart") 
		{
			lagwin <- bart
			q <- 1
		}
	if(window == "th") 
		{
			lagwin <- th
			q <- 2
		}
	if(window == "qs") 
		{
			lagwin <- qs
			q <- 2
		}

	a <- sqrt(n)
	if(change) c <- (log(a)+1)/(log(a)*r^q + 1)
		print(c)
	window <- (1/(1-c))*(lagwin(x,b)) - (c/(1-c))*(lagwin(x,b/r))
}


lug.old <-  lugsail(x,n = 1e2, r = 3, c = 1/2, change = FALSE)
lug.zero <-  lugsail(x,n = 1e2, r = 2, c = 1/2, change = FALSE)
lug2 <- lugsail(x,n = 1e2, r = 2)
lug3 <- lugsail(x,n = 1e3, r = 2)
lug4 <- lugsail(x,n = 1e4, r = 2)
lug5 <- lugsail(x, n = 1e5, r = 2)
lug6 <- lugsail(x, n = 1e6, r = 2)


####################################
### Figure 1: Lag windows
####################################

pdf("plots/lag_window_bart.pdf", width = 5, height = 5.5)
plot(x, bart(x,b), xlim = c(0,1.2),ylim = c(-.12,1.3), type = 'l', main = expression("Bartlett"), ylab = "", xlab = "k")
lines(x, lug.zero, lty = 1, col = "blue")
lines(x, lug2, lty = 2, col = "blue")
lines(x, lug3, lty = 3, col = "blue")
lines(x, lug4, lty = 4, col = "blue")
lines(x, lug5, lty = 5, col = "blue")
lines(x, lug.old, col = "red")
legend("topright", bty = "n", legend = c("Over Lugsail", "Adapt Lugsail", "Zero Lugsail", "Original Bartlett"), lty = c(1,2,1,1), col = c("red", "blue", "blue", "black"))
dev.off()


lug.old <-  lugsail(x,n = 1e2, r = 3, c = 1/5, change = FALSE, window = "th")
lug.zero <-  lugsail(x,n = 1e2, r = 2, c = 1/4, change = FALSE, window = "th")
lug2 <- lugsail(x,n = 1e2, r = 2, window = "th")
lug3 <- lugsail(x,n = 1e3, r = 2, window = "th")
lug4 <- lugsail(x,n = 1e4, r = 2, window = "th")
lug5 <- lugsail(x, n = 1e5, r = 2, window = "th")


pdf("plots/lag_window_th.pdf", width = 5, height = 5.5)
plot(x, th(x,b), xlim = c(0,1.2),ylim = c(-.12,1.3),  type = 'l', ylab = "", main = expression("Tukey-Hanning"), xlab = "k")
lines(x, lug.zero, lty = 1, col = "blue")
lines(x, lug2, lty = 2, col = "blue")
lines(x, lug3, lty = 3, col = "blue")
lines(x, lug4, lty = 4, col = "blue")
lines(x, lug5, lty = 5, col = "blue")
lines(x, lug.old, col = "red")
legend("topright", bty = "n",legend = c("Over Lugsail", "Adapt Lugsail", "Zero Lugsail", "Original TH"), lty = c(1,2,1,1), col = c("red", "blue", "blue", "black"))
dev.off()




lug.old <-  lugsail(x,n = 1e2, r = 3, c = 1/5, change = FALSE, window = "qs")
lug.zero <-  lugsail(x,n = 1e2, r = 2, c = 1/4, change = FALSE, window = "qs")
lug2 <- lugsail(x, n = 1e2, r = 2, window = "qs")
lug3 <- lugsail(x, n = 1e3, r = 2, window = "qs")
lug4 <- lugsail(x, n = 1e4, r = 2, window = "qs")
lug5 <- lugsail(x, n = 1e5, r = 2, window = "qs")


pdf("plots/lag_window_qs.pdf", width = 5, height = 5.5)
plot(x, qs(x,b), xlim = c(0,1.5),ylim = c(-.12,1.3),  type = 'l', ylab = "", main =  expression("QS"), xlab = "k")
lines(x, lug.zero, lty = 1, col = "blue")
lines(x, lug2, lty = 2, col = "blue")
lines(x, lug3, lty = 3, col = "blue")
lines(x, lug4, lty = 4, col = "blue")
lines(x, lug5, lty = 5, col = "blue")
lines(x, lug.old, col = "red")
legend("topright", bty = "n", legend = c("Over Lugsail", "Adapt Lugsail", "Zero Lugsail", "Original QS"), lty = c(1,2,1,1), col = c("red", "blue", "blue", "black"))
dev.off()


####################################
### Figure 2: True Bias plots
####################################

library(ts.extend)

# AR(p) process bias for batch means estimators
arp_bias <- function(a, b, phi, sigma2 = 1)
{
  n <- a*b
  p <- length(phi)
  gamma <- (sigma2)*ARMA.autocov(n = n, ar = phi, ma = 0,corr = FALSE)
  
  foo <- gamma[1]
  foo <- foo + 2*sum(gamma[2:b]) - 2*(a + 1)/(a*b) * sum( (1:(b-1))*gamma[2:b] ) - 2/(a - 1) * sum(  (1 - (b:(n-1))/n )*gamma[(b+1):n] )
  
  foo <- foo - sigma2/(1 - sum(phi))^2
  return(foo)
}

bias <- function(phi, a, b, r, alpha)
{
	rtn <- alpha*arp_bias(a = a, b = b, phi = phi) + (1-alpha)*arp_bias(a*r, b/r, phi)
	return(rtn)
}

truth <- function(phi) 1/(1- sum(phi))^2

plot_bias_bm <- function(phi.vec, N, comp, main = "")
{
	K <- length(N)
	bias.vec3 <- bias.vec2 <- bias.vec1 <- bias.vecJ1 <- numeric(length = K)
	bias.vecJ2 <- bias.vecJ3 <-  bias.vecJ4 <- numeric(length = K)
	c <- 1/2 
	alpha <- 1/(1-c)

	for(k in 1:K)
	{
		n <- N[k]
		b <- floor(n^(1/2))
		a <- floor(n/b)
		alpha <- 1/(1-c)	  
		
		bias.vec3[k] <- bias(phi.vec, a = a, b = b, r = 3,  alpha = alpha)/truth(phi.vec)
		bias.vec2[k] <- bias(phi.vec, a = a, b = b, r = 2,  alpha = alpha)/truth(phi.vec)
		
		bias.vec1[k] <- bias(phi.vec, a = a, b = b, r = 1,  alpha = alpha)/truth(phi.vec)

		r.J <- 1/2
		r.here <- 2
		alpha.J <- ( a + r.J)/(a*(1 - r.J))
		bias.vecJ1[k] <- bias(phi.vec, a = a, b = b, r = r.here, alpha = alpha.J)/truth(phi.vec)
		
		alpha.J <- ( log(a) + r.J)/(log(a) * (1 - r.J))
		bias.vecJ4[k] <- bias(phi.vec, a = a, b = b, r = r.here, alpha = alpha.J)/truth(phi.vec)
	}


	plot(N, bias.vec2, type = 'l', 
		log = "x", xaxt = "n", xlab = "n (in log scale)", ylab = "Relative Bias", main = main, font.main = 1, col = "blue", ylim = c(-.60, .08))
	axis(1,  at = N, labels = N, las=1)
	lines(N, bias.vec3, col = "red")
	lines(N, bias.vec1, col = "black", lty = 1)
	lines(N, bias.vecJ4, col = "blue", lty = 2)

	abline(h = 0, lty = 3)
	legend("bottomright", bty = "n", legend = c("Over Lugsail",  "Adapt Lugsail" ,"Zero Lugsail",  "Original BM"), col = c("red", "blue", "blue",  "black"), lty = c(1,2,1,1))
}



## AR(p) approximation of the time-series example
library(mvtnorm)
T <- 1e4
u.rho <- c(.50, .90, .98)
x.rho <- c(.50, .90, .98)
W.rho <- .995
p <- 4
foo <- matrix(W.rho, nrow= p, ncol=p)
W <- foo^(abs(col(foo)-row(foo)))
# W <- diag(1-x.rho^2, p)
beta.star <- rep(0, p)
w <- 1

main.title <- c("Moderate Correlation", "High Correlation", "Extreme Correlation")
for(j in 1:length(u.rho))
{
	X <- matrix(0, nrow  = T, ncol = p)
	X[1,] <- rnorm(p)
	for(t in 2:T)
	{
		X[t, ] <- x.rho[j] * X[t-1, ] + MASS::mvrnorm(1, mu = rep(0, p), Sigma = W)
	}

	u <- numeric(length = T)
	u[1] <- 0

	for(t in 2:T)
	{
		u[t] <- u.rho[j]*u[t-1] + rnorm(1, sd = sqrt(w))
	}

	y <- X%*%beta.star + u
	XtX <- t(X)%*%X
	beta.ols <- solve(XtX) %*% t(X)%*%y
	u.hat <- y - X%*%beta.ols
	v.hat <- apply(X, 2, function(x) x*u.hat)

	out <- v.hat
	arp <- numeric(length = dim(out)[2])
	for(i in 1:dim(out)[2])
	{
		arp[i] <- ar(out[,i])$order
	}

	N <- c(1e2, 1e3, 1e4, 1e5)
	phi <- ar(out[,1])$ar

	print(phi)[1]
	pdf(paste("plots/BM",j,".pdf",sep = ""), height = 6, width = 5)	
	plot_bias_bm(phi.vec = phi, N = N, comp = 1, main = main.title[j])
	dev.off()

}




####################################
### Examples
################################

### Table 2: Time-Series -- coverage

rm(list = ls())
true.var <- function(u.rho = .70, x.rho = .90, p = 5, W.rho = .50, w = 1)
{
	foo <- matrix(W.rho, nrow= p, ncol=p)
	W <- foo^(abs(col(foo)-row(foo)))

	lam <- w/(1 - u.rho^2)
	Lam <- W/(1 - x.rho^2)

	truth <- lam*Lam + (Lam + t(Lam))*lam*(u.rho*x.rho/(1 - u.rho*x.rho))

	return(truth)
}

T.vec <- c(100, 500, 1000)

u.rho <- c(.50, .70, .90)
x.rho <- c(.50, .70, .90)
W.rho <- .99
p <- 5
foo <- matrix(W.rho, nrow= p, ncol=p)
W <- foo^(abs(col(foo)-row(foo)))
beta.star <- rep(0, p)
w <- 1

truth.mat <- matrix(0, nrow = 9, ncol = 5)
k <- 1
for(i in 1:length(T.vec))
{
	for(j in 1:length(u.rho))
	{
		truth.mat[k,] <- diag(true.var(u.rho = u.rho[j], x.rho = x.rho[j], W.rho = W.rho, w = 1))
		k <- k + 1
	}
}

load("Examples/TimeSeries/bias_ts_other")
load("Examples/TimeSeries/mse_ts_other")
load("Examples/TimeSeries/coverage_ts_other")

cover_ts <- t(rbind(colMeans(coverage_bt), colMeans(coverage_bt2), colMeans(coverage_bt5), colMeans(coverage_bt3), colMeans(coverage_th), colMeans(coverage_th2), colMeans(coverage_th5), colMeans(coverage_th3), colMeans(coverage_qs), colMeans(coverage_qs2), colMeans(coverage_qs5), colMeans(coverage_qs3), colMeans(coverage_truth)))
colnames(cover_ts) <- c("BT", "BT2",  "BTlog","BT3", "TH", "TH2", "THlog", "TH3", "QS", "QS2", "QSlog", "QS3", "truth")
rownames(cover_ts) <- rep(c(".50", ".70", ".90"), 3)
cover_ts[-(1:3), ]


cover_se <- round(sqrt(cover_ts*(1-cover_ts)/1000), 4)
rownames(cover_se) <- rep(c(".50", ".70", ".90"), 3)
max(cover_se)




### Figure 3: Time-Series -- Bias
dum <- function(x)
{
	dum1 <- function(y)
	{
		foo <- matrix(y, ncol = 5, nrow = 5, byrow = TRUE)
		return(diag(foo))
	}
	rtn <- apply(x, 1, dum1)
	rtn <- rowMeans(rtn)
	return(rtn)  
}

bt <-  colMeans(sapply(bias_bt, dum)) /truth.mat[,1] 
bt2 <- colMeans(sapply(bias_bt2, dum))/truth.mat[,1] 
bt3 <- colMeans(sapply(bias_bt3, dum))/truth.mat[,1] 
bt4 <- colMeans(sapply(bias_bt4, dum))/truth.mat[,1]
bt5 <- colMeans(sapply(bias_bt5, dum))/truth.mat[,1]

bth <- colMeans(sapply(bias_th, dum)) /truth.mat[,1]
bth2 <- colMeans(sapply(bias_th2, dum))/truth.mat[,1] 
bth3 <- colMeans(sapply(bias_th3, dum))/truth.mat[,1] 
bth4 <- colMeans(sapply(bias_th4, dum))/truth.mat[,1]
bth5 <- colMeans(sapply(bias_th5, dum))/truth.mat[,1]

bqs <- colMeans(sapply(bias_qs, dum))/truth.mat[,1] 
bqs2 <- colMeans(sapply(bias_qs2, dum)) /truth.mat[,1]
bqs3 <- colMeans(sapply(bias_qs3, dum))/truth.mat[,1] 
bqs4 <- colMeans(sapply(bias_qs4, dum))/truth.mat[,1]
bqs5 <- colMeans(sapply(bias_qs5, dum))/truth.mat[,1]

bias_est <-   round(cbind(bt, bt2, bt4, bt5, bt3, bth, bth2, bth4, bth5, bth3, bqs, bqs2, bqs4, bqs5, bqs3), 5)

bias_table <- bias_est
rownames(bias_table) <- rep(c(".50", ".70", ".90"), 3)
bias_table


bias_bart <- bias_table[ , 1:5 ] # Barltlett
bias_tukey <- bias_table[ , 6:10 ] # Tukey
bias_qs <- bias_table[ , 11:15 ] # Quadratic Spectral 


pdf("plots/hac_bias_n1000.pdf", height = 4, width = 10)
par(mfrow = c(1,3))
plot(c(.50, .70, .90), bias_bart[7:9,1] , type = "b", ylim = range(bias_bart[7:9, c(1,2,4,5)]), ylab = "Average Relative Bias on Diagonals", xlab = expression(rho[x]), cex.lab = 1.2)
lines(c(.50, .70, .90), bias_bart[7:9,2] , type = "b", col = "blue", lty = 1)
lines(c(.50, .70, .90), bias_bart[7:9,4] , type = "b", col = "blue", lty = 2)
lines(c(.50, .70, .90), bias_bart[7:9,5] , type = "b", col = "red", lty = 1)
legend("bottomleft", bty = "n",col = c("black", "blue", "blue", "red"), lty = c(1,2,1,1), legend = c("BT", "Zero Lugsail", "Adapt Lugsail", "Over Lugsail"))

plot(c(.50, .70, .90),bias_tukey[7:9,1] , type = "b", ylim = range(bias_tukey[7:9, c(1,2,4,5)]), ylab = "Average Relative Bias on Diagonals", xlab = expression(rho[x]), pch = 2, cex.lab = 1.2)
lines(c(.50, .70, .90),bias_tukey[7:9,2] , type = "b", col = "blue", lty = 1, pch = 2)
lines(c(.50, .70, .90),bias_tukey[7:9,4] , type = "b", col = "blue", lty = 2, pch = 2)
lines(c(.50, .70, .90),bias_tukey[7:9,5] , type = "b", col = "red", lty = 1, pch = 2)
legend("bottomleft", bty = "n",col = c("black", "blue", "blue",  "red"), lty = c(1,2,1,1), legend = c("TH", "Zero Lugsail", "Adapt Lugsail", "Over Lugsail"))

plot(c(.50, .70, .90),bias_qs[7:9,1] , type = "b", ylim = range(bias_qs[7:9, c(1,2,4,5)]), ylab = "Average Relative Bias on Diagonals", xlab = expression(rho[x]), pch = 3, cex.lab = 1.2)
lines(c(.50, .70, .90),bias_qs[7:9,2] , type = "b", col = "blue", pch = 3, lty =1)
lines(c(.50, .70, .90),bias_qs[7:9,4] , type = "b", col = "blue", pch = 3, lty =2)
lines(c(.50, .70, .90),bias_qs[7:9,5] , type = "b", col = "red", pch = 3, lty =1)
legend("bottomleft", bty = "n",col = c("black", "blue", "blue", "red"), lty = c(1,2,1,1), legend = c("QS", "Zero Lugsail", "Adapt Lugsail", "Over Lugsail"))
dev.off()



### Table 3: Time-Series -- Time Comparisons


rm(list = ls())
load("Examples/TimeSeries/time_comp")
time <- matrix(0, nrow = 4, ncol = 9)
time[1, ] <- colMeans(time_bm)
time[2, ] <- colMeans(time_th) #bt was actually th. Unfortunate typographical error.
time[3, ] <- colMeans(time_bt)
time[4, ] <- colMeans(time_qs)
rownames(time) <- c("BM", "BT", "TH", "QS")
colnames(time) <- rep(c(5, 10, 50), 3)

round(time[ ,-c(1,4,7)], 3)






####################################
### VAR Example
####################################

load("Examples/VAR/var_coverage5e345p10_sq")


true.sig <- function(d = rep(1,p), p, omega = diag(1,p), phi = NULL)
{
	if(!is.matrix(phi))
	{
	dummy <- matrix(1:p^2, nrow = p, ncol = p)
	dummy <- qr.Q(qr(dummy))
	phi <- dummy %*% diag(d, nrow = p) %*% t(dummy)
	}
	
	variance <- matrix(solve(diag(1, p^2) - kronecker(phi, phi))%*%as.numeric(omega), nrow = p, ncol = p)
	final.cov <- solve(diag(1,p) - phi)%*%variance + variance%*%solve(diag(1,p) - t(phi)) - variance
	return(list(final.cov = final.cov, tar.var = variance, phi = phi, omega = omega))
}

p <- 10
phi <- diag(rep(.95,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))
true.var <- true.sig(p = p, omega = omega, phi = phi)$final.cov



### Figure 4: VAR  -- Bias
dum <- function(x)
{
	dum1 <- function(y)
	{
		foo <- matrix(y, ncol = 10, nrow = 10, byrow = TRUE)
		return(diag(foo))
	}
	rtn <- dum1(x)
	return(rtn)  #matrix(x[1,], ncol = 5, nrow = 5, byrow = TRUE)
}


b_bm <- numeric(length = 4)
b_wbm <- numeric(length = 4)
b_lug <- numeric(length = 4)
b_jack1 <- numeric(length = 4)
b_jack2 <- numeric(length = 4)

b_bm[1] <-  mean(rowMeans(apply(var.cov5e3$bias_bm, 1, dum))/diag(true.var)) 
b_bm[2] <- mean(rowMeans(apply(var.cov1e4$bias_bm, 1, dum))/diag(true.var)) 
b_bm[3] <- mean(rowMeans(apply(var.cov5e4$bias_bm, 1, dum))/diag(true.var)) 
b_bm[4] <- mean(rowMeans(apply(var.cov1e5$bias_bm, 1, dum))/diag(true.var)) 

b_wbm[1] <- mean(rowMeans(apply(var.cov5e3$bias_wbm, 1, dum))/diag(true.var)) 
b_wbm[2] <- mean(rowMeans(apply(var.cov1e4$bias_wbm, 1, dum))/diag(true.var)) 
b_wbm[3] <- mean(rowMeans(apply(var.cov5e4$bias_wbm, 1, dum))/diag(true.var)) 
b_wbm[4] <- mean(rowMeans(apply(var.cov1e5$bias_wbm, 1, dum))/diag(true.var)) 

b_lug[1] <- mean(rowMeans(apply(var.cov5e3$bias_lug, 1, dum))/diag(true.var)) 
b_lug[2] <- mean(rowMeans(apply(var.cov1e4$bias_lug, 1, dum))/diag(true.var)) 
b_lug[3] <- mean(rowMeans(apply(var.cov5e4$bias_lug, 1, dum))/diag(true.var)) 
b_lug[4] <- mean(rowMeans(apply(var.cov1e5$bias_lug, 1, dum))/diag(true.var)) 


b_jack1[1] <- mean(rowMeans(apply(var.cov5e3$bias_jacka, 1, dum))/diag(true.var)) 
b_jack1[2] <- mean(rowMeans(apply(var.cov1e4$bias_jacka, 1, dum))/diag(true.var)) 
b_jack1[3] <- mean(rowMeans(apply(var.cov5e4$bias_jacka, 1, dum))/diag(true.var)) 
b_jack1[4] <- mean(rowMeans(apply(var.cov1e5$bias_jacka, 1, dum))/diag(true.var)) 


b_jack2[1] <- mean(rowMeans(apply(var.cov5e3$"bias_jack-loga", 1, dum))/diag(true.var)) 
b_jack2[2] <- mean(rowMeans(apply(var.cov1e4$"bias_jack-loga", 1, dum))/diag(true.var)) 
b_jack2[3] <- mean(rowMeans(apply(var.cov5e4$"bias_jack-loga", 1, dum))/diag(true.var)) 
b_jack2[4] <- mean(rowMeans(apply(var.cov1e5$"bias_jack-loga", 1, dum))/diag(true.var)) 


foo <- cbind(b_bm, b_wbm, b_jack1, b_jack2, b_lug)
round(foo, 4)
sizes <- c(5e3, 1e4, 5e4, 1e5)
pdf("plots/var_bias.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(sizes, b_bm, type = "b", col = "black", ylim = c(-.3,.22), ylab = "Average bias on diagonals", xlab = "Chain length")
lines(sizes, b_wbm, col = "blue", type = "b", lty = 1)
# lines(sizes, b_jack1, col = "purple", type = "b", lty = 1)
lines(sizes, b_jack2, col = "purple", type = "b", lty = 1)
lines(sizes, b_lug, col = "red", type = "b", lty = 1)
abline(h = 0, lty = 2)
legend("topright", bty = "n",legend = c("BM", "Zero Lugsail", "Adapt Lugsail", "Over Lugsail"), col = c("black", "blue", "purple", "red"), lty = 1)
dev.off()



### Figure 4: VAR  -- ESS

rm(list = ls())
load("Examples/VAR/var_ess6_sq")
se_ess_bm <- apply(ess_track$BM, 2, sd)/sqrt(100)
se_ess_wbm <- apply(ess_track$WBM, 2, sd)/sqrt(100)
se_ess_lug <- apply(ess_track$LUG, 2, sd)/sqrt(100)
se_ess_jack1 <- apply(ess_track$Jack_lin, 2, sd)/sqrt(100)
se_ess_jack2 <- apply(ess_track$Jack_log, 2, sd)/sqrt(100)

pdf("plots/var_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track$BM), ylim = c(.0242, .031), 
	type = 'l', xlab = "Time", ylab = expression(hat(ESS)/n))
segments(x0 = subsize, y0 = colMeans(ess_track$BM) - 1.96*se_ess_bm, y1 = colMeans(ess_track$BM) + 1.96*se_ess_bm)

lines(subsize,colMeans(ess_track$WBM), col = "blue" )
segments(x0 = subsize, y0 = colMeans(ess_track$WBM) - 1.96*se_ess_wbm, y1 = colMeans(ess_track$WBM) + 1.96*se_ess_wbm, col = "blue")

lines(subsize,colMeans(ess_track$LUG), col = "red" )
segments(x0 = subsize, y0 = colMeans(ess_track$LUG) - 1.96*se_ess_lug, y1 = colMeans(ess_track$LUG) + 1.96*se_ess_lug, col = "red")

lines(subsize,colMeans(ess_track$Jack_log), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track$Jack_log) - 1.96*se_ess_jack2, y1 = colMeans(ess_track$Jack_log) + 1.96*se_ess_jack2, col = "purple")

abline(h = ess_true, lty = 2)
legend("topright", bty = "n",legend = c("BM", "Zero Lugsail", "Adapt Lugsail", "Over Lugsail"), col = c("black", "blue", "purple", "red"), lty = 1)
dev.off()







##############################################
####### Framingham Logistic regression ########
#############################################

### Figure 5: Framingham  -- ESS

rm(list = ls())
load("Examples/Framlogistic/fram_ess6_sq")
se_ess_bm <- apply(ess_track$BM, 2, sd)/sqrt(100)
se_ess_wbm <- apply(ess_track$WBM, 2, sd)/sqrt(100)
se_ess_lug <- apply(ess_track$LUG, 2, sd)/sqrt(100)
se_ess_jack1 <- apply(ess_track$Jack_lin, 2, sd)/sqrt(100)
se_ess_jack2 <- apply(ess_track$Jack_log, 2, sd)/sqrt(100)

pdf("plots/fram_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(subsize, colMeans(ess_track$BM), ylim = c(.017, .021), 
	type = 'l', xlab = "Time", ylab = expression(hat(ESS)/n))
segments(x0 = subsize, y0 = colMeans(ess_track$BM) - 1.96*se_ess_bm, y1 = colMeans(ess_track$BM) + 1.96*se_ess_bm)

lines(subsize,colMeans(ess_track$WBM), col = "blue" )
segments(x0 = subsize, y0 = colMeans(ess_track$WBM) - 1.96*se_ess_wbm, y1 = colMeans(ess_track$WBM) + 1.96*se_ess_wbm, col = "blue")

lines(subsize,colMeans(ess_track$LUG), col = "red" )
segments(x0 = subsize, y0 = colMeans(ess_track$LUG) - 1.96*se_ess_lug, y1 = colMeans(ess_track$LUG) + 1.96*se_ess_lug, col = "red")

lines(subsize,colMeans(ess_track$Jack_log), col = "purple" )
segments(x0 = subsize, y0 = colMeans(ess_track$Jack_log) - 1.96*se_ess_jack2, y1 = colMeans(ess_track$Jack_log) + 1.96*se_ess_jack2, col = "purple")


legend("topright", bty = "n",legend = c("BM", "Zero Lugsail", "Adapt Lugsail", "Over Lugsail"), col = c("black", "blue", "purple", "red"), lty = 1)
dev.off()


################################
## Supplementary material ######

## Bias and variance for SV estimators###

# Bartlett window #

biasvar <- function(method = "original", a = NULL)
{
	bias <- numeric(length = 3)
	vars <- numeric(length = 3)

	if(method == "original")
	{
		r <- 1
		c <- 1/2

		bias[1] <- (1- r*c)/(1-c) * 1
		bias[2] <- (1- r^2*c)/(1-c) * pi^2/4
		bias[3] <- (1- r^2*c)/(1-c) * 1.4212

		# Variance for bartlett , TH, QS
		vars[1] <- 2 * (1 + c^2/r - 3*c/r + c/r^2)/(3 * (1 - c)^2)
		vars[2] <- 3/4 # 3/(4 * (1 - c)^2) * (1 + c^2/r - 4*c/3 * (r^3 *sinpi(1/r) + pi*r^2 - pi)/(pi*r^3 - pi*r))
		vars[3] <- (1-c)^(-2) * (1 + c^2/r - c*((r-1)^3*(1 + 3*r + r^2) - (r+1)^3*(1 - 3*r+r^2))/(4*r^3))

	}

	if(method == "zero")
	{
		r <- 2
		c <- 1/2

		bias[1] <- (1- r*c)/(1-c) * 1
		vars[1] <- 2 * (1 + c^2/r - 3*c/r + c/r^2)/(3 * (1 - c)^2)		

		r <- 2
		c <- 1/4	
		bias[2] <- (1- r^2*c)/(1-c) * pi^2/4
		bias[3] <- (1- r^2*c)/(1-c) * 1.4212

		# Variance for bartlett , TH, QS
		vars[2] <-  3/(4 * (1 - c)^2) * (1 + c^2/r - 4*c/3 * (r^3 *sinpi(1/r) + pi*r^2 - pi)/(pi*r^3 - pi*r))
		vars[3] <- (1-c)^(-2) * (1 + c^2/r - c*((r-1)^3*(1 + 3*r + r^2) - (r+1)^3*(1 - 3*r+r^2))/(4*r^3))

	}

	if(method == "over")
	{
		r <- 3
		c <- 1/2

		bias[1] <- (1- r*c)/(1-c) * 1
		vars[1] <- 2 * (1 + c^2/r - 3*c/r + c/r^2)/(3 * (1 - c)^2)		

		r <- 3
		c <- 1/5		
		bias[2] <- (1- r^2*c)/(1-c) * pi^2/4
		bias[3] <- (1- r^2*c)/(1-c) * 1.4212

		# Variance for bartlett , TH, QS
		vars[2] <- 3/(4 * (1 - c)^2) * (1 + c^2/r - 4*c/3 * (r^3 *sinpi(1/r) + pi*r^2 - pi)/(pi*r^3 - pi*r))
		vars[3] <- (1-c)^(-2) * (1 + c^2/r - c*((r-1)^3*(1 + 3*r + r^2) - (r+1)^3*(1 - 3*r+r^2))/(4*r^3))

	}

	if(method == "adapt")
	{
		r <- 2
		c <- (log(a) + 1)/(r*log(a) + 1)

		bias[1] <- (1- r*c)/(1-c) * 1
		vars[1] <- 2 * (1 + c^2/r - 3*c/r + c/r^2)/(3 * (1 - c)^2)		

		c <- (log(a) + 1)/(r^2*log(a) + 1)	
		bias[2] <- (1- r^2*c)/(1-c) * pi^2/4
		bias[3] <- (1- r^2*c)/(1-c) * 1.4212

		# Variance for bartlett , TH, QS
		vars[2] <-  3/(4 * (1 - c)^2) * (1 + c^2/r - 4*c/3 * (r^3 *sinpi(1/r) + pi*r^2 - pi)/(pi*r^3 - pi*r))
		vars[3] <- (1-c)^(-2) * (1 + c^2/r - c*((r-1)^3*(1 + 3*r + r^2) - (r+1)^3*(1 - 3*r+r^2))/(4*r^3))

	}

	names(bias) <- c("BT", "TH", "QS")
	names(vars) <- c("BT", "TH", "QS")

	return(vars)
}

# Original 
out <- cbind(biasvar(), biasvar("zero"), biasvar("adapt", 100), biasvar("adapt", 1e5), biasvar("over", 10))
colnames(out) <- c("OG", "Zero", "ad1e2", "ad1e4", "over")
round(out, 3)



