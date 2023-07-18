source("VAR_bias.R")


# Parameter values
p <- 10
phi <- diag(rep(.95,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))
true.var <- true.sig(p = p, omega = omega, phi = phi)$final.cov

# Number of replications
B <- 1000

# Biases at different chain-length
var.cov5e3 <- var_coverage(N = 5e3, phi = phi, omega = omega, B = B)
var.cov1e4 <- var_coverage(N = 1e4, phi = phi, omega = omega, B = B)
var.cov5e4 <- var_coverage(N = 5e4, phi = phi, omega = omega, B = B)
var.cov1e5 <- var_coverage(N = 1e5, phi = phi, omega = omega, B = B)


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

# Memory allocation for biases for different estimation methods
b_bm <- numeric(length = 4)
b_wbm <- numeric(length = 4)
b_lug <- numeric(length = 4)
b_cc <- numeric(length = 4)


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


b_cc[1] <- mean(rowMeans(apply(var.cov5e3$bias_cc, 1, dum))/diag(true.var)) 
b_cc[2] <- mean(rowMeans(apply(var.cov1e4$bias_cc, 1, dum))/diag(true.var)) 
b_cc[3] <- mean(rowMeans(apply(var.cov5e4$bias_cc, 1, dum))/diag(true.var)) 
b_cc[4] <- mean(rowMeans(apply(var.cov1e5$bias_cc, 1, dum))/diag(true.var)) 




foo <- cbind(b_bm, b_wbm, b_cc, b_lug)
round(foo, 4)
sizes <- c(5e3, 1e4, 5e4, 1e5)
#pdf("plots/var_bias.pdf", height = 6, width = 6)
par(mfrow = c(1,1))
plot(sizes, b_cc, type = "b", col = "black", ylab = "Average bias on diagonals", xlab = "Chain length")
lines(sizes, b_wbm, col = "blue", type = "b", lty = 1)
# lines(sizes, b_jack1, col = "purple", type = "b", lty = 1)
lines(sizes, b_cc, col = "purple", type = "b", lty = 1)
lines(sizes, b_lug, col = "red", type = "b", lty = 1)
abline(h = 0, lty = 2)
legend("topright", bty = "n",legend = c("BM", "Zero Lugsail", "Cov-Corr", "Over Lugsail"), col = c("black", "blue", "purple", "red"), lty = 1)
#dev.off()

b_cc
b_bm
b_wbm
b_lug
