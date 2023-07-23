source("VAR_bias.R")


# Parameter values
p <- 10
phi <- diag(rep(.95,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))
true.var <- true.sig(p = p, omega = omega, phi = phi)$final.cov

# Number of replications
B <- 100

# Biases at different chain-length
var.cov1e3 <- var_coverage(N = 1e3, phi = phi, omega = omega, B = B)
var.cov5e3 <- var_coverage(N = 5e3, phi = phi, omega = omega, B = B)
var.cov1e4 <- var_coverage(N = 1e4, phi = phi, omega = omega, B = B)
var.cov5e4 <- var_coverage(N = 5e4, phi = phi, omega = omega, B = B)


dum <- function(x){
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
b_ise <- numeric(length = 4)


# Storing bias values for different methods at different chain lengths
b_bm[1] <-  mean(rowMeans(apply(var.cov1e3$bias_bm, 1, dum))/diag(true.var)) 
b_bm[2] <- mean(rowMeans(apply(var.cov5e3$bias_bm, 1, dum))/diag(true.var)) 
b_bm[3] <- mean(rowMeans(apply(var.cov1e4$bias_bm, 1, dum))/diag(true.var)) 
b_bm[4] <- mean(rowMeans(apply(var.cov5e4$bias_bm, 1, dum))/diag(true.var)) 

b_wbm[1] <- mean(rowMeans(apply(var.cov1e3$bias_wbm, 1, dum))/diag(true.var)) 
b_wbm[2] <- mean(rowMeans(apply(var.cov5e3$bias_wbm, 1, dum))/diag(true.var)) 
b_wbm[3] <- mean(rowMeans(apply(var.cov1e4$bias_wbm, 1, dum))/diag(true.var)) 
b_wbm[4] <- mean(rowMeans(apply(var.cov5e4$bias_wbm, 1, dum))/diag(true.var)) 

b_lug[1] <- mean(rowMeans(apply(var.cov1e3$bias_lug, 1, dum))/diag(true.var)) 
b_lug[2] <- mean(rowMeans(apply(var.cov5e3$bias_lug, 1, dum))/diag(true.var)) 
b_lug[3] <- mean(rowMeans(apply(var.cov1e4$bias_lug, 1, dum))/diag(true.var)) 
b_lug[4] <- mean(rowMeans(apply(var.cov5e4$bias_lug, 1, dum))/diag(true.var)) 


b_cc[1] <- mean(rowMeans(apply(var.cov1e3$bias_cc, 1, dum))/diag(true.var)) 
b_cc[2] <- mean(rowMeans(apply(var.cov5e3$bias_cc, 1, dum))/diag(true.var)) 
b_cc[3] <- mean(rowMeans(apply(var.cov1e4$bias_cc, 1, dum))/diag(true.var)) 
b_cc[4] <- mean(rowMeans(apply(var.cov5e4$bias_cc, 1, dum))/diag(true.var)) 


b_ise[1] <-  mean(rowMeans(apply(var.cov1e3$bias_ise, 1, dum))/diag(true.var)) 
b_ise[2] <- mean(rowMeans(apply(var.cov5e3$bias_ise, 1, dum))/diag(true.var)) 
b_ise[3] <- mean(rowMeans(apply(var.cov1e4$bias_ise, 1, dum))/diag(true.var)) 
b_ise[4] <- mean(rowMeans(apply(var.cov5e4$bias_ise, 1, dum))/diag(true.var)) 



# Saving bias data in Rdata file
save(b_bm, b_wbm, b_lug, b_cc, b_ise, file = "bias_data.Rdata")
