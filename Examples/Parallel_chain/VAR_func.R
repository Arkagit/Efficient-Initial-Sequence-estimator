########### Calling libraries
library(mcmcse)
library(Rcpp)
library(MASS)

########### Sourcing files
# source("Asymp_var.R")


# Generated Data from Vector Auto-covariance model
var1 <- function(p = 2, phi, nsim = 1e3,  omega = diag(rep(1, p)), start = rep(0,p))
{
  # sig <- diag(rep(rho, p))
  chain <- matrix(0, nrow = nsim, ncol = p)	
  decomp <- svd(omega)
  omega.root <- decomp$v %*% diag((decomp$d)^(1/2), nrow = p, ncol = p) %*% t(decomp$u)
  
  variance <- matrix(solve(diag(1, p^2) - kronecker(phi, phi))%*%as.numeric(omega), nrow = p, ncol = p)
  # chain[1, ] <- start
  chain[1, ] <- mvrnorm(1, mu = rep(0, p), Sigma = variance)
  # sigma.root%*%rnorm(p)
  for(i in 2:nsim)
  {
    chain[i, ] <- phi%*%chain[i-1, ] + omega.root %*% rnorm(p)
  }
  return(chain)
}


# True stationary covariance matrix
#true.sig <- function(d = rep(1,p), p, omega = diag(1,p), phi = NULL)
#{
#  if(!is.matrix(phi))
#  {
#    dummy <- matrix(1:p^2, nrow = p, ncol = p)
#    dummy <- qr.Q(qr(dummy))
#    phi <- dummy %*% diag(d, nrow = p) %*% t(dummy)
#  }
#  
#  variance <- matrix(solve(diag(1, p^2) - kronecker(phi, phi))%*%as.numeric(omega), nrow = p, ncol = p)
#  final.cov <- solve(diag(1,p) - phi)%*%variance + variance%*%solve(diag(1,p) - t(phi)) - variance
#  return(list(final.cov = final.cov, tar.var = variance, phi = phi, omega = omega))
#}


