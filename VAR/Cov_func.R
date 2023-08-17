# install.packages("devtools")
#devtools::install_github("hsong1/momentLS")
library(momentLS)
library(mcmcse)
library(mcmc)
library(Rcpp)
library(phonTools)
library(RcppArmadillo)
library(RcppEigen)



cppFunction('Rcpp::List Ise (Rcpp::NumericVector x){
  int i = 0;
  double y = - x[0];
  double z = x[2*i]+ x[2*i+1];
  while (z > 0) {
  y = y + 2*z;
  i = i +1;
  z = x[2*i]+ x[2*i+1];
}
  Rcpp::List out = Rcpp::List::create(y, i);
  return out;
}')


#new covariance matrix
cov.sig <- function(data, type = "geyer"){
  n = dim(data)[1]
  p = dim(data)[2]
  d = rep(0,p)
  lag = rep(0,p)

  if(type == "MomentLS"){
    for(i in 1:p){
      # compute the empirical autocovariances
      r = autocov(data[,i])
      # tune delta
      delta_tilde = tune_delta(data[,i], nSplits = 5, c_M_const = 0)$delta*0.8
      # fit MomentLS
      m = SR1(r, delta = delta_tilde) # fit
      d[i] = asympVariance(weights = m$weights, support = m$support)
    }
    covariance = sqrt(d) * cov2cor(mcse.multi(data, method = "bm", r = 1)$cov) * rep(sqrt(d), each = p)
    return(list("covariance" = covariance, "est" = colMeans(data)))

  } else {
    # compute the initial sequences
    for (i in 1:p) {
      Ise_data = Ise(fastacf(data[,i], show = "FALSE", window = "rectangular", lag.max = floor(n/2))$acf)
      d[i] = Ise_data[[1]]*var(data[,i])*(length(data[,i]) - 1)/length(data[,i])
      lag[i] = Ise_data[[2]]
    }
    covariance = sqrt(d) * cov2cor(mcse.multi(data, method = "bm", r = 1)$cov) * rep(sqrt(d), each = p)
    return(list("covariance" = covariance, "est" = colMeans(data), "stopping_lag" = lag))
  }
}


# Comparison with Geyer (1992)
#data = matrix(rnorm(2e2), nrow = 1e2)
#initseq(data[,1])$var.pos; initseq(data[,2])$var.pos
#cov.sig(data,"geyer")
#cov.sig(data, "MomentLS")

