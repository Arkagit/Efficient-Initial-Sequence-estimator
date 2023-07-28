library(mcmcse)
library(mcmc)
library(Rcpp)
library(phonTools)
library(RcppArmadillo)
library(RcppEigen)

#library(microbenchmark)

#sourceCpp("Mat_prod.cpp")

# Function for computing univariate ISE from acf vector

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
cov.sig <- function(data){
  n = dim(data)[1]
  p = dim(data)[2]
  d = rep(0,p)
  lag = rep(0,p)
  for (i in 1:p) {
    Ise_data = Ise(fastacf(data[,i], show = "FALSE", window = "rectangular", lag.max = floor(n/2))$acf)
    d[i] = Ise_data[[1]]*var(data[,i])*(length(data[,i]) - 1)/length(data[,i])
    lag[i] = Ise_data[[2]]
  }
  covariance = sqrt(d) * cov2cor(mcse.multi(data, method = "bm", r = 1)$cov) * rep(sqrt(d), each = p)
  return(list("covariance" = covariance, "est" = colMeans(data), "stopping_lag" = lag))
}


# Comparison with Geyer (1992)
#data = matrix(rnorm(2e2), nrow = 1e2)
#initseq(data[,1])$var.pos; initseq(data[,2])$var.pos
#cov.sig(data)

