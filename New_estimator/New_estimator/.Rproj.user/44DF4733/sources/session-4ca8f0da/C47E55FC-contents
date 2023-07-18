sourceCpp("Mat_prod.cpp")

library(mcmcse)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(microbenchmark)

# Function for computing univariate ISE from acf vector
cppFunction('double Ise (NumericVector x){
  int i = 1;
  double y = x[0];
  while (x[i] > 0) {
  y = y + 2*x[i];
  i++;
}
  return y;
}')


#new covariance matrix
cov.sig <- function(data){
  d = numeric(0)
  n = dim(data)[1]; p = dim(data)[2]
  for (i in 1:p) {
    acf_data = acf(data[,i], type = "covariance", plot = "FALSE")$acf
    d[i] = Ise(acf_data)
  }
  covariance = eigenMapMatMult(eigenMapMatMult(diag(d), cov2cor(mcse.multi(data, method = "bm")$cov)), diag(d))
  return(list("covariance" = covariance, "est" = colMeans(data)))
}

