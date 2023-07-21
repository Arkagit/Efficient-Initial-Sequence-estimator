library(MASS)
library(mcmcse)
library(ellipse)
library(mcmc)
library(HadamardR)

p=8
# Variance for Error
omega = diag(1, p)

#Asymptotic Variances
MAT = function(dimension, rho){
  #Product Vector
  phi = (1/dimension)*Hadamard_Matrix(dimension)%*%diag(rho^(1:dimension))%*%t(Hadamard_Matrix(dimension))
  V = solve(diag(1,dimension) - phi%*%phi)%*%omega
  sig = 2*solve(diag(1,dimension) - phi)%*%solve(diag(1,dimension) - phi%*%phi)%*%omega - solve(diag(1,dimension) - phi%*%phi)%*%omega
  return(list(V, sig))
}

#Data Generation
data_g = function(number, dimension, rho){
  #Product Vector
  phi = (1/dimension)*Hadamard_Matrix(dimension)%*%diag(rho^(1:dimension))%*%t(Hadamard_Matrix(dimension))
  #Error Vector
  err = mvrnorm(number, mu = rep(0,dimension), Sigma = omega)
  data = matrix(0, number, dimension) 
  data[1,] = rep(1,dimension)
  for (i in 2:number){ 
    data[i,] = t(phi)%*%data[i-1,] + err[i,]
  }
  return(data)
}


