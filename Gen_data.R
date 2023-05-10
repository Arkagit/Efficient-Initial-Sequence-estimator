library(MASS)
library(mcmcse)
library(ellipse)
library(mcmc)
library(HadamardR)

# datasize
n = 1e7
# dimensions
p = 8

# Product vector
rho = 0.95
phi = (1/p)*Hadamard_Matrix(p)%*%diag(rho^(1:p))%*%t(Hadamard_Matrix(p))

# Variance for Error
omega = diag(1, p)
#omega = matrix(1, ncol = p, nrow = p, byrow = TRUE)
#for (k in 1:p) {
#  for (l in 1:p) {
#    omega[k,l] = rho^(abs(k-l))
#  }
#}
#phi <- omega
#phi <- phi/max(eigen(phi)$values + .1)
#vec_omega = numeric(0)

#for (k in 1:p) {
#  vec_omega = append(vec_omega, omega[k,])
#}
#vec_V = solve(diag(rep(1, p^(2))) - kronecker(phi,phi))%*%vec_omega
#V = matrix(vec_V, nrow = p, byrow = TRUE)

# Asymptotic variance
#sig = (diag(1,p) + 2*phi%*%phi%*%solve(diag(1,p) - phi%*%phi))%*%solve(diag(1,p) - phi%*%phi)%*%omega
sig = 2*solve(diag(1,p) - phi%*%phi)%*%solve(diag(1,p) - phi%*%phi)%*%omega - solve(diag(1,p) - phi%*%phi)%*%omega

# VAR process
err = mvrnorm(n, mu = rep(0,p), Sigma = omega)
data = matrix(0, n, p) 
data[1,] = rep(1,p)
for (i in 2:n){ 
  data[i,] = t(phi)%*%data[i-1,] + err[i,]
}
head(data)

