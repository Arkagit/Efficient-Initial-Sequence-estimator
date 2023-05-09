library(MASS)
library(mcmcse)
library(ellipse)
library(mcmc)

# datasize
n = 5e4
# dimensions
p = 5

# Product vector
#phi = diag(c(0.9, 0.8, 0.1, 0.1, 0.1))

# Variance for Error
rho = 0.95
omega = matrix(1, ncol = p, nrow = p, byrow = TRUE)
for (k in 1:p) {
  for (l in 1:p) {
    omega[k,l] = rho^(abs(k-l))
  }
}
phi <- omega
phi <- phi/max(eigen(phi)$values + .1)
vec_omega = numeric(0)

for (k in 1:p) {
  vec_omega = append(vec_omega, omega[k,])
}
vec_V = solve(diag(rep(1, p^(2))) - kronecker(phi,phi))%*%vec_omega
V = matrix(vec_V, nrow = p, byrow = TRUE)

# Asymptotic variance
sig = solve(diag(rep(1,p)) - phi)%*%V + V%*%solve(diag(rep(1,p)) - phi) - V

# VAR process
err = mvrnorm(n, mu = c(0,0,0,0,0), Sigma = omega)
data = matrix(0, n, p) 
data[1,] = c(1,1,1,1,1)
for (i in 2:n){ 
  data[i,] = t(phi)%*%data[i-1,] + err[i,]
}
head(data)

