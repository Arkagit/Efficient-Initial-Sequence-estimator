library(MASS)
library(mcmcse)
library(ellipse)
library(mcmc)

# datasize
n = 1e7
# dimensions
p = 5
# Product vector
phi = diag(c(0.9, 0.5, 0.1, 0.1, 0.1))
# Variance for Error
rho = 0.9
vec_omega = c(1, rho, rho^(2), rho^(3), rho^(4),
              rho, 1, rho, rho^(2), rho^(3),
              rho^(2), rho, 1, rho, rho^(2),
              rho^(3), rho^(2), rho, 1, rho,
              rho^(4), rho^(3), rho^(2), rho, 1)

omega = matrix(vec_omega, nrow = 5, byrow = TRUE)
eigen(omega)
vec_V = solve(diag(rep(1, p^(2))) - kronecker(phi,phi))%*%vec_omega
V = matrix(vec_V, nrow = 5, byrow = TRUE)
# Asymptotic variance
sig = solve(diag(rep(1,p)) - phi)%*%V + V%*%solve(diag(rep(1,p)) - phi) - V

# VAR process
err = mvrnorm(n, mu = c(0,0,0,0,0), Sigma = omega)

data = matrix(0, n, p) 
data[1,] = c(1,1,1,1,1)
for (i in 2:n){ 
  data[i,] = phi%*%data[i-1,] + err[i,]
}
head(data)

