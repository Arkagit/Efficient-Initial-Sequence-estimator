#install.packages("HadamardR")
library(HadamardR)

########## Initialization
rho_max = 1.01
p = 12
omega = diag(rep(1, p))


########## Product matrix of VAR process
phi = (1/p)*Hadamard_Matrix(p)%*%diag(rho_max^(-(1:p)))%*%t(Hadamard_Matrix(p))


# True stationary covariance matrix
true.sig.gen <- function(p, omega, phi){
  omega = as.matrix(omega)
  phi = as.matrix(phi)
  variance <- matrix(solve(diag(1, p^2) - kronecker(phi, phi))%*%as.numeric(omega), nrow = p, ncol = p)
  final.cov <- solve(diag(1,p) - phi)%*%variance + variance%*%solve(diag(1,p) - t(phi)) - variance
  return(list(final.cov = final.cov, tar.var = variance, phi = phi, omega = omega))
}

#true.sig.gen(p, omega, phi)