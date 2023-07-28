source("Coverage_func.R")


B <- 1e3
p <- 10
rho <- 0.95
phi <- diag(rep(rho,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))
N=1e5
#length(colMeans(minichain))

chart <- coverage(N, phi, omega, B, 0.90)
save("coverage table" = chart, file = "Coverage.Rdata")