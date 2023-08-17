set.seed(100)
source("Coverage_func.R")


B <- 1e3
p <- 10
rho <- 0.99
phi <- diag(rep(rho,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))
subsize <- c(5e3, 1e4, 5e4)
level <- 0.95
Out = matrix(0, nrow = 3, ncol = 6)
#length(colMeans(minichain))

for(j in 1:B){
	print(j)
	Out = Out + coverage(subsize, phi, omega, level)
}
chart <- Out/B
chart
#save("coverage table" = chart, file = "Coverage.Rdata")