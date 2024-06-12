set.seed(100)
source("Coverage_func.R")

library(foreach)
library(doParallel)



B <- 1000
#p <- 10
#rho <- 0.99
#phi <- diag(rep(rho,p))
#omega <- matrix(.9, nrow = p, ncol = p)
#diag(omega) <- 1
#omega <- omega^(abs(col(omega)-row(omega)))
subsize <- c(5e3, 8e3, 1e4, 3e4, 5e4, 8e4, 1e5, 3e5, 5e5)
level <- 0.95
Coverage = matrix(0, nrow = length(subsize), ncol = 6)
Timer = matrix(0, nrow = length(subsize), ncol = 6)
#length(colMeans(minichain))

parallel::detectCores()
n.cores <- 50
doParallel::registerDoParallel(cores = n.cores)


cover = foreach(j = 1:B, .packages = c("mcmcse"))%dopar%{
	print(j)
	Out = list()
	time = list()
	cover = coverage(subsize, phi, omega, level)
	cover
}
#chart <- Out/B
#chart
#time
#cover

#for(i in 1:B){
#	Coverage = Coverage + cover[[i]][[1]]
#	Timer = Timer + cover[[i]][[2]]
#}

save(cover, subsize, B, file = "dat_matices.Rdata")

#save("coverage table" = chart, file = "Coverage.Rdata")