set.seed(100)
########## Sourcing fiiles
source("VAR_bias.R")

######### Calling libraries
library(mcmcse)
library(foreach)
library(doParallel)

### True Covariance Matrix
true.var <- true.sig.gen(p = p, omega = omega, phi = phi)$final.cov

# Number of replications
B <- 1000

dum <- function(x){
  dum1 <- function(y)
  {
    foo <- matrix(y, ncol = p, nrow = p, byrow = TRUE)
    return(diag(foo))
  }
  rtn <- dum1(x)
  return(rtn)  #matrix(x[1,], ncol = 5, nrow = 5, byrow = TRUE)
}

# Biases at different chain-length
subsize = c(1e4, 5e4, 1e5, 5e5, 1e6)
# Initializing output list
BS = list()

# Setting up parallelizing config
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n.cores)

# Parallel loop
BS = foreach(b = 1:length(subsize), .packages = "mcmcse")%dopar%{
	var.cov = var_coverage(N = subsize[b], phi = phi, omega = omega, B = B)
	print(subsize[b])
	b_bm <- mean(rowMeans(apply(var.cov$bias_bm, 1, dum))/diag(true.var))
	#b_lug <- mean(rowMeans(apply(var.cov$bias_lug, 1, dum))/diag(true.var)) 
	b_ise <- mean(rowMeans(apply(var.cov$bias_ise, 1, dum))/diag(true.var))
	b_cc <- mean(rowMeans(apply(var.cov$bias_cc, 1, dum))/diag(true.var))
	b_sve <- mean(rowMeans(apply(var.cov$bias_sve, 1, dum))/diag(true.var)) 
	b_mls <- mean(rowMeans(apply(var.cov$bias_mls, 1, dum))/diag(true.var))
	
	list(b_bm, b_ise, b_cc, b_sve, b_mls)
}


# Saving bias data in Rdata file
save(BS, subsize, file = "bias_data.Rdata")
