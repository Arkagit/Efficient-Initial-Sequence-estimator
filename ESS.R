source("Gen_data.R")

#Variables
num = 1e5
rho = 0.95
dim = 8
reps <- 100

# Initializing ESS
ess_bm = numeric(0)
ess_rbm = numeric(0)
ess_is = numeric(0)
ess_new = numeric(0)

# Initializing time
time_bm = numeric(0)
time_rbm = numeric(0)
time_is = numeric(0)
time_new = numeric(0)


for (j in 1:reps) {
  print(j)
  # VAR process
  data = data_g(num,dim,rho)
  
  print(paste(1))
  
  # Batch means ESS and time
  start <- Sys.time()
  BM = mcse.multi(data, method = "bm", r = 3)$cov
  time_rbm[j] = Sys.time() - start
  ess_rbm[j] = multiESS(data, method = "bm", r = 3)
  
  BM = mcse.multi(data, method = "bm", r = 1)$cov
  time_bm[j] = Sys.time() - start
  ess_bm[j] = multiESS(data, method = "bm", r = 1)
  
  print(paste(2))
  
  # ISE ESS and time
  start = Sys.time()
  IS = mcse.initseq(data)$cov
  time_is[j] = Sys.time() - start
  ess_is[j] = multiESS(data, covmat = IS)
  
  print(paste(3))
  
  # New variance ESS and time
  start <- Sys.time()
  BM = mcse.multi(data, method = "bm", r = 1)$cov
  corrMat <- cov2cor(BM)
    
  sds = apply(data, 2, function(l) sqrt(initseq(l)$var.pos))
  var_est = diag(sds)%*%corrMat%*%diag(sds)
  time_new[j] = Sys.time() - start
  
  if(det(var_est) > 0){
    ess_new[j] = multiESS(data, covmat = var_est)
  }else{
    ess_new[j] = 0
  }
  
  print(paste(4))
}

# Average ESS
mean(ess_bm); mean(ess_is); mean(ess_new)
foo <- cbind(ess_bm, ess_rbm,  ess_is, ess_new)
boxplot(foo, main = "Boxplot of estimated ESS")
abline(h = n*(det(MAT(8,0.95)[[1]])/det(MAT(8,0.95)[[2]]))^(1/p), col = "red")

# Average Time
mean(time_bm); mean(time_is); mean(time_new)

