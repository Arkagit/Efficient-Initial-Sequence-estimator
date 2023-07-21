source("Gen_data.R")

#Variables
num = 1e5
rho = 0.95
dim = 8
reps <- 100

# Initializing ESS
ess_bm = numeric(0)
ess_rbm = numeric(0)
ess_ise = numeric(0)
ess_cv = numeric(0)

# Initializing time
time_bm = numeric(0)
time_rbm = numeric(0)
time_ise = numeric(0)
time_cv = numeric(0)


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
  time_ise[j] = Sys.time() - start
  ess_ise[j] = multiESS(data, covmat = IS)
  
  print(paste(3))
  
  # New variance ESS and time
  start <- Sys.time()
  BM = mcse.multi(data, method = "bm", r = 1)$cov
  corrMat <- cov2cor(BM)
    
  sds = apply(data, 2, function(l) sqrt(initseq(l)$var.pos))
  var_est = diag(sds)%*%corrMat%*%diag(sds)
  time_cv[j] = Sys.time() - start
  
  if(det(var_est) > 0){
    ess_cv[j] = multiESS(data, covmat = var_est)
  }else{
    ess_cv[j] = 0
  }
  
  print(paste(4))
}

# Average ESS
mean(ess_bm); mean(ess_ise); mean(ess_cv)
foo <- cbind(ess_bm, ess_rbm,  ess_ise, ess_cv)
boxplot(foo, main = "Boxplot of estimated ESS")
abline(h = n*(det(MAT(8,0.95)[[1]])/det(MAT(8,0.95)[[2]]))^(1/dim), col = "red")

# Average Time
mean(time_bm); mean(time_ise); mean(time_cv)

pdf("/Users/arkabanerjee/Desktop/myplot.pdf", width=20, height=10,onefile=T)
par(mfrow = c(1,2))
dev.off()
