# Initializing ESS
ess_bm = numeric(0)
ess_is = numeric(0)
ess_new = numeric(0)

# Initializing time
time_bm = numeric(0)
time_is = numeric(0)
time_new = numeric(0)


for (j in 1:100) {
  source("Gen_data.R")
  n = 1e6
  data = data[1:n,]
  
  print(paste(1))
  
  # Batch means ESS and time
  start <- Sys.time()
  BM = mcse.multi(data, method = "bm", r = 3, size = "sqroot")$cov
  time_bm[j] = Sys.time() - start
  ess_bm[j] = multiESS(data, method = "bm", r = 3, size = "sqroot")
  
  print(paste(2))
  
  # ISE ESS and time
  start = Sys.time()
  IS = mcse.initseq(data)$cov
  time_is[j] = Sys.time() - start
  ess_is[j] = multiESS(data, covmat = IS)
  
  print(paste(3))
  
  # New variance ESS and time
  start <- Sys.time()
  BM = cor(mcse.multi(data, method = "bm", r = 1, size = "sqroot")$cov)
  sd = numeric(0)
  for (i in 1:p) {
    sd[i] = sqrt(initseq(data[,i])$var.pos)
  }
  var_est = diag(sd)%*%corrm%*%diag(sd)
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

# Average Time
mean(time_bm); mean(time_is); mean(ess_new)
