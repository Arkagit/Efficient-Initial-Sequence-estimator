source("Gen_data.R")

# datasize
n = 1e5
dim = 8
rho = 0.95

data = data_g(n, dim, rho)

# Batch means estimation
BME = mcse.multi(data, method = "bm", r = 3)
BME$cov

# ISE estimation 
IS = mcse.initseq(data)
IS$cov

# correlation matrix by batch means
BM = cov2cor(mcse.multi(data, method = "bm", r = 1)$cov)
# Standard deviation bt ISE
sd = numeric(0)
for (i in 1:p) {
  sd[i] = sqrt(initseq(data[,i])$var.pos)
}
# new variance estimate
var_est = diag(sd)%*%BM%*%diag(sd)


# Plots for 1st and 3rd component
new_v13 = var_est[c(1,dim), c(1,dim)]
pl3 = ellipse(new_v13, centre = c(0,0), which = c(1, 2), npoints = 100)
plot(pl3, type = "l", col = "blue", main = "Confidence Ellipsoid",
     xlab = "1st Component", ylab = paste(dim,"th Component"))

ise13 = IS$cov[c(1,dim), c(1,dim)]
pl2 = ellipse(ise13, centre = c(0,0), which = c(1, 2), npoints = 100)
lines(pl2, type = "l", col = "green", main = "Confidence Ellipsoid")

b13 = BME$cov[c(1,dim), c(1,dim)]
pl1 = ellipse(b13, centre = c(0,0), which = c(1, 2), npoints = 100)
lines(pl1, type = "l", col = "red", main = "Confidence Ellipsoid")

tr13 = (MAT(dim, rho)[[2]])[c(1,dim),c(1,dim)] 
pl0 = ellipse(tr13, centre = c(0,0), which = c(1, 2), npoints = 100)
lines(pl0, type = "l", col = "black", main = "Confidence Ellipsoid")

legend(x = "topleft", lty = c(1,1,1,1), text.font = 0.5, 
       col= c("blue","green", "red", "black"),text.col = "black", 
       legend=c("New", "ISE", "Lugsail BM", "True"))

