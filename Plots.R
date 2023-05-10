source("Gen_data.R")

# datasize
n = 1e5
data = data[1:n,]

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
new_v13 = matrix(c(var_est[1,p],var_est[1,p],var_est[p,1],var_est[p,p]),
                 byrow = TRUE, nrow = 2)
pl3 = ellipse(new_v13, centre = c(0,0), which = c(1, 2), npoints = 100)
plot(pl3, type = "l", col = "blue", main = "Confidence Ellipsoid",
     xlab = "1st Component", ylab = "pth Component")

ise13 = matrix(c(IS$cov[1,1],IS$cov[1,p],IS$cov[p,1],IS$cov[p,p]),
               byrow = TRUE, nrow = 2)
pl2 = ellipse(ise13, centre = c(0,0), which = c(1, 2), npoints = 100)
lines(pl2, type = "l", col = "green", main = "Confidence Ellipsoid")

b13 = matrix(c(BME$cov[1,1],BME$cov[1,p],BME$cov[p,1],BME$cov[p,p]),
             byrow = TRUE, nrow = 2)
pl1 = ellipse(b13, centre = c(0,0), which = c(1, 2), npoints = 100)
lines(pl1, type = "l", col = "red", main = "Confidence Ellipsoid")

tr13 = matrix(c(sig[1,1], sig[1,p], sig[p,1], sig[p,p]),
              byrow = TRUE, nrow = 2)
pl0 = ellipse(tr13, centre = c(0,0), which = c(1, 2), npoints = 100)
lines(pl0, type = "l", col = "black", main = "Confidence Ellipsoid")

legend(x = "topleft", lty = c(1,1,1,1), text.font = 1, 
       col= c("blue","green", "red", "black"),text.col = "blue", 
       legend=c("New", "ISE", "Lugsail BM", "True"))

