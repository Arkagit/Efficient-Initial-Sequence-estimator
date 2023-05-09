source("Gen_data.R")

# datasize
n = 1e6
data = data[1:n,]

# Batch means estimation
BM = mcse.multi(data, method = "bm", r = 3, size = "sqroot")
BM$cov

# ISE estimation 
IS = mcse.initseq(data)
IS$cov

# correlation matrix by batch means
BM = mcse.multi(data, method = "bm", r = 1, size = "sqroot")
corrm = cor(BM$cov)
# Standard deviation bt ISE
sd = numeric(0)
for (i in 1:p) {
  sd[i] = sqrt(initseq(data[,i])$var.pos)
}
# new variance estimate
var_est = diag(sd)%*%corrm%*%diag(sd)


# Plots for 1st and 3rd component
new_v13 = matrix(c(var_est[1,1],var_est[1,3],var_est[3,1],var_est[3,3]),
                 byrow = TRUE, nrow = 2)
pl3 = ellipse(new_v13, centre = c(0,0), which = c(1, 2), npoints = 100)
plot(pl3, type = "l", col = "blue", main = "Confidence Ellipsoid")

ise13 = matrix(c(IS$cov[1,1],IS$cov[1,3],IS$cov[3,1],IS$cov[3,3]),
               byrow = TRUE, nrow = 2)
pl2 = ellipse(ise13, centre = c(0,0), which = c(1, 2), npoints = 100)
lines(pl2, type = "l", col = "green", main = "Confidence Ellipsoid")

b13 = matrix(c(BM$cov[1,1],BM$cov[1,3],BM$cov[3,1],BM$cov[3,3]),
             byrow = TRUE, nrow = 2)
pl1 = ellipse(b13, centre = c(0,0), which = c(1, 2), npoints = 100)
lines(pl1, type = "l", col = "brown", main = "Confidence Ellipsoid")

legend(x = "topleft", lty = c(1,1,1), text.font = 1, 
       col= c("blue","green", "brown"),text.col = "blue", 
       legend=c("New", "ISE", "Lugsail"))
