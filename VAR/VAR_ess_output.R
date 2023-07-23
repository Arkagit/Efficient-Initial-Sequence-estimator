source("VAR_ess.R")

p <- 10
phi <- diag(rep(.95,p))
omega <- matrix(.9, nrow = p, ncol = p)
diag(omega) <- 1
omega <- omega^(abs(col(omega)-row(omega)))

truth <- true.sig(p = p, omega = omega, phi = phi)
ess_true <- (det(truth$tar.var)/det(truth$final.cov))^(1/p)
ess_track <- var_track(N = 5e4, phi = phi, omega = omega, B = 10)

nloops <- 50
subsize <- seq(5e3, 5e4, length = nloops)
# plot(subsize, colMeans(ess_track$BM), type= 'l', ylim = range(ess_track))
# lines(subsize, colMeans(ess_track$WBM), col = "red")
# lines(subsize, colMeans(ess_track$LUG), col = "blue")
# abline(h = ess_true)
save(subsize, ess_true, ess_track, file = "ESS_data.Rdata")

