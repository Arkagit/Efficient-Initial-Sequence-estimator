gamma <- function(data, t){
	mu = colMeans(data)
	n = dim(data)[1]
	foo = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
	for(i in 1:(n-t)){
		foo = foo + (data[i,] - mu)%*%t(data[i+t,] - mu)
	}
	foo = (foo + t(foo))/2
	return(foo/n)
}

mise <- function(data){
	n <- dim(data)[1]
	p <- dim(data)[2]
	m <- 0
	s_n <- 0
	Sigma_sn <- gamma(data, 2*s_n) + 2*gamma(data, 2*s_n+1)

	while((min(eigen(Sigma_sn, only.values = TRUE)$values) <= 0) == TRUE){
		s_n = s_n + 1
		Sigma_sn = Sigma_sn + 2*gamma(data, 2*s_n) + 2*gamma(data, 2*s_n+1)
	}

	for(j in (s_n + 1):floor(n/2 - 1)){
		dummy <- Sigma_sn
		trun = j
		Sigma_sn = Sigma_sn + 2*gamma(data, 2*j) + 2*gamma(data, 2*j+1)
		if(det(dummy) > det(Sigma_sn)){
			break
		}
	}
	return(list("Variance" = dummy, "truncation" = trun))
}
