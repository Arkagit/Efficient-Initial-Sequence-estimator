#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

arma :: vec rrinvgauss(int n, double mu, double lambda){
  
  arma ::vec random_vector(n);
  double z,y,x,u;
  if(mu > 1e6)
  { random_vector = pow(Rcpp::rnorm(1, 0, 1), -2);}
  else{
    for(int i=0; i<n; ++i){
      z=R::rnorm(0,1);
      y=z*z;
      x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
      u=R::runif(0,1);
      if(u <= mu/(mu+x)){
        random_vector(i)=x;
      }else{
        random_vector(i)=mu*mu/x;
      };
    }
  }
  return(random_vector);
}


arma::vec reptest2(arma ::vec x, arma::vec y) {
  int n = y.size();
  arma::vec myvector(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    int p = y[i];
    std::fill(myvector.begin()+ind, myvector.begin()+ind+p, x[i]);
    ind += p;
  }
  return myvector;
}


arma::vec backsolve_c(arma::mat U, arma::vec b){
  int n = b.n_elem;
  arma::vec x = zeros(n);
  for(int i = 0; i<n; i++)
  {
    x(n-1-i) = (b(n-1-i) - dot(U.row(n-1-i), x)) / U((n-1-i), (n-1-i));
  }
  return(x);
}


arma::vec forwardsolve_c(arma::mat U, arma::vec b){
  int n = b.n_elem;
  arma::vec x = zeros(n);
  for(int i = 0; i<n; i++)
  {
    x(i) = (b(i) - dot(U.row(i), x)) / U((i), (i));
  }
  return(x);
}

arma::vec iter_bgl_original(vec beta, double sigma2, arma ::mat X, arma ::vec Y, double alpha,
                            double xi, arma::vec XTY, int n, int p, double lambda, 
                            arma::vec group_size)
{
  int k = group_size.n_elem;
  arma ::vec beta_l2(k);
  int start = 0;
  int  end;
  for(int i = 0; i < k; i++){
    end = start + group_size(i) - 1;
    vec beta_i;
    beta_i = beta.subvec(start, end);
    beta_l2(i) = sum(square(beta_i));
    start = start + group_size(i);
  }
  arma :: vec tau_inv(k);
  for(int i = 0; i < k; i++){
    double mu, betai;
    mu = sqrt(pow(lambda, 2) * sigma2 / beta_l2(i));
    betai = as_scalar(rrinvgauss(1, mu, pow(lambda, 2)));
    tau_inv(i) = betai;
  }
  arma :: vec tau_inv_f(p);
  tau_inv_f = reptest2(tau_inv, group_size);
  arma :: vec tau = 1 / tau_inv_f;
  arma ::vec u(p);
  arma ::vec v_norm(n);
  for(int i = 0; i < p; i++){
    u(i) = sqrt(tau(i)) * (R::rnorm(0,1) );
  }
  for(int i = 0; i < n; i++){
    v_norm(i) = R::rnorm(0, 1);
  }
  arma ::vec v(n);
  v = X * u + v_norm;
  arma ::mat mat_t;
  mat_t = arma::eye<arma::mat>(n,n) + X * diagmat(tau) * X.t();
  arma ::mat U =chol(mat_t);
  arma :: vec w(n);
  w = backsolve_c(U, forwardsolve_c(U.t(), Y / sqrt(sigma2) - v));
  
  arma ::vec beta_new(p);
  beta_new = sqrt(sigma2) * (u + tau % (X.t() * w));
  double sigma2_new, part1, part2;
  part1 = sum(square(Y - X * beta_new));
  part2 = sum(square(beta_new) % tau_inv_f ) ;
  sigma2_new = (part1 + part2 + 2*xi)/ (R::rchisq(n+p + 2*alpha));
  arma::vec res((p + 1));
  res.subvec(0, (p-1)) = beta_new;
  res(p) = sigma2_new;
  return(res);
}


arma::vec iter_bgl_fast(vec beta, double sigma2, arma ::mat X, arma ::vec Y, double alpha,
                        double xi, arma::vec XTY, int n, int p, double lambda, 
                        arma::vec group_size)
{
  int k = group_size.n_elem;
  arma ::vec beta_l2(k);
  int start = 0;
  int  end;
  for(int i = 0; i < k; i++){
    end = start + group_size(i) - 1;
    vec beta_i;
    beta_i = beta.subvec(start, end);
    beta_l2(i) = sum(square(beta_i));
    start = start + group_size(i);
  }
  arma :: vec tau_inv(k);
  for(int i = 0; i < k; i++){
    double mu, betai;
    mu = sqrt(pow(lambda, 2) * sigma2 / beta_l2(i));
    betai = as_scalar(rrinvgauss(1, mu, pow(lambda, 2)));
    tau_inv(i) = betai;
  }
  arma :: vec tau_inv_f(p);
  tau_inv_f = reptest2(tau_inv, group_size);
  arma :: vec tau = 1 / tau_inv_f;
  arma ::vec u(p);
  arma ::vec v_norm(n);
  for(int i = 0; i < p; i++){
    u(i) = sqrt(tau(i)) * (R::rnorm(0,1) );
  }
  for(int i = 0; i < n; i++){
    v_norm(i) = R::rnorm(0, 1);
  }
  arma ::vec v(n);
  v = X * u + v_norm;
  arma ::mat mat_t;
  mat_t = arma::eye<arma::mat>(n,n) + X * diagmat(tau) * X.t();
  arma ::mat U =chol(mat_t);
  arma :: vec w(n);
  w = backsolve_c(U, forwardsolve_c(U.t(),  -v));
  
  arma ::vec beta_tilde;
  beta_tilde = X.t() * backsolve_c(U, forwardsolve_c(U.t(), Y));
  beta_tilde = tau % beta_tilde;
  double sigma2_new, part1, part2;
  part1 = sum(square(Y)) - sum(XTY % beta_tilde);
  part2 = R::rchisq(n + 2*alpha);
  sigma2_new = part1 / part2;
  arma ::vec beta_new(p);
  beta_new = beta_tilde + sqrt(sigma2_new) * (u + tau % (X.t() * w));
  
  arma::vec res((p + 1));
  res.subvec(0, (p-1)) = beta_new;
  res(p) = sigma2_new;
  return(res);
}

//' Bayesian group Lasso.
//' 
//' Inference for the group lasso model by two block Gibbs sampling from the 
//' Bayesian posterior distribution.
//' 
//' @param X matrix.
//' @param Y vector of response.
//' @param group_size Integer vector representing the size of the groups of 
//' the design matrix.
//' @param beta initial setting of the regression coefficients.
//' @param sigam2 initial setting of sigma square. The default setting for 
//' sigma2 is 1.
//' @param lambda group lasso penalty parameter. Default is 1.
//' @param alpha shape parameter for the prior of sigma square.
//'  The default setting for alpha is 0.
//' @param xi rate parameter for the prior of sigma squaresquare. 
//' The default setting for xi is 0.
//' @param K total number of MCMC samplers to be collected. The default is 10000.
//' @return A list with the following elements:
//' \itemize{
//' \item{betas}{a K by p matrix contains mcmc samples of beta.}
//' \item{sigma2s}{a vector contains mcmc sample of sigma2.}
//' }
//' @examples{
//' n = 100
//' p = 20
//' r <- .5
//' Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
//' X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
//' X <- scale(X.raw)*sqrt(n/(n-1))
//'  X <- matrix(as.vector(X),n,p)
//'  beta_t <- c(.3, -1, 0, .5, .01, rep(0, 5), rep(.8, 5), rep(0, 5))
//'  Y.raw <- drop(X%*%beta_t+rnorm(n,0,3))
//'  Y <- Y.raw-mean(Y.raw)
//'  fit_bgl <- bgl(X, Y, rep(5,4), rep(1, p), 1, 1)

//' }
//' 
// [[Rcpp::export]]
List bgl(arma::mat X, arma::vec Y, arma::vec group_size, arma::vec beta, 
         double sigma2, double lambda =1, double alpha = 0, double xi = 0,
         int K = 10000){
  arma ::vec XTY = X.t() * Y;
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat betas(K, p);
  arma::vec sigma2s(K);
  arma::vec iter;
  for(int i = 0; i<K; i++){
    iter = iter_bgl_fast(beta, sigma2, X, Y, alpha, xi, XTY, n, p, lambda, group_size);
    beta = iter.subvec(0, (p-1));
    sigma2 = iter(p);
    sigma2s(i) = sigma2;
    betas.row(i) = beta.t();
  }
  List chains = List::create(Rcpp ::Named("betas") = betas, 
                             Rcpp ::Named("sigma2s")= sigma2s);
  return(chains);
}

//' Bayesian group Lasso.
//' 
//' Inference for the group lasso model by three step Gibbs sampling from 
//' the Bayesian posterior distribution.
//' 
//' @param X A matrix respresenting the design matrix of the linear regression model.
//' @param Y vector of response.
//' @param group_size Integer vector representing the size of the groups of 
//' the design matrix.
//' @param beta initial setting of the regression coefficients.
//' @param sigam2 initial setting of sigma square. The default setting for 
//' sigma2 is 1.
//' @param lambda group lasso penalty parameter.
//' @param alpha shape parameter for the prior of sigma square.
//'  The default setting for alpha is 0.
//' @param xi rate parameter for the prior of sigma squaresquare. 
//' The default setting for xi is 0.
//' @param K total number of MCMC samplers to be collected. Default is 10000.
//' @return A list with the following elements:
//' \itemize{
//' \item{betas:}{ a K by p matrix contains mcmc samples of beta.}
//' \item{sigma2s:}{ a vector contains mcmc sample of sigma2.}
//' }
//' @examples{
//' n = 100
//' p = 20
//' r <- .5
//' Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
//' X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
//' X <- scale(X.raw)*sqrt(n/(n-1))
//'  X <- matrix(as.vector(X),n,p)
//'  beta_t <- c(.3, -1, 0, .5, .01, rep(0, 5), rep(.8, 5), rep(0, 5))
//'  Y.raw <- drop(X%*%beta_t+rnorm(n,0,3))
//'  Y <- Y.raw-mean(Y.raw)
//'  fit_bgl <- bgl_original(X, Y, rep(5,4), rep(1, p), 1, 1)

//' }
//' 
// [[Rcpp::export]]
List bgl_original(arma::mat X, arma::vec Y, arma::vec group_size, arma::vec beta, 
                  double sigma2, double lambda =1, double alpha = 0, double xi = 0,
                  int K = 10000){
  arma ::vec XTY = X.t() * Y;
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat betas(K, p);
  arma::vec sigma2s(K);
  arma::vec iter;
  for(int i = 0; i<K; i++){
    iter = iter_bgl_original(beta, sigma2, X, Y, alpha, xi, XTY, n, p, lambda, group_size);
    beta = iter.subvec(0, (p-1));
    sigma2 = iter(p);
    sigma2s(i) = sigma2;
    betas.row(i) = beta.t();
  }
  List chains = List::create(Rcpp ::Named("betas") = betas, 
                             Rcpp ::Named("sigma2s")= sigma2s);
  return(chains);
}

