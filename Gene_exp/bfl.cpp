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


arma::vec rep(arma ::vec x, arma::vec y) {
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

arma::mat tridiag(arma::vec upper, arma::vec lower, arma::vec main){
  int dim = main.n_elem;
  arma::mat tri_diag(dim, dim); 
  tri_diag = zeros(dim, dim);
  tri_diag.diag() = main;
  int k = upper.n_elem;
  for( int i=0; i<k; i++){
    tri_diag((i + 1), i) = lower(i);
    tri_diag(i, (i + 1)) = upper(i);
  }
  return(tri_diag);
}

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat U) {
  int ncols = U.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * U;
}


arma::vec iter_bfl_original(vec beta, double sigma2, arma ::mat X, arma ::vec Y, double alpha,
                       double xi, arma::vec XTY, arma::mat XTX, int n, int p, 
                       double lambda1, double lambda2){
  arma::vec tau_inv(p);
  for(int i = 0; i < p; i++){
    double mu, betai;
    mu = sqrt(pow(lambda1, 2) * sigma2 / pow(beta(i), 2));
    betai = as_scalar(rrinvgauss(1, mu, pow(lambda1,2)));
    tau_inv(i) = betai;
  }
  arma :: vec omega_inv(p-1), beta_diff(p-1);
  beta_diff = beta.subvec(1, p-1) - beta.subvec(0, p-2);
  for(int i = 0; i <(p-1); i++)
  {
    double mu, betai;
    mu = sqrt(pow(lambda2, 2) * sigma2 / pow(beta_diff(i),2));
    betai = as_scalar(rrinvgauss(1, mu, pow(lambda2, 2)));
    omega_inv(i) = betai;
  }
  arma ::vec part1(p), part2(p);
  part1 = zeros(p);
  part2 = zeros(p);
  part1.subvec(0, (p-2)) = omega_inv;
  part2.subvec(1, (p-1)) = omega_inv;
  arma ::vec diag_omega = part1 + part2;
  arma ::vec diag_sum = tau_inv + diag_omega;
  arma ::mat inv_omega_tau =  tridiag((-1)*omega_inv, (-1)*omega_inv, diag_sum);
  arma ::mat A_tau_omega = XTX  + inv_omega_tau;
  arma ::mat U = chol(A_tau_omega);
  arma ::vec b_mean = backsolve_c(U, forwardsolve_c(U.t(), XTY));
  arma ::mat beta_new;
  arma ::mat V(p,p);
  arma ::mat I = arma::eye<arma::mat>(p,p);
  for(int i = 0; i<p; i++){
    V.col(i) = backsolve_c(U, I.col(i));
  }
  beta_new = mvrnormArma(1, b_mean, sqrt(sigma2) * V.t()).t();
  double rate, sigma2_new;
  rate = sum(square(Y-X*beta_new))+as_scalar(beta_new.t()*inv_omega_tau*beta_new)+2*xi;
  sigma2_new = rate / R::rchisq(n + p + 2*alpha);
  arma::vec res((p + 1));
  res.subvec(0, (p-1)) = beta_new;
  res(p) = sigma2_new;
  return(res);
}


arma::vec iter_bfl_fast(vec beta, double sigma2, arma ::mat X, arma ::vec Y, double alpha,
                            double xi, arma::vec XTY, arma::mat XTX, int n, int p, 
                            double lambda1, double lambda2){
  arma::vec tau_inv(p);
  for(int i = 0; i < p; i++){
    double mu, betai;
    mu = sqrt(pow(lambda1, 2) * sigma2 / pow(beta(i), 2));
    betai = as_scalar(rrinvgauss(1, mu, pow(lambda1,2)));
    tau_inv(i) = betai;
  }
  arma :: vec omega_inv(p-1), beta_diff(p-1);
  beta_diff = beta.subvec(1, p-1) - beta.subvec(0, p-2);
  for(int i = 0; i <(p-1); i++)
  {
    double mu, betai;
    mu = sqrt(pow(lambda2, 2) * sigma2 / pow(beta_diff(i),2));
    betai = as_scalar(rrinvgauss(1, mu, pow(lambda2, 2)));
    omega_inv(i) = betai;
  }
  arma ::vec part1(p), part2(p);
  part1 = zeros(p);
  part2 = zeros(p);
  part1.subvec(0, (p-2)) = omega_inv;
  part2.subvec(1, (p-1)) = omega_inv;
  arma ::vec diag_omega = part1 + part2;
  arma ::vec diag_sum = tau_inv + diag_omega;
  arma ::mat inv_omega_tau =  tridiag((-1)*omega_inv, (-1)*omega_inv, diag_sum);
  arma ::mat A_tau_omega = XTX  + inv_omega_tau;
  arma ::mat U = chol(A_tau_omega);
  double rate, sigma2_new;
  rate = sum(square(Y))-as_scalar(XTY.t()*backsolve_c(U, forwardsolve_c(U.t(), XTY)))+2*xi;
  sigma2_new = rate / R::rchisq(n + 2*alpha);
  arma :: vec b_mean = backsolve_c(U, forwardsolve_c(U.t(), XTY));
  arma ::mat V(p,p);
  arma ::mat I = arma::eye<arma::mat>(p,p);
  for(int i = 0; i<p; i++){
    V.col(i) = backsolve_c(U, I.col(i));
  }
  arma :: mat beta_new = mvrnormArma(1, b_mean, sqrt(sigma2_new) * V.t()).t();
  arma::vec res((p + 1));
  res.subvec(0, (p-1)) = beta_new;
  res(p) = sigma2_new;
  return(res);
}  
// [[Rcpp::export]]
List bfl_original(arma::mat X, arma::vec Y, arma::vec beta, double sigma2, double lambda1 =1, 
                  double lambda2 = 1, double alpha = 0, double xi = 0, int K = 20000){
  arma ::vec XTY = X.t() * Y;
  arma ::mat XTX = X.t() * X;
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat betas(K, p);
  arma::vec sigma2s(K);
  arma::vec iter;
  for(int i = 0; i<K; i++){
    iter = iter_bfl_original(beta, sigma2, X, Y, alpha, xi, XTY, XTX, n, p, lambda1, lambda2);
    beta = iter.subvec(0, (p-1));
    sigma2 = iter(p);
    sigma2s(i) = sigma2;
    betas.row(i) = beta.t();
  }
  List chains = List::create(Rcpp ::Named("betas") = betas, 
                             Rcpp ::Named("sigma2s")= sigma2s);
  return(chains);
}

// [[Rcpp::export]]
List bfl(arma::mat X, arma::vec Y, arma::vec beta, double sigma2, double lambda1 =1, 
                  double lambda2 = 1, double alpha = 0, double xi = 0, int K = 20000){
  arma ::vec XTY = X.t() * Y;
  arma ::mat XTX = X.t() * X;
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat betas(K, p);
  arma::vec sigma2s(K);
  arma::vec iter;
  for(int i = 0; i<K; i++){
    iter = iter_bfl_fast(beta, sigma2, X, Y, alpha, xi, XTY, XTX, n, p, lambda1, lambda2);
    beta = iter.subvec(0, (p-1));
    sigma2 = iter(p);
    sigma2s(i) = sigma2;
    betas.row(i) = beta.t();
  }
  List chains = List::create(Rcpp ::Named("betas") = betas, 
                             Rcpp ::Named("sigma2s")= sigma2s);
  return(chains);
}

