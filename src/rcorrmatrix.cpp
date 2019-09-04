// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include <math.h> 

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

// Implementation borrowed from cPCG package
arma::vec cgsolve(arma::mat A, arma::vec b, float tol = 1e-6, int maxIter = 1000) {
  /* Function for solving linear equations Ax = b using conjugate gradient
  !!!todo: preconditioning,c++
  Input:
  A: matrix.
  b: vector
  Output
  x: vector
  */
  // get number of columns of A
  int C = A.n_cols ;
  int R = A.n_rows ;
  // initiate solution x as zeros
  arma::vec x(C) ;
  x.zeros() ; 
  arma::vec oneVec(C);
  oneVec.ones() ;
  
  arma::vec r = b - A * x;
  arma::vec p = r;
  double rs_old = as_scalar( r.t() * r );
  double rs_new=1;
  arma::vec rs_ratio(1);
  
  arma::vec Ap(R);
  double alpha;
  // vector version of alpha
  arma::vec alphaVec(1);
  
  for(int iter = 0; (iter < maxIter) && (rs_new > tol); iter++){
    Ap = A * p;
    alpha = rs_old / as_scalar(p.t() * Ap);
    alphaVec.fill(alpha); 
    x += (oneVec * alphaVec) % p;
    r -= (oneVec * alphaVec) % Ap;
    rs_new = as_scalar( r.t() * r );
    rs_ratio.fill(rs_new / rs_old); 
    
    p = r + (oneVec * rs_ratio) % p;
    rs_old = rs_new;
    if (iter >= maxIter){
      Rcout << "cg did not converge." << endl;
    }
  }
  
  return x;
  
} 

arma::mat init_and_fill_in_off_diagonal(double alp, int d){
  
  arma::mat R(d, d);
  for (int i=0; i < d - 1; i++){
    R(i, i) = 1;
    double rcorr = 2*R::rbeta(alp, alp) - 1;
    R(i, i+1) = rcorr;
    R(i+1, i) = rcorr;
  }
  R(d-1, d-1) = 1;
  return(R);
}

double rjm(arma::mat rsub, double alp, double tol, int maxiter_cg){
  
  int b = rsub.n_rows;
  
  arma::vec r1 = rsub.submat(1, 0, b-2, 0);
  arma::vec r3 = rsub.submat(1, b-1, b-2, b-1);
  arma::mat R2 = rsub.submat(1, 1, b-2, b-2);
  arma::vec m1 = cgsolve(R2, r1, tol, maxiter_cg);
  arma::vec m3 = cgsolve(R2, r3, tol, maxiter_cg);  
  

  double rcond = 2*R::rbeta(alp, alp) - 1;
  double tem13 = dot(r1, m3);
  double tem11 = dot(r1, m1);
  double tem33 = dot(r3, m3);
  return(tem13 + rcond*sqrt((1 - tem11)*(1 - tem33)));
}

arma::mat fill_in_rest(arma::mat rr, double alpha_d, double tol, int maxiter_cg){
  
  int d = rr.n_rows;

  for (int m =1; m < d-1; m++){
    for (int j =0; j <= d-m -2; j++){
      arma::mat rsub = rr.submat(j, j, j + m + 1, j + m + 1);
      double alp = alpha_d + (d - 1 - (m + 1))/2.0;
      rr(j, j + m + 1) = rjm(rsub, alp, tol, maxiter_cg);
      rr(j + m + 1, j) = rr(j, j + m + 1);
    }    
  }
  return(rr);
}

arma::mat rcorrmatrix_cpp_(double alphad, int d, float tol = .0000001, int maxiter_cg=1000){
  
  arma::mat rr = init_and_fill_in_off_diagonal(alphad + (d-2)/2.0, d);
  rr = fill_in_rest(rr, alphad, tol, maxiter_cg);
  return(rr);
  
}
