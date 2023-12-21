//#define TMB_LIB_INIT R_init_hake
#include <TMB.hpp>
//#include "utils_AOS.h"

//using namespace R_inla;
using namespace density;
using namespace Eigen;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]
// ------------------ Main TMB template ----------------------------------------
// 
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(Y_i);      // response
  DATA_MATRIX(X_f);   // model matrix for fixed effects
  PARAMETER_VECTOR(beta); //
  PARAMETER(lnSigma);    // log- observation sdd
  // int n_i = y_i.rows();   // number of observations
  vector<Type> mu ;
  mu = X_f * beta;
  Type nll = -sum(dnorm(Y_i, mu, exp(lnSigma),true));
  return nll ;
}

// template<class Type>
// Type objective_function<Type>::operator() ()
// {
//   DATA_VECTOR(Y);
//   DATA_MATRIX(X);
//   PARAMETER_VECTOR(b);
//   PARAMETER(logSigma);
//   Type nll = sum(dnorm(Y, X*b , exp(logSigma), true));
//   return -nll;
// }