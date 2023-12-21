//#define TMB_LIB_INIT R_init_hake
#include <TMB.hpp>
#include "utils.h"

//using namespace R_inla;
using namespace density;
using namespace Eigen;
using namespace tmbutils;

// List of matrices
template <class Type>
struct LOM : vector<matrix<Type> > {
  LOM(SEXP x) {  // x = list passed from R
    (*this).resize(LENGTH(x));
    for (int i = 0; i < LENGTH(x); i++) {
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};


// ------------------ Main TMB template ----------------------------------------
// 
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(Y_i);      // response
  
  // Fixed Effects Matrix
  DATA_MATRIX(Xf);   // model matrix for fixed effects
  // int n_i = y_i.rows();   // number of observations
  
  DATA_STRUCT(Zs, LOM); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix
  
  //DATA_STRUCT(proj_Zs, sdmTMB_AOS::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  //DATA_MATRIX(proj_Xs); // smoother linear effect matrix
  
  // optional stuff for penalized regression splines
  DATA_INTEGER(has_smooths);  // whether or not smooths are included
  DATA_IVECTOR(b_smooth_start); /// allows for different basis function nubmers among smoothes.
  DATA_INTEGER(n_smooth) ;
  
  // Fixed effects
  PARAMETER_VECTOR(betaf); // fixed coefficients.
  
  // Smooth Effects
  PARAMETER_VECTOR(bs); // smoother linear effects
  
  // PARAMETER_VECTOR(ln_tau_O);    // spatial process
  // PARAMETER_ARRAY(ln_tau_Z);    // optional spatially varying covariate process
  // PARAMETER_VECTOR(ln_tau_E);    // spatio-temporal process
  // PARAMETER_ARRAY(ln_kappa);    // Matern parameter
  
  PARAMETER(lnSigma);    // log- observation sdd
  
  // Random effects
 PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
 //PARAMETER_VECTOR(beta_s);
 PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included
  
  /////// END PARAMETERS /////////
  
  // Joint Negative log-likelihood.
  Type jnll = 0. ;
  
  /////////////
  int n_i = Y_i.rows();   // number of observations
  // int n_t = ;   // number of years

  // Spline penalties
  // p-splines/smoothers
  // vector<Type>  eta_smooth_i; //(n_i);
  // eta_smooth_i.setZero();
  // if (has_smooths) {
  //     for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
  //       vector<Type> beta_s(Zs(s).cols());
  //       beta_s.setZero();
  //       for (int j = 0; j < beta_s.rows(); j++) {
  //         beta_s(j) = b_smooth(b_smooth_start(s) + j);
  //         PARALLEL_REGION jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
  //         //if (sim_re(5)) SIMULATE{beta_s(j) = rnorm(Type(0), exp(ln_smooth_sigma(s)));}
  //       }
  //       eta_smooth_i += Zs(s) * beta_s ; // iterate over the Zs list of smooth matrices
  //     }
       // eta_smooth_i += Xs * bs; // only a single Xs for all of the smoothes
  //   
  //   REPORT(b_smooth);     // smooth coefficients for penalized splines
  //   REPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
  // }
  
// 
   matrix<Type>  eta_smooth_i(n_i,n_smooth);
   eta_smooth_i.setZero();

  for(int s = 0; s < n_smooth; s++){
    vector<Type> beta_s(Zs(s).cols());
    beta_s.setZero();
         for (int j = 0; j < beta_s.rows(); j++) {
           beta_s(j) = b_smooth(b_smooth_start(s) + j);
           PARALLEL_REGION jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(s)), true);
         }
    eta_smooth_i.col(s) = Zs(s) * beta_s ;
    //eta_smooth_i = eta_smooth_i(ID+1) * beta_s ;
  }
  // eta_smooth_i += Zs(ID) * beta_s ; // iterate over the Zs list of smooth matrices
  // eta_smooth_i += Xs * bs; // only a single Xs for all of the smoothes

  vector<Type> eta_smooth_all;
  eta_smooth_all = eta_smooth_i.rowwise().sum();
  vector<Type> Di; // latent variable at the station-depth level.
  Di =  Xf*betaf + Xs*bs + eta_smooth_all ;  // 
  
  jnll -= sum(dnorm(Y_i,  Di , exp(lnSigma),true));
  return jnll ;
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
