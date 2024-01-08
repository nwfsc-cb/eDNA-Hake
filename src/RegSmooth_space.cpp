//#define TMB_LIB_INIT R_init_hake
#include <TMB.hpp>
#include "utils.h"

using namespace R_inla;
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

// function for calculating spde tau from kappa, fixing sigma_omega at 1.
template<class Type>
Type calc_tau(Type ln_kappa) {
  Type tau = 1 / sqrt(Type(4.) * M_PI * exp(Type(2.) * ln_kappa));
  return tau;
}


// ------------------ Main TMB template ----------------------------------------
// 
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(Y_i);      // response
  //DATA_VECTOR();      // log offset for volume filtered
  
  // observation indexes
  DATA_IVECTOR(year_idx);  // integer index for years
  DATA_IVECTOR(depth_idx); // integer index for discrete depths.
  
  // Fixed Effects Matrix
  DATA_MATRIX(Xf);   // model matrix for fixed effects
  // int n_i = y_i.rows();   // number of observations
  
  DATA_STRUCT(Zs, LOM); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix
  
  //DATA_STRUCT(proj_Zs, sdmTMB_AOS::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  //DATA_MATRIX(proj_Xs); // smoother linear effect matrix
  
  // optional stuff for penalized regression splines
  DATA_INTEGER(has_smooths);  // whether or not smooths are included
  DATA_IVECTOR(b_smooth_start); /// allows for different basis function numbers among smoothes.
  DATA_INTEGER(n_smooth) ;
  
  DATA_IVECTOR(b_L_smooth_start);
  //DATA_INTEGER(n_L_smooth) ;  
  
  // SPDE objects from R-INLA
  //DATA_STRUCT(spde_aniso, spde_aniso_t);
  DATA_STRUCT(spde, spde_t);
  //PARAMETER_ARRAY(ln_H_input);
  //DATA_INTEGER(anisotropy);
  
  //DATA_IVECTOR(n_fy); /// vector indicating the year for each factor
  DATA_INTEGER(n_i);  /// number of total observations
  DATA_INTEGER(n_f);  /// number of factors per year for weight matrix L. (currently)
  DATA_INTEGER(n_d);  /// number of distinct depths in model
  DATA_INTEGER(n_y);  /// number of years in model
  DATA_INTEGER(n_s);  /// number of knot locations (constant across )
  DATA_INTEGER(n_fy);  /// number of factor-year combinations
  
  // Weight matrix smooth design matrices
  DATA_STRUCT(Z_L, LOM); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_STRUCT(X_L, LOM); // smoother linear effect matrix
  
  // factor indexes for weight matrix 
  DATA_IVECTOR(F_L_idx);  // integer index for years (of length n_fy)
  DATA_IVECTOR(Y_L_idx); // integer index for years (of length n_fy)
  DATA_IVECTOR(FY_start_idx); // integer index for years (of length n_fy)
  
  DATA_SPARSE_MATRIX(A_st); // INLA 'A' projection matrix for unique stations
  DATA_IVECTOR(A_ID_idx); // Vector of station ID to match up observation and A_st output
  
  
  
  // Projections
  // DATA_SPARSE_MATRIX(proj_mesh);
  // DATA_STRUCT(proj_X_ij, sdmTMB::LOM_t);
  // DATA_MATRIX(proj_X_rw_ik);
  // DATA_FACTOR(proj_year);
  // DATA_MATRIX(proj_z_i);
  // DATA_IVECTOR(proj_spatial_index);
  // 
  // DATA_IVECTOR(spatial_only); // !spatial_only means include spatiotemporal(!)
  // DATA_INTEGER(spatial_covariate); // include SVC?
  
  ////////////////////////////////////////////////////////////////
  // Fixed effects
  PARAMETER_VECTOR(betaf); // fixed coefficients.
  
  // Smooth Effects
  PARAMETER_VECTOR(bs); // smoother linear effects
  PARAMETER_VECTOR(bs_L); // smoother linear effects for weight matrix
  
  // PARAMETER_VECTOR(ln_tau_O);    // spatial process
  // PARAMETER_ARRAY(ln_tau_Z);    // optional spatially varying covariate process
  // PARAMETER_VECTOR(ln_tau_E);    // spatio-temporal process
  // PARAMETER_ARRAY(ln_kappa);    // Matern parameter
  
  PARAMETER(lnSigma);    // log- observation sdd
  
  // Random effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included

  PARAMETER_VECTOR(b_L_smooth);  // P-spline smooth parameters for L matrix
  PARAMETER_VECTOR(ln_L_smooth_sigma);  // variances of spline L REs if included
  
  //PARAMETER_VECTOR(ln_tau_O);    // spatial process
  PARAMETER_VECTOR(ln_kappa);    // Matern parameter
  
  //Random effects
  PARAMETER_ARRAY(omega_s);    // spatial effects; years x number of knots x number of factors 
  
  /////// END PARAMETERS /////////
  
  // Joint Negative log-likelihood.
  Type jnll = 0. ;
  

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

   matrix<Type>  eta_smooth_i(n_i,n_smooth);
   eta_smooth_i.setZero();
   vector<Type>  eta_smooth_all;
   eta_smooth_all.setZero() ;
   
if(has_smooths){  
  for(int k = 0; k < n_smooth; k++){
    vector<Type> beta_s(Zs(k).cols());
    beta_s.setZero();
         for (int j = 0; j < beta_s.rows(); j++) {
           beta_s(j) = b_smooth(b_smooth_start(k) + j);
           PARALLEL_REGION jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(k)), true);
         }
    eta_smooth_i.col(k) = Zs(k) * beta_s ;
    //eta_smooth_i = eta_smooth_i(ID+1) * beta_s ;
  }
  // eta_smooth_i += Zs(ID) * beta_s ; // iterate over the Zs list of smooth matrices
  // eta_smooth_i += Xs * bs; // only a single Xs for all of the smoothes

  //vector<Type> eta_smooth_all;
  eta_smooth_all = eta_smooth_i.rowwise().sum();

  //REPORT(b_smooth);     // smooth coefficients for penalized splines
  //REPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
  REPORT(b_smooth);     // smooth coefficients for penalized splines
  REPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
}

// ------------------ Geospatial ---------------------------------------------

// Matern Range calculations
array<Type> range(1,n_f);
  for (int f = 0; f < n_f; f++) {
    range(f) = sqrt(Type(8.)) / exp(ln_kappa(f));
  }
ADREPORT(range);

// calculate tau for each factor from kappa assuming sigma == 1
array<Type> tau_O(n_f);
for(int f=0; f < n_f; f++){
  tau_O(f) = calc_tau(ln_kappa(f)) ;
}
REPORT(tau_O);

// Make spde model. 
Eigen::SparseMatrix<Type> Q_f; // Precision matrix for factors (one for each factor)

// f = factorcombination, supplied from data (looks like 1,1,1,2,2,2,3,3,3, etc.)
  for(int fy = 0; fy < n_fy; fy++){ // loop over 
    Q_f = R_inla::Q_spde(spde, exp(ln_kappa(F_L_idx(fy))));
    PARALLEL_REGION jnll += SCALE(GMRF(Q_f), 1. / tau_O(F_L_idx(fy)))(omega_s.col(fy));
  }
      // if (sim_re(0)) {
      //   vector<Type> omega_s_tmp(omega_s.rows());
      //   SIMULATE {
      //     GMRF(Q_temp, s).simulate(omega_s_tmp);
      //     omega_s.col(m) = omega_s_tmp / exp(ln_tau_O(m));
      //   }
      // }
  //}

  // Make spatial factor matrix
  vector<matrix<Type>> L(n_y); // vector of weight matrices

  /// make L matrix from smoothes.
  for(int y = 0; y < n_y; y++){
    matrix<Type> TMP_mat(n_d,n_f) ;
    TMP_mat.setZero();
    for(int f=0; f < n_f;f++){
      vector<Type> beta_L(Z_L(y).cols());
      beta_L.setZero();
      for (int j = 0; j < beta_L.rows(); j++) {
        beta_L(j) = b_L_smooth(b_L_smooth_start(y) + j);
        PARALLEL_REGION jnll -= dnorm(beta_L(j), Type(0), exp(ln_L_smooth_sigma(f)), true);
      }
      L(y).col(f) = vector<Type>(X_L(y) * bs_L(y)) + Z_L(y) * beta_L ; 
    } // end f loop
  } //end y loop
  // 

  // Make a epsilon matrix for each year with n_d (# depth_cat) columns, n_s (# knots) rows.
  vector<matrix<Type>> epsilon_s(n_y);
   for(int y = 0; y < n_y; y++){
      matrix<Type> epsilon_temp(n_s,n_d);
      epsilon_temp.setZero() ;
      matrix<Type>omega_s_temp(n_s,n_f)  ;
      omega_s_temp =  omega_s.block(0,FY_start_idx(y),n_s,n_f) ;
      epsilon_temp = L(y) * omega_s_temp.transpose() ;
      epsilon_s(y) = epsilon_temp.transpose().matrix() ;
  }
  // 
  // // ------------------ INLA projections ---------------------------------------
  // 
  // Here we are projecting the spatiotemporal and spatial random effects to the
  // locations of the data using the INLA 'A' matrices.

  // Make the spatial projection from knots to station locations.
  vector<matrix<Type>> epsilon_i(n_y);
  for(int y = 0; y < n_y; y++){
    matrix<Type>epsilon_temp(n_i,n_d) ;
    for(int d=0; d<n_d; d++){
      epsilon_temp.col(d) = A_st * epsilon_s(y).col(d) ;
    }
    epsilon_i(y) = epsilon_temp;
  }

  vector<Type> D_i(n_i); // latent variable at the year-station-depth level for each observation.
  D_i =  Xf*betaf + Xs*bs + eta_smooth_all  ;  // all fixed and smooth effects.
   for(int i=0; i < n_i; i++){
     D_i(i) = D_i(i) + epsilon_i(year_idx(i))(A_ID_idx(i),depth_idx(i)); // Add spatial effects.
   }
  // 

  jnll -= sum(dnorm(Y_i,  D_i , exp(lnSigma),true));
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
