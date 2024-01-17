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
Type calc_tau(Type logkappa) {
  Type tau = 1 / sqrt(Type(4.) * M_PI * exp(Type(2.) * logkappa));
  return tau;
}

template<class Type>
Type calc_sigma(Type logtau, Type logkappa) {
  Type sigma = 1. / sqrt(4. * M_PI * exp(2. * logtau) * exp(2. * logkappa));
  return sigma;
}

// https://github.com/hrue/r-inla/blob/devel/r-inla.org/doc/prior/pc.matern.pdf
template <class Type>
Type pc_prior_matern(Type logtau, Type logkappa, Type matern_range,
                     Type matern_SD, Type range_prob, Type SD_prob,
                     int give_log = 0, int share_range = 0, int stan_flag = 0) {
  Type d = 2.;  // dimension
  Type dhalf = d / 2.;
  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(SD_prob) / matern_SD;
  Type range = sqrt(8.) / exp(logkappa);
  Type sigma = 1. / sqrt(4. * M_PI * exp(2. * logtau) * exp(2. * logkappa));
  Type range_ll = log(dhalf) + log(lam1) + log(pow(range, -1. - dhalf)) -
    lam1 * pow(range, -dhalf);
  Type sigma_ll = log(lam2) - lam2 * sigma;
  Type penalty = sigma_ll;
  if (!share_range) penalty += range_ll;
  
  // Note: these signs are + (and different from inst/jacobian-pcprior-tests)
  // because the jnll is accumulated
  if (stan_flag) {
    penalty += log(sqrt(8.)) - log(pow(range, 2.)); // P(sigma)
    Type C = sqrt(exp(lgamma(1. + dhalf)) * pow(4. * M_PI, dhalf));
    penalty += log(C) + logkappa;
  }
  // std::cout << "PC penalty: " << penalty << "\n";
  if (give_log)
    return penalty;
  else
    return exp(penalty);
}

// ------------------ Main TMB template ----------------------------------------
// 
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(Y_i);      // response
  
  DATA_VECTOR(Y_stand);      // positive Ct for standards
  DATA_IVECTOR(Z_stand);      // binomial Ct for standards
  DATA_VECTOR(ln_copies);    // know number of copies for standards
  DATA_IVECTOR(plate_stand_idx); /// index for plate id for standards data frame
  
  DATA_VECTOR(Y_samp);      // positive Ct for samples
  DATA_IVECTOR(Z_samp);      // binomial Ct for samples
  DATA_IVECTOR(plate_samp_idx); /// index for plate id for sample data frame
  
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
  DATA_INTEGER(n_stand);  /// number of total observations in the standards data frame
  DATA_INTEGER(n_f);  /// number of factors per year for weight matrix L. (currently)
  DATA_INTEGER(n_d);  /// number of distinct depths in model
  DATA_INTEGER(n_y);  /// number of years in model
  DATA_INTEGER(n_s);  /// number of knot locations (constant across )
  DATA_INTEGER(n_st);  /// number of unique station locations (among all years / surveys.)
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
  
  // PRIORS
  // Standards
  DATA_VECTOR(std_mu); 
  DATA_VECTOR(std_sigma); 
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
  
  
  // standard parameters
  PARAMETER_VECTOR(std_phi_0) ;
  PARAMETER_VECTOR(std_phi_1) ;
  PARAMETER_VECTOR(std_beta_0) ; 
  PARAMETER_VECTOR(std_beta_1) ;
  PARAMETER_VECTOR(ln_std_sigma) ; // Hyper prior SD for phi_0, phi_1, beta_0, beta_1.
  
  PARAMETER_VECTOR(ln_std_tau); // SD of observations on the positive model component
  
  // Fixed effects
  PARAMETER_VECTOR(betaf); // fixed coefficients.
  
  // Smooth Effects
  PARAMETER_VECTOR(bs); // smoother linear effects
  
  // PARAMETER_VECTOR(ln_tau_O);    // spatial process
  // PARAMETER_ARRAY(ln_tau_Z);    // optional spatially varying covariate process
  // PARAMETER_VECTOR(ln_tau_E);    // spatio-temporal process
  // PARAMETER_ARRAY(ln_kappa);    // Matern parameter
  
  
  // Random effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included

  // PARAMETER_ARRAY(bs_L); // smoother linear effects for weight matrix (year x factor)
  // PARAMETER_ARRAY(b_L_smooth);  // P-spline smooth parameters for L matrix (dim(year,factor,smooth knots))
  // PARAMETER_ARRAY(ln_L_smooth_sigma);  // variances of spline L REs if included (dim(year,factor))
  // 
  // //PARAMETER_VECTOR(ln_tau_O);    // spatial process
  PARAMETER(ln_kappa);    // Matern parameter determining range.
  // 
  // //Random effects
  PARAMETER_MATRIX(omega_s);    // spatial effects; (number of knots x (years by number of factors)

  PARAMETER(lnSigma);    // log- observation sdd
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
  
  REPORT(bs);
  REPORT(b_smooth);     // smooth coefficients for penalized splines
  REPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
  REPORT(eta_smooth_i) ;
  REPORT(eta_smooth_all) ;
}
REPORT(betaf) ;

// ------------------ Geospatial ---------------------------------------------
// 
// Matern Range calculations
Type kappa;
Type range;
  //for (int f = 0; f < n_f; f++) {
    kappa = exp(ln_kappa) ; 
    range = sqrt(Type(8.)) / kappa;
    
  //}
ADREPORT(range);
ADREPORT(exp(ln_kappa));

// calculate tau for each factor from kappa assuming sigma == 1
Type tau_O;
Type siggy;
for(int f=0; f < n_f; f++){
  tau_O = calc_tau(ln_kappa) ;
  siggy = calc_sigma(log(tau_O),ln_kappa);
}
REPORT(tau_O);
REPORT(siggy);


// Make spde model.
Eigen::SparseMatrix<Type> Q_f; // Precision matrix for factors (one for each factor)

Q_f = R_inla::Q_spde(spde, exp(ln_kappa));

// f = factorcombination, supplied from data (looks like 1,1,1,2,2,2,3,3,3, etc.)
  for(int fy = 0; fy < n_fy; fy++){ // loop over
    
    // CANNOT MAKE THE KAPPA LOOP OVER factors.  No idea why it doesn't work!!!!!
    // Q_f = R_inla::Q_spde(spde, exp(ln_kappa(fy)));
    // Q_f = R_inla::Q_spde(spde, exp(ln_kappa(F_L_idx(fy))));
    PARALLEL_REGION jnll += SCALE(GMRF(Q_f), 1. / tau_O)(omega_s.col(fy));
    //PARALLEL_REGION jnll += GMRF(Q_f,true)(omega_s.col(fy));
  }
  REPORT(range);
  REPORT(tau_O);
  REPORT(omega_s) ;
      // if (sim_re(0)) {
      //   vector<Type> omega_s_tmp(omega_s.rows());
      //   SIMULATE {
      //     GMRF(Q_temp, s).simulate(omega_s_tmp);
      //     omega_s.col(m) = omega_s_tmp / exp(ln_tau_O(m));
      //   }
      // }
  //}

  // Make spatial factor matrix
  // vector<matrix<Type>> L_mat(n_y); // vector of weight matrices
  // 
  // /// make L matrix from smoothes.
  // for(int y = 0; y < n_y; y++){
  //   matrix<Type> tmp_mat(n_d,n_f) ;
  //   tmp_mat.setZero();
  //   //vector<Type> A(n_d);
  //   for(int f=0; f < n_f;f++){ // loop over factors in each year
  //     vector<Type> beta_L(Z_L(y).cols());
  //     beta_L.setZero();
  //     
  //     /// THIS IS FOR ALLOWING THE WEIGHT MATRICES TO VARY FROM YEAR TO YEAR.
  //     for (int j = 0; j < beta_L.size(); j++) {
  //       beta_L(j) = b_L_smooth(y,f,j);
  //       PARALLEL_REGION jnll -= dnorm(beta_L(j), Type(0), exp(ln_L_smooth_sigma(y,f)), true);
  //     }
  //     tmp_mat.col(f) = matrix<Type>(vector<Type>(X_L(y) *  bs_L(y,f)) + Z_L(y) * beta_L) ; 
  //     // for (int j = 0; j < beta_L.rows(); j++) {
  //     //   beta_L(j) = b_L_smooth(j);
  //     //   PARALLEL_REGION jnll -= dnorm(beta_L(j), Type(0), exp(ln_L_smooth_sigma(f)), true);
  //     // } // end j loop
  //     // tmp_mat.col(f) = matrix<Type>(vector<Type>(X_L(y) * bs_L) + Z_L(y) * beta_L) ;
  //   } // end f loop
  //   L_mat(y) = tmp_mat ;
  // } //end y loop
  // 
  // REPORT(bs_L) ;
  // REPORT(b_L_smooth) ;
  // REPORT(ln_L_smooth_sigma) ;
  // REPORT(L_mat) ;

  // Make a epsilon matrix for each year with n_d (# depth_cat) columns, n_s (# knots) rows.
  vector<matrix<Type>> epsilon_s(n_y);
  for(int y = 0; y < n_y; y++){
       matrix<Type> epsilon_temp(n_s,n_d);
       epsilon_temp.setZero() ;
       for(int d = 0; d< n_d; d++){
         epsilon_temp.col(d) = omega_s.col(FY_start_idx(y));
       }
      //matrix<Type>omega_s_temp(n_s,n_f)  ;
      epsilon_s(y) = epsilon_temp;
      
      //  epsilon_s(y) = omega_s.block(0,FY_start_idx(y),n_s,n_f) * L_mat(y).transpose() ;
  }
   // REPORT(L_mat(0)); 
   // REPORT(epsilon_s(0)) ;
  //
  // // ------------------ INLA projections ---------------------------------------
  //
  // Here we are projecting the spatiotemporal and spatial random effects to the
  // locations of the data using the INLA 'A' matrices.

  // Make the spatial projection from knots to station locations.
  array<Type> epsilon_st(n_st,n_d,n_y);
  for(int y = 0; y < n_y; y++){
    for(int d=0; d<n_d; d++){
       epsilon_st.col(y).col(d) = (A_st * epsilon_s(y).col(d)).array() ;
    }
  }
  REPORT(epsilon_s(0)) ; // epsilon at each knot location.
  REPORT(epsilon_st) ;  // epsilon at each station.
  
  vector<Type> D_i(n_i); // latent variable at the year-station-depth level for each observation.
  vector<Type> E_i(n_i); // latent variable at the year-station-depth level for each observation.
  D_i =  Xf*betaf + Xs*bs + eta_smooth_all  ;  // all fixed and smooth effects.
   for(int i=0; i < n_i; i++){
     D_i(i) = D_i(i) + epsilon_st(A_ID_idx(i),depth_idx(i),year_idx(i)); // Add spatial effects.
     
   }
   
   // This adds offsets for wash error, dilution offsets, etc.
   E_i(i) = D_i(i) + DILUTION OFFSET + WASH ERROR.
 
 REPORT(D_i);
 REPORT(exp(D_i));
 //
 // vector<Type> D_i(n_i); // latent variable at the year-station-depth level for each observation.
 // D_i =  Xf*betaf + Xs*bs + eta_smooth_all  ;  // all fixed and smooth effects.

  jnll -= sum(dnorm(Y_i,  exp(D_i) , exp(lnSigma),true));

  
  // contribution of the standards to the likeihood
  // section modified slightly from Eric's sdmTMB code.
  
    // Presence-Absence component of model.
    vector<Type> theta_stand(n_stand);
    vector<Type> std_sigma = exp(ln_std_sigma);
    
    // prior for binomial part of the model.
    for(int i = 0; i < std_phi_0.size(); i++){
      jnll -= dnorm(std_phi_0(i), std_mu(0), std_sigma(0), true);//phi_0 ~ normal(2, 2)
      jnll -= dnorm(std_phi_1(i), std_mu(1), std_sigma(1), true);//phi_1 ~ normal(4, 2)
    }
    // Positive component of model.
    for(int i = 0; i < std_beta_0.size(); i++){
      jnll -= dnorm(std_beta_0(i), std_mu(2), std_sigma(2), true);//beta_0 ~ normal(std_mu(2))
      jnll -= dnorm(std_beta_1(i), std_mu(3), std_sigma(3), true);//beta_1 ~ normal()
    }
    
    // Likelihoods:
    vector<Type> kappa_stand(n_stand);    
    for(int i = 0; i < n_stand;; i++){ // likelihood
      theta_stand(i) = std_phi_0(plate_stand_idx(i)) + std_phi_1(plate_stand_idx(i)) * (ln_copies(i));// - stand_offset);
      jnll -= dbinom_robust(Z_stand(i), Type(1), theta_stand(i), true); // likelihood
      if(Z_stand(i) == 1){
        kappa_stand(i) = std_beta_0(plate_stand_idx(i)) + std_beta_1(plate_stand_idx(i)) * (ln_copies(i));// - stand_offset);
        jnll -= dnorm(Y_stand(i), kappa_stand(i), exp(ln_std_tau(0)), true); // likelihood  
      }

      
      // Repeat for unknown Samles 
      //Presence-Absence component of model.
      vector<Type> theta_samp(n_i);
      vector<Type> kappa_samp(n_i);
      vector<Type> std_sigma = exp(ln_std_sigma);
      
      for(int i = 0; i < n_stand;; i++){ // likelihood
        theta_samp(i) = std_phi_0(plate_samp_idx(i)) + std_phi_1(plate_samp_idx(i)) * (D(i));// - stand_offset);
        jnll -= dbinom_robust(Z_samp(i), Type(1), theta_samp(i), true); // likelihood
        if(Z_stand(i) == 1){
          kappa_samp(i) = std_beta_0(plate_samp_idx(i)) + std_beta_1(plate_samp_idx(i)) * (D(i));// - stand_offset);
          jnll -= dnorm(Y_samp(i), kappa_samp(i), exp(ln_std_tau(0)), true); // likelihood  
        }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
  //Prior on bs_L
  // for(int y = 0; y < n_y; y++){
  //   jnll -= dnorm(bs_L(y), Type(0), Type(1), true);
  // }

  
  // MATERN PC PRIOR  
  // jnll -= sdmTMB::pc_prior_matern(
  //   ln_tau_O(m), ln_kappa(0,m),
  //   priors(0), priors(1), priors(2), priors(3),
  //   true, /* log */
  //   false, /* share range */
  //   stan_flag);
  
  
  
  
  
  return jnll ;
}


// Likelihood for standard curves.

// case stdcurve_family: {
// if(y_i(i,m) > 0) {
//   tmp_ll = dnorm(y_i(i,m), std_beta_0(pcr_idx(i)) + std_beta_1(pcr_idx(i)) * mu_i(i,m), phi(m), true);
//   tmp_ll += dbinom_robust(Type(1), size(i), std_phi_0(pcr_idx(i)) + std_phi_1(pcr_idx(i)) * mu_i(i,m), true);
// } else {
//   tmp_ll = dbinom_robust(Type(0), size(i), std_phi_0(pcr_idx(i)) + std_phi_1(pcr_idx(i)) * mu_i(i,m), true);
// }
// SIMULATE{y_i(i,m) = rbinom(size(i), invlogit(std_phi_0(pcr_idx(i)) + std_phi_1(pcr_idx(i)) * mu_i(i,m))) * rnorm(std_beta_0(pcr_idx(i)) + std_beta_1(pcr_idx(i)) * mu_i(i,m), phi(m));}
// break;
// }
//   default:
//   error("Family not implemented.");
//   }
//   tmp_ll *= weights_i(i);
//   jnll_obs(i) -= tmp_ll; // for cross validation
//   jnll -= tmp_ll; // * keep
//   }



// Standard Curve section from Eric.
// if(stdcurve_flag == 1) {
//   // Presence-Absence component of model.
//   vector<Type> theta_stand(N_stand_bin);
//   vector<Type> std_sigma = exp(ln_std_sigma);
//   
//   for(int i = 0; i < std_phi_0.size(); i++){
//     jnll -= dnorm(std_phi_0(i), std_mu(0), std_sigma(0), true);//phi_0 ~ normal(2, 2)
//     jnll -= dnorm(std_phi_1(i), std_mu(1), std_sigma(1), true);//phi_1 ~ normal(4, 2)
//   }
//   for(int i = 0; i < N_stand_bin; i++){ // likelihood
//     theta_stand(i) = std_phi_0(pcr_stand_bin_idx(i)) + std_phi_1(pcr_stand_bin_idx(i)) * (D_bin_stand(i));// - stand_offset);
//     jnll -= dbinom_robust(bin_stand(i), Type(1), theta_stand(i), true); // likelihood
//   }
//   // Positive component of model.
//   for(int i = 0; i < std_beta_0.size(); i++){
//     jnll -= dnorm(std_beta_0(i), std_mu(2), std_sigma(2), true);//beta_0 ~ normal(40,5)
//     jnll -= dnorm(std_beta_1(i), std_mu(3), std_sigma(3), true);//beta_1 ~ normal(-3.32,0.1)
//   }
//   vector<Type> kappa_stand(N_stand_pos);
//   //sigma_all_stand = exp(log_sigma_all_stand);
//   for(int i = 0; i < N_stand_pos; i++){
//     kappa_stand(i) = std_beta_0(pcr_stand_pos_idx(i)) + std_beta_1(pcr_stand_pos_idx(i)) * (D_pos_stand(i));// - stand_offset);
//     jnll -= dnorm(pos_stand(i), kappa_stand(i), phi(0), true); // likelihood
//   }


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
