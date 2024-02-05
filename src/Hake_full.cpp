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
  Type sigma = 1. / sqrt(4. * M_PI * exp(Type(2.) * logtau) * exp(Type(2.) * logkappa));
  return sigma;
}

// https://github.com/hrue/r-inla/blob/devel/r-inla.org/doc/prior/pc.matern.pdf
// template <class Type>
// Type pc_prior_matern(Type logtau, Type logkappa, Type matern_range,
//                      Type matern_SD, Type range_prob, Type SD_prob,
//                      int give_log = 0, int share_range = 0, int stan_flag = 0) {
//   Type d = 2.;  // dimension
//   Type dhalf = d / 2.;
//   Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
//   Type lam2 = -log(SD_prob) / matern_SD;
//   Type range = sqrt(8.) / exp(logkappa);
//   Type sigma = 1. / sqrt(4. * M_PI * exp(2. * logtau) * exp(2. * logkappa));
//   Type range_ll = log(dhalf) + log(lam1) + log(pow(range, -1. - dhalf)) -
//     lam1 * pow(range, -dhalf);
//   Type sigma_ll = log(lam2) - lam2 * sigma;
//   Type penalty = sigma_ll;
//   if (!share_range) penalty += range_ll;
//   
//   // Note: these signs are + (and different from inst/jacobian-pcprior-tests)
//   // because the jnll is accumulated
//   if (stan_flag) {
//     penalty += log(sqrt(8.)) - log(pow(range, 2.)); // P(sigma)
//     Type C = sqrt(exp(lgamma(1. + dhalf)) * pow(4. * M_PI, dhalf));
//     penalty += log(C) + logkappa;
//   }
//   // std::cout << "PC penalty: " << penalty << "\n";
//   if (give_log)
//     return penalty;
//   else
//     return exp(penalty);
// }

template <class Type>
matrix<Type> MakeH(vector<Type> x) {
  matrix<Type> H(2, 2);
  H(0, 0) = exp(x(0));
  H(1, 0) = x(1);
  H(0, 1) = x(1);
  H(1, 1) = (1 + x(1) * x(1)) / exp(x(0));
  return H;
}




// ------------------ Main TMB template ----------------------------------------
// 
template <class Type>
Type objective_function<Type>::operator()()
{
  //DATA_VECTOR(Y_i);      // response
  
  DATA_VECTOR(Y_stand);      // positive Ct for standards
  DATA_VECTOR(Z_stand);      // binomial of Ct for standards
  DATA_IVECTOR(Z_stand_int);      // binomial of Ct for standards
  DATA_VECTOR(ln_copies);    // know number of copies for standards
  DATA_IVECTOR(plate_stand_idx); /// index for plate id for standards data frame
  
  DATA_VECTOR(Y_samp);      // positive Ct for samples
  DATA_VECTOR(Z_samp);      // binomial of Ct for samples
  DATA_VECTOR(Z_samp_int);      // binomial of Ct for samples
  
  DATA_IVECTOR(plate_samp_idx); /// index for plate id for sample data frame
  
  // observation indexes
  DATA_IVECTOR(year_idx);  // integer index for years
  DATA_IVECTOR(depth_idx); // integer index for discrete depths.
  // prediction indexes
  DATA_IVECTOR(year_pred_idx);  // integer index for years
  DATA_IVECTOR(depth_pred_idx); // integer index for discrete depths.
    
  // Offsets
  DATA_VECTOR(ln_dilution_offset);
  DATA_VECTOR(ln_vol_offset);
  DATA_VECTOR(ln_expand_offset);
  
  // Wash effects
  DATA_VECTOR(wash_cov);
  
  // Fixed Effects Matrix
  DATA_MATRIX(Xf);   // model matrix for fixed effects
  // int n_i = y_i.rows();   // number of observations
  
  //Design matrix for Bottle random effects and design matrices
  DATA_VECTOR(Xs_ln_tau);
  DATA_VECTOR(Zs_ln_tau);
  DATA_VECTOR(Xs_ln_tau_sm);
  DATA_VECTOR(Zs_ln_tau_sm);
  DATA_INTEGER(tau_has_smooths);  // whether or not smooths are included
  
  DATA_MATRIX(X_bottle) ;
  DATA_IVECTOR(bottle_RE_idx) ; 
  //DATA_IVECTOR(drop_bottle_RE) ;
  
  // DATA_IVECTOR(b_tau_smooth_start); /// allows for different basis function numbers among smoothes.
  // DATA_INTEGER(n_smooth) ;
  
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
  DATA_STRUCT(spde, spde_t);
  DATA_STRUCT(spde_aniso, spde_aniso_t);
  //PARAMETER_ARRAY(ln_H_input);
  DATA_INTEGER(anisotropy);
  
  DATA_INTEGER(n_i);  /// number of total observations
  DATA_INTEGER(n_stand);  /// number of total observations in the standards data frame
  DATA_INTEGER(n_f);  /// number of factors per year for weight matrix L. (currently)
  DATA_INTEGER(n_d);  /// number of distinct depths in model
  DATA_INTEGER(n_y);  /// number of years in model
  DATA_INTEGER(n_s);  /// number of knot locations (constant across )
  DATA_INTEGER(n_st);  /// number of unique station locations (among all years / surveys.)
  DATA_INTEGER(n_fy);  /// number of factor-year combinations
  DATA_INTEGER(n_bottle);  /// number of factor-year combinations
  DATA_INTEGER(n_pred_locs); // number of spatial points predicted from A_pred.
  DATA_INTEGER(n_pred);  /// number of total prediction points (space x depth x years)
  DATA_INTEGER(n_d_pred);  /// number of factor-year combinations
  
  // Weight matrix releated smooth design matrices
  DATA_STRUCT(Z_L, LOM); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_STRUCT(X_L, LOM); // smoother linear effect matrix
  
  // PREDICTION MATRICES 
    // fixed matrices
  DATA_MATRIX(Xf_new); // fixed linear effect matrix for unique combinations of fixed effects.
  DATA_MATRIX(Xf_pred); // fixed linear effect matrix for new predictions of fixed effects.
    //smoothes
  DATA_STRUCT(Zs_pred, LOM); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs_pred); // smoother linear effect matrix
    // smoothes for weight matrices
  DATA_STRUCT(Z_L_pred, LOM); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_STRUCT(X_L_pred, LOM); // smoother linear effect matrix
  
  // factor indexes for weight matrix 
  DATA_IVECTOR(F_L_idx);  // integer index for years (of length n_fy)
  DATA_IVECTOR(Y_L_idx); // integer index for years (of length n_fy)
  DATA_IVECTOR(FY_start_idx); // integer index for years (of length n_fy)
  
  DATA_SPARSE_MATRIX(A_st); // INLA 'A' projection matrix for unique stations
  DATA_IVECTOR(A_ID_idx); // Vector of station ID to match up observation and A_st output

  DATA_SPARSE_MATRIX(A_pred); // INLA 'A' projection matrix for all 5km grid locations
  DATA_IVECTOR(A_ID_pred_idx); // Vector of station ID to match up observation and A_st output
  
  
    
  // PRIORS
  // Standards
  // DATA_VECTOR(std_mu_prior); 
  // DATA_VECTOR(std_sigma_prior); 
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
  
  // Hyper parameters for standards to Ct relationship.
  PARAMETER_VECTOR(std_mu); //[4]
  PARAMETER_VECTOR(ln_std_sigma) ; // Hyper prior SD for phi_0, phi_1, beta_0, beta_1.
  
  // standard parameters
  PARAMETER_VECTOR(std_phi_0) ;
  PARAMETER_VECTOR(std_phi_1) ;
  PARAMETER_VECTOR(std_beta_0) ; 
  PARAMETER_VECTOR(std_beta_1) ;
  
  PARAMETER_VECTOR(ln_std_tau); // SD of observations on the positive model component  (log-sd in Ct space)
  PARAMETER(ln_sample_tau); // extra variability for field samples.
  
  // Fixed effects
  PARAMETER_VECTOR(betaf); // fixed coefficients.
  PARAMETER(gamma_wash) ;
  
  // Bottle random effect
  PARAMETER_VECTOR(gamma_bottle);
  PARAMETER_VECTOR(ln_tau_bottle); // currently only allows a linear or 3 param smooth in log-space.
  //PARAMETER(ln_tau_bottle_test);
  
  
  // Smooth Effects
  PARAMETER_VECTOR(bs); // smoother linear effects
  PARAMETER_MATRIX(bs_L); // smoother linear effects for weight matrix (year x factor)
  
  // PARAMETER_VECTOR(ln_tau_O);    // spatial process
  // PARAMETER_ARRAY(ln_tau_Z);    // optional spatially varying covariate process
  // PARAMETER_VECTOR(ln_tau_E);    // spatio-temporal process
  // PARAMETER_ARRAY(ln_kappa);    // Matern parameter
  
  // Random effects
  PARAMETER_VECTOR(b_smooth);  // P-spline smooth parameters
  PARAMETER_VECTOR(ln_smooth_sigma);  // variances of spline REs if included

  PARAMETER_ARRAY(b_L_smooth);  // P-spline smooth parameters for L matrix (dim(year,factor,smooth knots))
  PARAMETER(ln_L_smooth_sigma);  // variances of spline L REs if included (dim(year,factor))
  // 
  PARAMETER(ln_kappa);    // Matern parameter determining range.
  // 
  // //Random effects
  PARAMETER_MATRIX(omega_s);    // spatial effects; (number of knots x (years by number of factors)

    
    /////// END PARAMETERS /////////
  
  // Joint Negative log-likelihood.
  Type jnll = 0. ;

  // Spline penalties
  // p-splines/smoothers

   matrix<Type>  eta_smooth_i(n_i,n_smooth);
   eta_smooth_i.setZero();
   matrix<Type>  eta_smooth_pred(n_pred,n_smooth);
   eta_smooth_pred.setZero();
   vector<Type>  eta_smooth_all;
   eta_smooth_all.setZero() ;
   vector<Type>  eta_smooth_all_pred;
   eta_smooth_all_pred.setZero() ;
   
if(has_smooths){  
  for(int k = 0; k < n_smooth; k++){
    vector<Type> beta_s(Zs(k).cols());
    beta_s.setZero();
         for (int j = 0; j < beta_s.rows(); j++) {
           beta_s(j) = b_smooth(b_smooth_start(k) + j);
           PARALLEL_REGION jnll -= dnorm(beta_s(j), Type(0), exp(ln_smooth_sigma(k)), true); // MADE THIS SIMPLE BECAUSE COVARIATE SMOOTHES ARE ALWAYS LOW DIMENSIONAL.
         }
    eta_smooth_i.col(k) = Zs(k) * beta_s ;
    eta_smooth_pred.col(k) = Zs_pred(k) * beta_s ;
  }

  eta_smooth_all = eta_smooth_i.rowwise().sum();
  eta_smooth_all_pred = eta_smooth_pred.rowwise().sum();

  REPORT(bs);
  REPORT(b_smooth);     // smooth coefficients for penalized splines
  REPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
  REPORT(eta_smooth_i) ;
  REPORT(eta_smooth_all) ;
  
  ADREPORT(bs);
  ADREPORT(b_smooth);     // smooth coefficients for penalized splines
  ADREPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
  //REPORT(eta_smooth_i) ;
  //REPORT(eta_smooth_all) ;
}
REPORT(betaf) ;
ADREPORT(betaf) ;

// ------------------ Geospatial ---------------------------------------------
// 
// Matern Range calculations
Type kappa;
Type range;
  //for (int f = 0; f < n_f; f++) {
    kappa = exp(ln_kappa) ; 
    range = sqrt(Type(8.)) / kappa;
    
  //}
  REPORT(range);
  REPORT(kappa);
  REPORT(ln_kappa);
  
  ADREPORT(range);
  ADREPORT(kappa);
  ADREPORT(ln_kappa);

// calculate tau for each factor from kappa assuming sigma == 1
Type tau_O;
Type siggy;
for(int f=0; f < n_f; f++){
  tau_O = calc_tau(ln_kappa) ;
  siggy = calc_sigma(log(tau_O),ln_kappa);
}
REPORT(tau_O);
ADREPORT(tau_O);
REPORT(siggy);

// Make spde model.
Eigen::SparseMatrix<Type> Q_f; // Precision matrix for factors (one for each factor)

  if(anisotropy==1) {
     Q_f = R_inla::Q_spde(spde, exp(ln_kappa));
  }
  if(anisotropy==0) {
    // matrix<Type> H = MakeH(vector<Type>(ln_H_input.col(0)));
    // Q_f = R_inla::Q_spde(spde_aniso, exp(ln_kappa), H);
    // REPORT(H);
  }
  for(int fy = 0; fy < n_fy; fy++){ // loop over
    
    // CANNOT MAKE THE KAPPA LOOP OVER factors.  No idea why it doesn't work!!!!!
    // Q_f = R_inla::Q_spde(spde, exp(ln_kappa(fy)));
    // Q_f = R_inla::Q_spde(spde, exp(ln_kappa(F_L_idx(fy))));
    PARALLEL_REGION jnll += SCALE(GMRF(Q_f), 1. / tau_O)(omega_s.col(fy));
    //PARALLEL_REGION jnll += GMRF(Q_f,true)(omega_s.col(fy));
  }
  
/////////////////////////////////////

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
  vector<matrix<Type>> L_mat(n_y); // vector of weight matrices.  One for each year.
  vector<matrix<Type>> L_mat_pred(n_y); // vector of weight matrices.  One for each year.
  
  // make L matrix from smooth design matrix
  //  note that each factor has one less entry than the factor before it (identifibility)
  for(int y = 0; y < n_y; y++){
    matrix<Type> tmp_mat(n_d,n_f) ;
    tmp_mat.setZero();
    matrix<Type> tmp_mat_pred(n_d_pred,n_f) ;
    tmp_mat_pred.setZero();
    //vector<Type> A(n_d);
    for(int f=0; f < n_f;f++){ // loop over factors in each year
      vector<Type> beta_L(Z_L(f).cols());
      beta_L.setZero();

      for (int j = 0; j < beta_L.size(); j++) {
        beta_L(j) = b_L_smooth(f,j) ;
        PARALLEL_REGION jnll -= dnorm(beta_L(j), Type(0), exp(ln_L_smooth_sigma), true);
      }

      tmp_mat.col(f) = matrix<Type>(vector<Type>(X_L(f) * vector<Type>(bs_L.row(f))) + vector<Type>(Z_L(f) * beta_L)) ;
      tmp_mat_pred.col(f) = matrix<Type>(vector<Type>(X_L_pred(f) * vector<Type>(bs_L.row(f))) + vector<Type>(Z_L_pred(f) * beta_L)) ;
      
      // matrix<Type> m1_block = m1.block(0,0,1,1);
      // REPORT(m1_block); // block starting (0,0) taking 1 row & 1 col
      // //
      // int END = n_d-f ;
      // matrix<Type>(THIS)(END,1);
      //
      // THIS = matrix<Type>(vector<Type>(X_L(f) * bs_L(y,f)) + Z_L(f)*beta_L) ;
      // //tmp_mat.block(f,f,END,1) = matrix<Type>(vector<Type>(X_L(f) * bs_L(y,f)) + Z_L(f)*beta_L) ;


      // for (int j = 0; j < beta_L.rows(); j++) {
      //   beta_L(j) = b_L_smooth(j);
      //   PARALLEL_REGION jnll -= dnorm(beta_L(j), Type(0), exp(ln_L_smooth_sigma(f)), true);
      // } // end j loop
      // tmp_mat.col(f) = matrix<Type>(vector<Type>(X_L(y) * bs_L) + Z_L(y) * beta_L) ;
    } // end f loop
    L_mat(y) = tmp_mat ;
    L_mat_pred(y) = tmp_mat_pred ;
  } //end y loop

   
   // // Add prior for ln_L_smooth_sigma 
   // for(int y=0; y<n_y; y++){
   //   for(int f=0; f<n_f; f++){
   //     jnll -= dnorm(bs_L(y,f), Type(0), Type(3),true);
   //   }
   // }
   

  REPORT(bs_L) ;
  REPORT(b_L_smooth) ;
  REPORT(ln_L_smooth_sigma) ;  
  ADREPORT(bs_L) ;
  ADREPORT(b_L_smooth) ;
  ADREPORT(ln_L_smooth_sigma) ; 
  REPORT(L_mat) ;
  REPORT(L_mat_pred) ;
  ADREPORT(L_mat_pred(0)) ;
  // Make a epsilon matrix for each year with n_d (# depth_cat) columns, n_s (# knots) rows.
  vector<matrix<Type>> epsilon_s(n_y);
  vector<matrix<Type>> epsilon_s_pred(n_y);
  for(int y = 0; y < n_y; y++){
       // matrix<Type> epsilon_temp(n_s,n_d);
       // epsilon_temp.setZero() ; 
       // matrix<Type> epsilon_temp_pred(n_s,n_d_pred);
       // epsilon_temp_pred.setZero() ;   
       
       // if(n_f == 1){ // if there is a single, shared field for all depths.
       //    for(int d = 0; d< n_d; d++){
       //      epsilon_temp.col(d) = omega_s.col(FY_start_idx(y));
       //    }
       //    epsilon_s(y) = epsilon_temp;
       // }else{ // if there is a more than one factor (>1 spatial field among depths.)
          
          if(n_f == 1){
            epsilon_s(y) = matrix<Type>(omega_s.col(FY_start_idx(y)) * L_mat(y).transpose()) ;
            epsilon_s_pred(y) = matrix<Type>(omega_s.col(FY_start_idx(y)) * L_mat_pred(y).transpose()) ;
          }else{
            epsilon_s(y) = omega_s.block(0,FY_start_idx(y),n_s,n_f) * L_mat(y).transpose() ;
            epsilon_s_pred(y) = omega_s.block(0,FY_start_idx(y),n_s,n_f) * L_mat_pred(y).transpose() ;
          }
  }
  
  // // ------------------ INLA projections ---------------------------------------
  //
  // Here we are projecting the spatiotemporal and spatial random effects to the
  // locations of the data using the INLA 'A' matrices.

  // Make the spatial projection from knots to station locations and predicted field locations.
  array<Type> epsilon_st(n_st,n_d,n_y);
  array<Type> epsilon_pred(n_pred_locs,n_d_pred,n_y);
  
  for(int y = 0; y < n_y; y++){
    for(int d = 0; d < n_d; d++){
       epsilon_st.col(y).col(d) = (A_st * epsilon_s(y).col(d)).array() ;
       //Prediction to 5km grid
    }
    for(int d = 0; d < n_d_pred; d++){
      //Prediction to 5km grid at all new depths.
       epsilon_pred.col(y).col(d) = (A_pred * epsilon_s_pred(y).col(d)).array() ;
    }
  }
  REPORT(epsilon_s) ; // epsilon at each knot location.
  REPORT(epsilon_st) ;  // epsilon at each station and depth.
  REPORT(epsilon_pred) ;  // epsilon at predicted spatial location and depth.
  
  vector<Type> D_i(n_i); // latent variable at the year-station-depth level for each observation.
  vector<Type> D_pred(n_pred); // latent variable at the year-station-depth level for each observation.
  vector<Type> E_i(n_i); // latent variable at the year-station-depth level for each observation.
  
  
  vector<Type> F_new;
  F_new = Xf_new * betaf ;
  //S_new = Xs*bs + eta_smooth_all ;
  
  D_i =  Xf*betaf + Xs*bs + eta_smooth_all  ;  // all fixed and smooth effects.
  D_pred =  Xf_pred*betaf + Xs_pred*bs + eta_smooth_all_pred  ;  // all fixed and smooth effects for prediction set
   
   for(int i=0; i < n_i; i++){
     D_i(i) = D_i(i) + epsilon_st(A_ID_idx(i),depth_idx(i),year_idx(i)); // Add spatial effects.
   }
   for(int i=0; i < n_pred; i++){
     D_pred(i) = D_pred(i) + epsilon_pred(A_ID_pred_idx(i),depth_pred_idx(i),year_pred_idx(i)); // Add spatial effects.
   }
   // So D_i is the average station-level log DNA in a common currency (2.5L sampled, undiluted, washed properly, etc.)
   
  /// MAKE OVERALL INDEX OF ABUNDANCE BASED ON THE PREDICTION SET.
    // Type INDEX ;
    // //<Type> INDEX_trim ;
    // 
    // INDEX = sum(exp(D_i)) ;
    // 
    // REPORT(INDEX) ;
    // ADREPORT(INDEX) ;
   
   ////////////////////////////////////////////////////////
   // bottle random effects
   
   // First make the smooth to define SD as a function of depth_category.
   vector<Type>tau_bottle(n_bottle) ;
   tau_bottle.setZero();
   vector<Type>tau_bottle_sm ;
   tau_bottle_sm.setZero();
   // vector<Type>zero_bottle(n_bottle) ;
   // zero_bottle.setZero() ;
   
   if(tau_has_smooths==0){  // 3 knot smooth in log space.
      tau_bottle = vector<Type>(exp(ln_tau_bottle(0) +
                         vector<Type>(Xs_ln_tau*ln_tau_bottle(1)) +
                         vector<Type>(Zs_ln_tau*ln_tau_bottle(2)))) ;
      tau_bottle_sm = vector<Type>(exp(ln_tau_bottle(0) +
                         vector<Type>(Xs_ln_tau_sm*ln_tau_bottle(1)) +
                         vector<Type>(Zs_ln_tau_sm*ln_tau_bottle(2)))) ;
   }
   if(tau_has_smooths==1){ // linear in log-space
     tau_bottle = vector<Type>(exp(ln_tau_bottle(0) +
                      vector<Type>(Xs_ln_tau*ln_tau_bottle(1)))) ;
     tau_bottle_sm = vector<Type>(exp(ln_tau_bottle(0) +
                      vector<Type>(Xs_ln_tau_sm*ln_tau_bottle(1)))) ;
   }
   // // // Random effect for bottle deviation. 
   // This is not in vector form because I want to drop the instances where there is a single bottle 
   // from a given station-depth combination (the if statement.... likely more efficient to vectorize it.)
  
   for(int b = 0; b < n_bottle; b++){
      jnll -= dnorm(gamma_bottle(bottle_RE_idx(b)) , Type(0.), tau_bottle(bottle_RE_idx(b)),true);
   }
   // 
   // VECTOR VERSION WITHOUT THE IF STATEMENT. 
   //  jnll -= sum(dnorm(gamma_bottle , zero_bottle  , tau_bottle,true)); 
    
    jnll -= dnorm( ln_tau_bottle(0), Type(-2), Type(0.5),true);

   // Constrain the bottle SD to be small-ish... this is probably statistically illegal....
    // for(int d = 0; tau_bottle_sm.size();d++){
    //   jnll -= dnorm( Type(tau_bottle_sm(d)), Type(0), Type(0.5),true);
    // }

   REPORT(ln_tau_bottle) ;
   
   REPORT(tau_bottle_sm);
   ADREPORT(tau_bottle_sm);
   REPORT(tau_bottle);     // smooth coefficients for penalized splines
   //ADREPORT(tau_bottle);     // smooth coefficients for penalized splines
   REPORT(gamma_bottle); // standard deviations of smooth random effects, in log-space
   
  ///////////////////////////////////////////////////////////// 
   // This adds offsets for wash error, dilution offsets, volume offsets, bottle random effect....
   E_i = D_i +  ln_dilution_offset + ln_vol_offset + ln_expand_offset + wash_cov * gamma_wash + X_bottle * gamma_bottle ;

   
 ADREPORT(F_new);   
 REPORT(D_i);
 REPORT(D_pred);
 //ADREPORT(D_i);
 REPORT(D_pred);
 REPORT(exp(D_i));
 REPORT(E_i)

  ///////////////////// PRIORS FOR ALL PARAMETERS
  // Add prior for ln_L_smooth_sigma 
  //for(int y=0; y<n_y; y++){
    // for(int f=0; f < n_f; f++){
    //   jnll -= dnorm(exp(ln_L_smooth_sigma(f,0)), Type(0), Type(1),true);
    // }
  //}
  
  //Kappa prior
  jnll -= dnorm(ln_kappa, Type(-4), Type(2),true);
  
  ///////////////////////////////////////////////////////////////////  
  // contribution of the standards to the likeihoods
  // section modified slightly from Eric's sdmTMB code.
  
    // Presence-Absence component of model.
    vector<Type> theta_stand(n_stand);
    vector<Type> std_sigma = exp(ln_std_sigma);
    vector<Type> kappa_stand(n_stand);    
    
    // hierarchical contribution for logit scale regression coefficients.
    for(int i = 0; i < std_phi_0.size(); i++){
      jnll -= dnorm(std_phi_0(i), std_mu(0), std_sigma(0), true);//phi_0 ~ normal(2, 2)
      jnll -= dnorm(std_phi_1(i), std_mu(1), std_sigma(1), true);//phi_1 ~ normal(4, 2)
    }
    // hierarchical contribution for log-scale regression coefficients.
    for(int i = 0; i < std_beta_0.size(); i++){
      jnll -= dnorm(std_beta_0(i), std_mu(2), std_sigma(2), true);//beta_0 ~ normal(std_mu(2))
      jnll -= dnorm(std_beta_1(i), std_mu(3), std_sigma(3), true);//beta_1 ~ normal()
    }
    
    // Likelihoods:
    for(int i = 0; i < n_stand; i++){ // likelihood
      theta_stand(i) = std_phi_0(plate_stand_idx(i)) + std_phi_1(plate_stand_idx(i)) * (ln_copies(i));// - stand_offset);
      jnll -= dbinom_robust(Z_stand(i), Type(1), theta_stand(i), true); // likelihood
      if(Z_stand_int(i)==1){
        kappa_stand(i) = std_beta_0(plate_stand_idx(i)) + std_beta_1(plate_stand_idx(i)) * (ln_copies(i));// - stand_offset);
        jnll -= dnorm(Y_stand(i), kappa_stand(i), exp(ln_std_tau(0) + ln_std_tau(1)*ln_copies(i)), true); // likelihood  
      }
    }

    // Repeat for likelihood contributions of unknown Samples 
    //Presence-Absence component of model.
     vector<Type> theta_samp(n_i);
     vector<Type> kappa_samp(n_i);

     for(int i = 0; i < n_i; i++){ // likelihood
        theta_samp(i) = std_phi_0(plate_samp_idx(i)) + std_phi_1(plate_samp_idx(i)) * (E_i(i));// - stand_offset);
        jnll -= dbinom_robust(Z_samp(i), Type(1), theta_samp(i), true); // likelihood
        if(Z_samp_int(i)==1){
          kappa_samp(i) = std_beta_0(plate_samp_idx(i)) + std_beta_1(plate_samp_idx(i)) * (E_i(i));// - stand_offset);
          jnll -= dnorm(Y_samp(i), kappa_samp(i), exp(ln_sample_tau + ln_std_tau(0) + ln_std_tau(1)*E_i(i)), true); // likelihood  
        }
     }
    
    
    // PRIOR ON REFERENCE YEAR INTERCEPT (2019)
    for(int y=0; y < n_y; y++){
      if(y==0){
        jnll -= dnorm(betaf(y),Type(1), Type(2),true);
      }else{jnll -= dnorm(betaf(y),Type(0), Type(2),true);}
    }
    
    // Extra variablity in field samples.  Prior to make this small (normal on the log scale)
          jnll -= dnorm(ln_sample_tau,Type(-2),Type(0.5),true) ;
     

     REPORT(std_mu);
     REPORT(std_sigma);
     REPORT(std_phi_0);
     REPORT(std_phi_1);
     REPORT(std_beta_0);
     REPORT(std_beta_1);
     
     ADREPORT(std_mu);
     ADREPORT(std_sigma);
     ADREPORT(std_phi_0);
     ADREPORT(std_phi_1);
     ADREPORT(std_beta_0);
     ADREPORT(std_beta_1);
     
     REPORT(ln_std_tau);
     ADREPORT(ln_std_tau);
     REPORT(exp(ln_std_tau));
     REPORT(ln_sample_tau);     
     
     REPORT(theta_stand);
     REPORT(kappa_stand);

     REPORT(theta_samp);
     REPORT(kappa_samp);
  
     REPORT(jnll);   
     
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
