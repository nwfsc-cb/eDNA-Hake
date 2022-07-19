// This is a model file for calculating the qPCR results for eDNA sampling in 2019.

data { /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Number of observations in various categories
    int N_sample; // Number of unique bottles used combinations observed.
    int N_control_sample; // Number of samples used for negative controls.
    int N_station_depth; // Number of station-depth combinations observed.
    int N_pcr ;    // Number of PCR plates
    int N_depth ;    // Number of depth categories
    
    int N_stand_bin ;   // Number of observations for binomial part of the standards model
    int N_stand_pos ; // Number of observations for count part of the standards model
    int N_obs_bin   ;  // Number of observations for binomial part of the sample model
    int N_obs_pos ;  // Number of observations for count part of the sample model
    int N_control_bin   ;  // Number of observations for binomial part of the sample model
    int N_control_pos ;  // Number of observations for count part of the sample model

    // Observations
    int bin_stand[N_stand_bin]     ;
    vector[N_stand_pos] pos_stand ;
    int bin_obs[N_obs_bin]      ; 
    vector[N_obs_pos] pos_obs   ;
    int bin_control[N_control_bin]      ; 
    vector[N_control_pos] pos_control   ;
    
    // Covariates (standards) (log counts)
    vector[N_stand_bin] D_bin_stand     ;
    vector[N_stand_pos] D_pos_stand ;
    
    // Covariates (samples) (log10(volume filtered))
    vector[N_sample] log_vol_obs ;
    
    vector[N_control_bin] bin_log_vol_control ;
    vector[N_control_pos] pos_log_vol_control ;

    vector[N_obs_bin] bin_log_dilution_obs ;
    vector[N_obs_pos] pos_log_dilution_obs ;
    vector[N_control_bin] bin_log_dilution_control ;
    vector[N_control_pos] pos_log_dilution_control ;
    
    // Standard Indices
    int pcr_stand_bin_idx[N_stand_bin] ;
    int pcr_stand_pos_idx[N_stand_pos] ;
    int pcr_obs_bin_idx[N_obs_bin] ;
    int pcr_obs_pos_idx[N_obs_pos] ;
    int pcr_control_bin_idx[N_control_bin] ;
    int pcr_control_pos_idx[N_control_pos] ;

    // LOO Index
    int loo_pos_idx[N_obs_pos] ;

    // Station-depth and sample indices
    int sample_idx[N_sample]    ;
    int samp_station_depth_idx[N_sample] ;
    real wash_idx[N_sample] ;
    real singleton_idx[N_sample] ;
    real Ct_bin_idx[N_sample] ;
    int sample_control_idx[N_control_sample]    ;
    int station_depth_idx[N_station_depth] ;
    int depth_idx[N_sample] ;

    // Sample related indices
    int sample_bin_idx[N_obs_bin]      ;
    int sample_pos_idx[N_obs_pos]  ;
    int sample_control_bin_idx[N_control_bin]      ;
    int sample_control_pos_idx[N_control_pos]  ;
    real wash_control_bin_idx[N_control_bin]      ;
    real wash_control_pos_idx[N_control_pos]  ;
    
    //int station_depth_bin_idx[N_obs_bin]      ;
    //int station_depth_pos_idx[N_obs_pos]  ;

    int single_n_control_pos; // This is an indicator for dealing with single values where program is looking for a vector

    //Offset
    real OFFSET ;
    real wash_offset_prior[2];
    
  ///////////////////////////////////////////////////////////////  
  ///// SMOOTH COMPONTENTS, extracted from BRMS package
  ///////////////////////////////////////////////////////////////
  // int<lower=1> N_station_depth;  // total number of STATION.DEPTH COMBINATIONS>
  int<lower=1> K;  // number of population-level effects
  matrix[N_station_depth, K] X;  // population-level design matrix
  // data for splines
  int Ks;  // number of linear effects
  matrix[N_station_depth, Ks] Xs;  // design matrix for the linear effects
  
  // data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)0
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N_station_depth, knots_1[1]] Zs_1_1;
  matrix[N_station_depth, knots_1[2]] Zs_1_2;
  matrix[N_station_depth, knots_1[3]] Zs_1_3;
  // data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)50
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N_station_depth, knots_2[1]] Zs_2_1;
  matrix[N_station_depth, knots_2[2]] Zs_2_2;
  matrix[N_station_depth, knots_2[3]] Zs_2_3;
  // data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)100
  int nb_3;  // number of bases
  int knots_3[nb_3];  // number of knots
  // basis function matrices
  matrix[N_station_depth, knots_3[1]] Zs_3_1;
  matrix[N_station_depth, knots_3[2]] Zs_3_2;
  matrix[N_station_depth, knots_3[3]] Zs_3_3;
  // data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)150
  int nb_4;  // number of bases
  int knots_4[nb_4];  // number of knots
  // basis function matrices
  matrix[N_station_depth, knots_4[1]] Zs_4_1;
  matrix[N_station_depth, knots_4[2]] Zs_4_2;
  matrix[N_station_depth, knots_4[3]] Zs_4_3;
  // data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)300
  int nb_5;  // number of bases
  int knots_5[nb_5];  // number of knots
  // basis function matrices
  matrix[N_station_depth, knots_5[1]] Zs_5_1;
  matrix[N_station_depth, knots_5[2]] Zs_5_2;
  matrix[N_station_depth, knots_5[3]] Zs_5_3;
  // data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)500
  int nb_6;  // number of bases
  int knots_6[nb_6];  // number of knots
  // basis function matrices
  matrix[N_station_depth, knots_6[1]] Zs_6_1;
  matrix[N_station_depth, knots_6[2]] Zs_6_2;
  matrix[N_station_depth, knots_6[3]] Zs_6_3;
  
  
  // data for spline s(bottom.depth,k=N.knots.bd)
  int nb_7;  // number of bases
  int knots_7[nb_7];  // number of knots
  // basis function matrices
  matrix[N_station_depth, knots_7[1]] Zs_7_1;
  
}
transformed data{
  /// SMOOTHES TRANSFORMED DATA
  int Kc = K - 1;
  matrix[N_station_depth, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  
}
parameters { /////////////////////////////////////////////////////////////////////////////////////////////

    // Contamination
       real mu_contam;             // log mean contamination derived from the control samples.
       real<lower=0> sigma_contam; // log sd contamination

      real wash_offset ;

 // Variance parameters for observation and random effects
       real<lower=0> sigma_stand_int ;
       real<lower=0> tau_sample[N_depth] ;
       real<lower=0> sigma_pcr ; 
      // real<lower=2> nu ;
       
    // Standards regression and logit coeffs
      real beta_0[N_pcr] ;
      real beta_1[N_pcr] ;
      
      real phi_0[N_pcr] ;
      real phi_1[N_pcr] ;

     //real sigma_stand_slope ;
     // real sigma_stand_slope2 ;
      
      // real sigma_stand_int_bar ;
      // real<lower=0> sigma_stand_int_sd ;
      // real sigma_stand_slope_bar ;
      // real<lower=0> sigma_stand_slope_sd ;
      // 
      //vector<lower=0>[N_pcr] sigma_stand_int ;

    // // Effects of station-depths and samples
       //real D[N_station_depth] ;
       //real D_error[N_sample] ;
       real D_control[N_control_sample] ;
       real delta[N_sample] ;
       
   ////// SMOOTHES ////////////////////////////
    vector[K] b;  // population-level effects
    //real Intercept;  // temporary intercept for centered predictors
    vector[Ks] bs;  // spline coefficients
    // parameters for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)0
    // standarized spline coefficients
    vector[knots_1[1]] zs_1_1;
    vector[knots_1[2]] zs_1_2;
    vector[knots_1[3]] zs_1_3;
    real<lower=0> sds_1_1;  // standard deviations of spline coefficients
    real<lower=0> sds_1_2;  // standard deviations of spline coefficients
    real<lower=0> sds_1_3;  // standard deviations of spline coefficients
  
    // parameters for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)50
    // standarized spline coefficients
    vector[knots_2[1]] zs_2_1;
    vector[knots_2[2]] zs_2_2;
    vector[knots_2[3]] zs_2_3;
    real<lower=0> sds_2_1;  // standard deviations of spline coefficients
    real<lower=0> sds_2_2;  // standard deviations of spline coefficients
    real<lower=0> sds_2_3;  // standard deviations of spline coefficients
    // parameters for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)100
    // standarized spline coefficients
    vector[knots_3[1]] zs_3_1;
    vector[knots_3[2]] zs_3_2;
    vector[knots_3[3]] zs_3_3;
    real<lower=0> sds_3_1;  // standard deviations of spline coefficients
    real<lower=0> sds_3_2;  // standard deviations of spline coefficients
    real<lower=0> sds_3_3;  // standard deviations of spline coefficients
  
    // parameters for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)150
    // standarized spline coefficients
    vector[knots_4[1]] zs_4_1;
    vector[knots_4[2]] zs_4_2;
    vector[knots_4[3]] zs_4_3;
    real<lower=0> sds_4_1;  // standard deviations of spline coefficients
    real<lower=0> sds_4_2;  // standard deviations of spline coefficients
    real<lower=0> sds_4_3;  // standard deviations of spline coefficients
    // parameters for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)300
    // standarized spline coefficients
    vector[knots_5[1]] zs_5_1;
    vector[knots_5[2]] zs_5_2;
    vector[knots_5[3]] zs_5_3;
    real<lower=0> sds_5_1;  // standard deviations of spline coefficients
    real<lower=0> sds_5_2;  // standard deviations of spline coefficients
    real<lower=0> sds_5_3;  // standard deviations of spline coefficients
    // parameters for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)500
    // standarized spline coefficients
    vector[knots_6[1]] zs_6_1;
    vector[knots_6[2]] zs_6_2;
    vector[knots_6[3]] zs_6_3;
    real<lower=0> sds_6_1;  // standard deviations of spline coefficients
    real<lower=0> sds_6_2;  // standard deviations of spline coefficients
    real<lower=0> sds_6_3;  // standard deviations of spline coefficients
  
    // parameters for spline s(bottom.depth,k=N.knots.bd)
    // standarized spline coefficients
    vector[knots_7[1]] zs_7_1;
    real<lower=0> sds_7_1;  // standard deviations of spline coefficients
}
transformed parameters { ////////////////////////////////////////////////////////////////////////////////
      // Variance Parameters
      real sigma_all_stand ;
      real sigma_all_samp ;
      
      // Latent variable for contaminated log-concentration of DNA that is observed by PCR.
      vector[N_station_depth] D;
      real D_contam[N_sample] ;
      real D_delta[N_sample] ;
      //vector[N_station_depth] mu_smooth;

  // actual spline coefficients
    vector[knots_1[1]] s_1_1;
    vector[knots_1[2]] s_1_2;
    vector[knots_1[3]] s_1_3;
    vector[knots_2[1]] s_2_1;
    vector[knots_2[2]] s_2_2;
    vector[knots_2[3]] s_2_3;
    vector[knots_3[1]] s_3_1;
    vector[knots_3[2]] s_3_2;
    vector[knots_3[3]] s_3_3;
    vector[knots_4[1]] s_4_1;
    vector[knots_4[2]] s_4_2;
    vector[knots_4[3]] s_4_3;
    vector[knots_5[1]] s_5_1;
    vector[knots_5[2]] s_5_2;
    vector[knots_5[3]] s_5_3;
    vector[knots_6[1]] s_6_1;
    vector[knots_6[2]] s_6_2;
    vector[knots_6[3]] s_6_3;
    vector[knots_7[1]] s_7_1;
  
  // compute actual spline coefficients
    s_1_1 = sds_1_1 * zs_1_1;
    s_1_2 = sds_1_2 * zs_1_2;
    s_1_3 = sds_1_3 * zs_1_3;
    s_2_1 = sds_2_1 * zs_2_1;
    s_2_2 = sds_2_2 * zs_2_2;
    s_2_3 = sds_2_3 * zs_2_3;
    s_3_1 = sds_3_1 * zs_3_1;
    s_3_2 = sds_3_2 * zs_3_2;
    s_3_3 = sds_3_3 * zs_3_3;
    s_4_1 = sds_4_1 * zs_4_1;
    s_4_2 = sds_4_2 * zs_4_2;
    s_4_3 = sds_4_3 * zs_4_3;
    s_5_1 = sds_5_1 * zs_5_1;
    s_5_2 = sds_5_2 * zs_5_2;
    s_5_3 = sds_5_3 * zs_5_3;
    s_6_1 = sds_6_1 * zs_6_1;
    s_6_2 = sds_6_2 * zs_6_2;
    s_6_3 = sds_6_3 * zs_6_3;
    s_7_1 = sds_7_1 * zs_7_1;
  
    D = // Intercept + rep_vector(0.0, N_station_depth) + 
                    X * b +  // Factor level effects
                    Xs * bs + // linear effects of smoothes
                    Zs_1_1 * s_1_1 + Zs_1_2 * s_1_2 + Zs_1_3 * s_1_3 + 
                    Zs_2_1 * s_2_1 + Zs_2_2 * s_2_2 + Zs_2_3 * s_2_3 + 
                    Zs_3_1 * s_3_1 + Zs_3_2 * s_3_2 + Zs_3_3 * s_3_3 + 
                    Zs_4_1 * s_4_1 + Zs_4_2 * s_4_2 + Zs_4_3 * s_4_3 + 
                    Zs_5_1 * s_5_1 + Zs_5_2 * s_5_2 + Zs_5_3 * s_5_3 + 
                    Zs_6_1 * s_6_1 + Zs_6_2 * s_6_2 + Zs_6_3 * s_6_3 + 
                    Zs_7_1 * s_7_1;
  
    // ADD IN A LATENT VARIABLE TO ACCOUNT FOR CONTAMINATION OF THE FIELD SAMPLES.
    for(i in 1:N_sample){
      D_delta[i] = D[samp_station_depth_idx[i]] + // This is the log-DNA concentration in each Niskin bottle.
                        delta[i];// * singleton_idx[i] * Ct_bin_idx[i]; //
      D_contam[i] = D_delta[i]  +
                        log_vol_obs[i] +
                        wash_offset * wash_idx[i]; 
            // + exp(D_error[i])) ;      //D_contam[i] = log(exp(D[i]) + exp(D_error[i])) ;
    }
    
    // Derive Variances for the different components of the model.
    //for( i in 1:N_count_stand){
      sigma_all_stand = sigma_stand_int;
                                //+ sigma_stand_slope * (D_count_stand[i]- OFFSET)),-2)   ;
                                //+ sigma_stand_slope2 * D_count_stand[i]^2),-2) ;

    //}  
    //for( i in 1:N_count_samp){
      sigma_all_samp = sigma_pcr;
                        // pow(sigma_stand_int^2 +
                        //     //sigma_stand_slope * (D[bottle_count_idx[i]]- OFFSET)) +
                        //     //+ sigma_stand_slope2 * D[bottle_count_idx[i]]^2) + 
                        //     sigma_pcr^2,-2) ;
    //} 
    
}
model {////////////////////////////////////////////////////////////////////////////////////////////////////
    { //LOCAL VARIABLES DECLARATION START
      vector[N_stand_bin] theta_stand ;
      vector[N_obs_bin] theta_obs ;
      vector[N_control_bin] theta_control ;

      vector[N_stand_pos] kappa_stand ;
      vector[N_obs_pos] kappa_obs ;
      vector[N_control_pos] kappa_control ;
      
    // Presence-Absence component of model.
    for(i in 1:N_stand_bin){
       theta_stand[i] = phi_0[pcr_stand_bin_idx[i]] + phi_1[pcr_stand_bin_idx[i]] *  (D_bin_stand[i] - OFFSET) ;
    }
    
    for(i in 1:N_obs_bin){
       theta_obs[i]  = phi_0[pcr_obs_bin_idx[i]] + 
                      phi_1[pcr_obs_bin_idx[i]] * (D_contam[sample_bin_idx[i]] + 
                                                      bin_log_dilution_obs[i] - 
                                                      OFFSET);
    }

    for(i in 1:N_control_bin){
       theta_control[i]  = phi_0[pcr_control_bin_idx[i]] + 
                      phi_1[pcr_control_bin_idx[i]] * (D_control[sample_control_bin_idx[i]] +
                                                        wash_offset * wash_control_bin_idx[i] +
                                                        bin_log_dilution_control[i] +
                                                        bin_log_vol_control[i] - OFFSET) ;
    }


    // Positive Comonent of the model
    for(i in 1:N_stand_pos){
       kappa_stand[i] = beta_0[pcr_stand_pos_idx[i]] + beta_1[pcr_stand_pos_idx[i]] * (D_pos_stand[i] - OFFSET) ;
    }
    
    for(i in 1:N_obs_pos){
      kappa_obs[i]  =  beta_0[pcr_obs_pos_idx[i]] + 
                        beta_1[pcr_obs_pos_idx[i]] * (D_contam[sample_pos_idx[i]] + 
                                                      pos_log_dilution_obs[i] - 
                                                      OFFSET) ;
    }
    
    for(i in 1:(N_control_pos-single_n_control_pos)){
      kappa_control[i]  = beta_0[pcr_control_pos_idx[i]] + 
                        beta_1[pcr_control_pos_idx[i]] * (D_control[sample_control_pos_idx[i]] + 
                                                          wash_offset * wash_control_pos_idx[i] +
                                                          pos_log_dilution_control[i] +
                                                          pos_log_vol_control[i] - OFFSET) ;
    }

    // Likelihood components
    bin_stand     ~ bernoulli( inv_logit(theta_stand) ) ;
    bin_obs       ~ bernoulli( inv_logit(theta_obs) ) ;
    bin_control   ~ bernoulli( inv_logit(theta_control) ) ;
    
    pos_stand   ~ normal(kappa_stand, sigma_all_stand) ;
    pos_obs     ~ student_t(3, kappa_obs, sigma_all_samp) ;
    //pos_obs     ~ normal( kappa_obs, sigma_all_samp) ;
    if(single_n_control_pos==1){
      pos_control[1] ~ student_t(3, kappa_control[1], sigma_all_samp) ;
      //pos_control[1] ~ normal( kappa_control[1], sigma_all_samp) ;
    }else{
      pos_control ~ student_t(3, kappa_control, sigma_all_samp) ;
      //pos_control ~ normal( kappa_control, sigma_all_samp) ;
    }
    
    } //LOCAL VARIABLES DECLARATION END 


    // Parameters for the gBLOCKS:
     // Random effects
      //D ~ normal(0,3);
      D_control ~ normal(mu_contam,sigma_contam); // Estimating contamination.
      //D_error   ~ normal(mu_contam,sigma_contam); // Unobserved contamination... D_error is an estimated latent variable.
      
      // Priors
      //sigma_stand_int ~ gamma(1,1) ;
      target += normal_lpdf(sigma_stand_int | 0, 2) - 1 * normal_lccdf(0 | 0, 2);
      //sigma_pcr ~ gamma(1,1) ;
      target += normal_lpdf(sigma_pcr | 0, 2) - 1 * normal_lccdf(0 | 0, 2);
      //target += normal_lpdf(nu | 2, 10) - 1 * normal_lccdf(2 | 2, 10);

      beta_0 ~ normal(38,4) ;
      beta_1 ~ normal(-3.32,0.1) ;

      phi_0 ~ normal(2, 2) ;
      phi_1 ~ normal(4, 2) ;

      for(i in 1:N_sample){
        delta[i] ~ normal(0,tau_sample[depth_idx[i]]) ;
      }
      target += normal_lpdf(tau_sample | 0, 1) - 1 * normal_lccdf(0 | 0, 1);

      mu_contam ~ normal(0,3) ;
      //sigma_contam ~ gamma(2,2) ;
      target += normal_lpdf(sigma_contam | 0, 2) - 1 * normal_lccdf(0 | 0, 2);
      
      wash_offset ~ normal(wash_offset_prior[1],wash_offset_prior[2]) ;
      
      // Priors for smooth effects
      
      // priors including all constants
  //target += student_t_lpdf(Intercept | 3, 0, 2.5);
  
  target += normal_lpdf(b|0,5);
  target += normal_lpdf(bs|0,5);
  target += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_1_2);
  target += std_normal_lpdf(zs_1_3);
  target += student_t_lpdf(sds_2_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_2_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_2_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_2_1);
  target += std_normal_lpdf(zs_2_2);
  target += std_normal_lpdf(zs_2_3);
  target += student_t_lpdf(sds_3_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_3_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_3_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_3_1);
  target += std_normal_lpdf(zs_3_2);
  target += std_normal_lpdf(zs_3_3);
  target += student_t_lpdf(sds_4_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_4_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_4_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_4_1);
  target += std_normal_lpdf(zs_4_2);
  target += std_normal_lpdf(zs_4_3);
  target += student_t_lpdf(sds_5_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_5_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_5_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_5_1);
  target += std_normal_lpdf(zs_5_2);
  target += std_normal_lpdf(zs_5_3);
  target += student_t_lpdf(sds_6_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_6_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_6_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_6_1);
  target += std_normal_lpdf(zs_6_2);
  target += std_normal_lpdf(zs_6_3);
  target += student_t_lpdf(sds_7_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_7_1);
}
generated quantities {
  // log_likelihoods for use with the loo package a
  vector[N_obs_bin] log_lik;
  vector[N_obs_bin] pred_bin;
  vector[N_obs_pos] pred_pos;
  
  
    for(i in 1:N_obs_bin){
      pred_bin[i] = phi_0[pcr_obs_bin_idx[i]] + 
                      phi_1[pcr_obs_bin_idx[i]] * (D_contam[sample_bin_idx[i]] + bin_log_dilution_obs[i]) ;
      log_lik[i]  = bernoulli_logit_lpmf(bin_obs[i] | pred_bin[i]) ;
    }
    for(i in 1:N_obs_pos){
      pred_pos[i] = beta_0[pcr_obs_pos_idx[i]] + 
                      beta_1[pcr_obs_pos_idx[i]] * (D_contam[sample_pos_idx[i]] + pos_log_dilution_obs[i]) ;
      log_lik[loo_pos_idx[i]]  =  log_lik[loo_pos_idx[i]] + 
                                             student_t_lpdf(pos_obs[i] |3, pred_pos[i], sigma_all_samp) ;
    }
  
}