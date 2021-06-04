// This is a model file for calculating the acoustic results for  sampling in 2019.

data { /////////////////////////////////////////////////////////////////////////////////////////////////////

    int N_obs_bin   ;  // Number of observations for binomial part of the sample model
    int N_obs_pos ;  // Number of observations for count part of the sample model

    // Observations
    int bin_weight_dens[N_obs_bin]      ; 
    vector[N_obs_pos] pos_weight_dens   ;

    // LOO Index
    // int loo_pos_idx[N_obs_pos] ;

  ///////////////////////////////////////////////////////////////  
  ///// SMOOTH COMPONTENTS, extracted from BRMS package
  ///////////////////////////////////////////////////////////////

  // data for splines
  int Ks;  // number of linear effects
  matrix[N_obs_bin, Ks] Xs;  // design matrix for the linear effects
  
  // data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)0
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N_obs_bin, knots_1[1]] Zs_1_1;
  matrix[N_obs_bin, knots_1[2]] Zs_1_2;
  matrix[N_obs_bin, knots_1[3]] Zs_1_3;
 
  // data for spline s(bottom.depth,k=N.knots.bd)
  // int nb_2;  // number of bases
  // int knots_2[nb_2];  // number of knots
  // basis function matrices
  // matrix[N_obs_pos, knots_2[1]] Zs_2_1;

  // Priors and offset on cloglog components
  real thresh_mt_nm2;
  
  real phi_0_fix ;  
  // real phi_0_mean ;
  // real<lower=0> phi_0_sd ;

  real phi_1_mean ;
  real<lower=0> phi_1_sd ;
}
transformed data{
  real log_thresh_mt_nm2;
  
  log_thresh_mt_nm2 = log(thresh_mt_nm2) ;
}
parameters { //////////////////////////////////////////////////////////////////////////
  // Overall spatial intercept.
      real Intercept ;

  // Variance parameter for observation 
    real<lower=0> sigma ;
       
    // Logit coeffs
      //real phi_0 ; Use phi_0_fix instead
    //real<lower=0> phi_1 ;

    // Effects of station-depths and samples

   ////// SMOOTHES ////////////////////////////
    vector[Ks] bs;  // spline coefficients
    // parameters for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)0
    // standarized spline coefficients
    vector[knots_1[1]] zs_1_1;
    vector[knots_1[2]] zs_1_2;
    vector[knots_1[3]] zs_1_3;
    real<lower=0> sds_1_1;  // standard deviations of spline coefficients
    real<lower=0> sds_1_2;  // standard deviations of spline coefficients
    real<lower=0> sds_1_3;  // standard deviations of spline coefficients

    // parameters for spline s(bottom.depth,k=N.knots.bd)
    // standarized spline coefficients
    // vector[knots_2[1]] zs_2_1;
    // real<lower=0> sds_2_1;  // standard deviations of spline coefficients
}
transformed parameters { ////////////////////////////////////////////////////////////////////////////////
    // Latent variable for contaminated log-concentration of DNA that is observed by PCR.
    //  vector[N_obs_bin] D;

    // actual spline coefficients
      vector[knots_1[1]] s_1_1;
      vector[knots_1[2]] s_1_2;
      vector[knots_1[3]] s_1_3;
      //vector[knots_2[1]] s_2_1;

    // compute actual spline coefficients
      s_1_1 = sds_1_1 * zs_1_1;
      s_1_2 = sds_1_2 * zs_1_2;
      s_1_3 = sds_1_3 * zs_1_3;
      //s_2_1 = sds_2_1 * zs_2_1;
}
model {////////////////////////////////////////////////////////////////////////////////////////////////////
    { //LOCAL VARIABLES DECLARATION START
        vector[N_obs_bin] D;
        //vector[N_obs_bin] theta_bin;
    // Presence-Absence component of model.
    // Positive Comonent of the model
      D =  Intercept + 
                    // X * b +  // Factor level effects
                    Xs * bs + // linear effects of smoothes
                    Zs_1_1 * s_1_1 + Zs_1_2 * s_1_2 + Zs_1_3 * s_1_3; // + 
                    //Zs_2_1 * s_2_1;
    // for(i in (N_obs_pos+1):N_obs_bin){
    //   if(D[i] < -7){
    //     D[i] = -7;
    //   }
    // }
        
        
    //theta_bin = inv_logit(phi_0_fix + 20 * (exp(D) - thresh_mt_nm2));
    
    // print("exp_D ", exp(D[(N_obs_pos+1)])) ;
    // print("theta_bin before ", theta_bin[(N_obs_pos+1)]) ;
    // 
    // print("pow1 ",(1- 0.00000001));
    // print("pow2 ",( 0.00000001));
    
    // for(i in 1:N_obs_bin){
    //   if(theta_bin[i] > (1- 0.00000001)){
    //      theta_bin[i] = (1- 0.00000001);
    //   }
      // if(theta_bin[i] < 0.00000001){
      //    theta_bin[i] = 0.00000001 ;
      // }
    //}
    // print("theta_bin after", theta_bin[(N_obs_pos+1)]) ;

    // Likelihood components
    //bin_weight_dens     ~ bernoulli( theta_bin) ;
    pos_weight_dens     ~ lognormal(D[1:N_obs_pos] - 0.5 * sigma^2, sigma) ;

    // print("D1 ",D[1:10]);
    // print("bernD1 ",theta_bin[1:10]);
    // print("log_posD1 ",log(pos_weight_dens[1:10]));
    // print("---");
    // print("D2 ",D[(N_obs_pos-1):(N_obs_pos+1)]);
    // print("bernD2 ",theta_bin[(N_obs_pos-1):(N_obs_pos+1)]);
    // print("log_posD2 ",log(pos_weight_dens[(N_obs_pos-1):(N_obs_pos)]));


    } //LOCAL VARIABLES DECLARATION END 

      // Priors
      //sigma_stand_int ~ gamma(1,1) ;
      Intercept ~ normal(0,3); 
      target += normal_lpdf(sigma | 0, 0.05) - 1 * normal_lccdf(0 | 0, 0.05);
      
      //phi_0 ~ normal(phi_0_mean,phi_0_sd) ;
      //phi_1 ~ normal(phi_1_mean,phi_1_sd) ;

      // Priors for smooth effects
      
      // priors including all constants
  //target += student_t_lpdf(Intercept | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += student_t_lpdf(sds_1_2 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += student_t_lpdf(sds_1_3 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_1_2);
  target += std_normal_lpdf(zs_1_3);

  // target += student_t_lpdf(sds_2_1 | 3, 0, 2) - 1 * student_t_lccdf(0 | 3, 0, 2);
  // target += std_normal_lpdf(zs_2_1);
}
generated quantities {
  // log_likelihoods for use with the loo package
  //
  // vector[N_obs_bin] log_lik;
  // 
  //   for(i in 1:N_obs_bin){
  //    log_lik[i]  = bernoulli_logit_lpmf(bin_obs[i] | phi_0[pcr_obs_bin_idx[i]] + 
  //                                                     phi_1[pcr_obs_bin_idx[i]] * (D_contam[sample_bin_idx[i]] + 
  //                                                                           bin_log_dilution_obs[i])) ;
  //   }
  //   for(i in 1:N_obs_pos){
  //     log_lik[loo_pos_idx[i]]  =  log_lik[loo_pos_idx[i]] + 
  //                                            student_t_lpdf(pos_obs[i] |3, beta_0[pcr_obs_pos_idx[i]] + 
  //                                                       beta_1[pcr_obs_pos_idx[i]] * (D_contam[sample_pos_idx[i]] + 
  //                                                                 pos_log_dilution_obs[i]), sigma_all_samp);
  //   }
  // 
}