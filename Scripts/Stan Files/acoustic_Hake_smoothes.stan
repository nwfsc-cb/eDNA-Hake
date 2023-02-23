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
  int Ks_pos;  // number of linear effects
  int Ks_bin;
  matrix[N_obs_pos, Ks_pos] Xs_pos;  // design matrix for the linear effects
  matrix[N_obs_bin, Ks_bin] Xs_bin;  // design matrix for the linear effects
  
  int nb_1_pos;  // number of bases
  int nb_1_bin;  // number of bases
  int knots_1_pos[nb_1_pos];  // number of knots
  int knots_1_bin[nb_1_bin];  // number of knots
  
  // basis function matrices
  matrix[N_obs_pos, knots_1_pos[1]] Zs_1_1_pos;
  matrix[N_obs_pos, knots_1_pos[2]] Zs_1_2_pos;
  matrix[N_obs_pos, knots_1_pos[3]] Zs_1_3_pos;

  matrix[N_obs_bin, knots_1_bin[1]] Zs_1_1_bin;
  matrix[N_obs_bin, knots_1_bin[2]] Zs_1_2_bin;
  matrix[N_obs_bin, knots_1_bin[3]] Zs_1_3_bin;

  // data for spline s(bottom.depth,k=N.knots.bd)
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N_obs_pos, knots_2[1]] Zs_2_1_pos;
  matrix[N_obs_bin, knots_2[1]] Zs_2_1_bin;

  // Priors and offset on cloglog components
  real thresh_mt_km2;
  
  real phi_0_fix ;  
  // real phi_0_mean ;
  // real<lower=0> phi_0_sd ;

  real phi_1_mean ;
  real<lower=0> phi_1_sd ;
}
transformed data{
  real log_thresh_mt_km2;
  
  log_thresh_mt_km2 = log(thresh_mt_km2) ;
}
parameters { //////////////////////////////////////////////////////////////////////////
  // Overall spatial intercept.
      real Intercept_pos ;
      real Intercept_bin ;
  // Variance parameter for observation 
     real<lower=0> sigma ;
       
    //real v_0;   
    //real v_1;
    
    // Logit coeffs
    //real phi_0 ; Use phi_0_fix instead
    //real<lower=0> phi_1 ;

    // Effects of station-depths and samples

   ////// SMOOTHES ////////////////////////////
    vector[Ks_pos] bs_pos;  // spline coefficients
    vector[Ks_bin] bs_bin;  // spline coefficients
    
    // parameters for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)0
    // standarized spline coefficients
    vector[knots_1_pos[1]] zs_1_1_pos;
    vector[knots_1_pos[2]] zs_1_2_pos;
    vector[knots_1_pos[3]] zs_1_3_pos;
    
    vector[knots_1_bin[1]] zs_1_1_bin;
    vector[knots_1_bin[2]] zs_1_2_bin;
    vector[knots_1_bin[3]] zs_1_3_bin;
    
    vector[knots_2[1]] zs_2_1_pos;
    vector[knots_2[1]] zs_2_1_bin;
    
    real<lower=0> sds_1_1_pos;  // standard deviations of spline coefficients
    real<lower=0> sds_1_2_pos;  // standard deviations of spline coefficients
    real<lower=0> sds_1_3_pos;  // standard deviations of spline coefficients

    real<lower=0> sds_1_1_bin;  // standard deviations of spline coefficients
    real<lower=0> sds_1_2_bin;  // standard deviations of spline coefficients
    real<lower=0> sds_1_3_bin;  // standard deviations of spline coefficients

    real<lower=0> sds_2_1_pos;  // standard deviations of spline coefficients
    real<lower=0> sds_2_1_bin;  // standard deviations of spline coefficients
    

    // parameters for spline s(bottom.depth,k=N.knots.bd)
    // standarized spline coefficients
    // vector[knots_2[1]] zs_2_1;
    // real<lower=0> sds_2_1;  // standard deviations of spline coefficients
}
transformed parameters { ////////////////////////////////////////////////////////////////////////////////
    // Latent variable for contaminated log-concentration of DNA that is observed by PCR.
    //  vector[N_obs_bin] D;

    // actual spline coefficients
      vector[knots_1_pos[1]] s_1_1_pos;
      vector[knots_1_pos[2]] s_1_2_pos;
      vector[knots_1_pos[3]] s_1_3_pos;
      vector[knots_2[1]] s_2_1_pos;

      vector[knots_1_bin[1]] s_1_1_bin;
      vector[knots_1_bin[2]] s_1_2_bin;
      vector[knots_1_bin[3]] s_1_3_bin;
      vector[knots_2[1]] s_2_1_bin;

      //vector[knots_2[1]] s_2_1;

    // compute actual spline coefficients
      s_1_1_pos = sds_1_1_pos * zs_1_1_pos;
      s_1_2_pos = sds_1_2_pos * zs_1_2_pos;
      s_1_3_pos = sds_1_3_pos * zs_1_3_pos;
      
      s_1_1_bin = sds_1_1_bin * zs_1_1_bin;
      s_1_2_bin = sds_1_2_bin * zs_1_2_bin;
      s_1_3_bin = sds_1_3_bin * zs_1_3_bin;

      s_2_1_pos = sds_2_1_pos * zs_2_1_pos;
      s_2_1_bin = sds_2_1_bin * zs_2_1_bin;
}
model {////////////////////////////////////////////////////////////////////////////////////////////////////
    { //LOCAL VARIABLES DECLARATION START
        vector[N_obs_pos] D_pos;
        vector[N_obs_bin] theta_bin;
        //vector[N_obs_bin] theta_bin;
        //vector[N_obs_pos] sigma2;
    
    // Presence-Absence component of model.
    // Positive Comonent of the model
      D_pos =  Intercept_pos + 
                    // X * b +  // Factor level effects
                    Xs_pos * bs_pos + // linear effects of smoothes
                    Zs_1_1_pos * s_1_1_pos + Zs_1_2_pos * s_1_2_pos + Zs_1_3_pos * s_1_3_pos + // + 
                    Zs_2_1_pos * s_2_1_pos;
      theta_bin =  Intercept_bin + 
                    // X * b +  // Factor level effects
                    Xs_bin * bs_bin + // linear effects of smoothes
                    Zs_1_1_bin * s_1_1_bin + Zs_1_2_bin * s_1_2_bin + Zs_1_3_bin * s_1_3_bin + // + 
                    Zs_2_1_bin * s_2_1_bin;
                    //Zs_2_1 * s_2_1;
    
    // for(i in (N_obs_pos+1):N_obs_bin){
    //   if(D[i] < -7){
    //     D[i] = -7;
    //   }
    // }
        
    //theta_bin = inv_logit(phi_0_fix + 20 * (exp(D_bin) - thresh_mt_nm2));
    //sigma = exp(0.182 + v_1 * D_pos) ;
    //sigma2 = sigma .* sigma ;
    // print("exp_D ", exp(D[(N_obs_pos+1)])) ;
    // print("theta_bin before ", theta_bin[(N_obs_pos+1)]) ;
    // 
    // print("pow1 ",(1- 0.00000001));
    // print("pow2 ",( 0.00000001));
    
    // for(i in 1:N_obs_bin){
    //    if(theta_bin[i] > (1- 0.00000001)){
    //       theta_bin[i] = (1- 0.00000001);
    //    }
    //   // if(theta_bin[i] < 0.00000001){
    //   //    theta_bin[i] = 0.00000001 ;
    //   // }
    // }
    // print("theta_bin after", theta_bin[(N_obs_pos+1)]) ;

    // Likelihood components
    bin_weight_dens     ~ bernoulli(inv_logit(theta_bin)) ;
    pos_weight_dens     ~ lognormal(D_pos - 0.5 * sigma^2, sigma) ;

    // print("D1 ",D_pos[1:10]);
    // print("bernD1 ",theta_bin[1:10]);
    // print("bernD2 ",inv_logit(theta_bin[1:10]));
    // print("log_posD1 ",log(pos_weight_dens[1:10]));
    // print("---");
    // print("D2 ",theta_bin[(N_obs_pos-1):(N_obs_pos+8)]);
    // print("bernD2 ",inv_logit(theta_bin[(N_obs_pos-1):(N_obs_pos+8)]));
    // print("bern_obs ",bin_weight_dens[(N_obs_pos-1):(N_obs_pos+8)]);
    // print("log_posD2 ",log(pos_weight_dens[(N_obs_pos-1):(N_obs_pos)]));

    } //LOCAL VARIABLES DECLARATION END 

      // Priors
      //sigma_stand_int ~ gamma(1,1) ;
      Intercept_pos ~ normal(0,3); 
      Intercept_bin ~ normal(-1,3); 
      target += normal_lpdf(sigma | 0, 1) - 1 * normal_lccdf(0 | 0, 1);
      
      bs_pos ~ normal(0,5);
      bs_bin ~ normal(0,5);
      
      //phi_0 ~ normal(phi_0_mean,phi_0_sd) ;
      //phi_1 ~ normal(phi_1_mean,phi_1_sd) ;

      // Priors for smooth effects
      
      // priors including all constants
  //target += student_t_lpdf(Intercept | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_1_pos | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_2_pos | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_3_pos | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1_pos);
  target += std_normal_lpdf(zs_1_2_pos);
  target += std_normal_lpdf(zs_1_3_pos);
  target += student_t_lpdf(sds_2_1_pos | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_2_1_pos);

  target += student_t_lpdf(sds_1_1_bin | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_2_bin | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sds_1_3_bin | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1_bin);
  target += std_normal_lpdf(zs_1_2_bin);
  target += std_normal_lpdf(zs_1_3_bin);
  target += student_t_lpdf(sds_2_1_bin | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_2_1_bin);

  // target += student_t_lpdf(sds_2_1 | 3, 0, 2) - 1 * student_t_lccdf(0 | 3, 0, 2);
  // target += std_normal_lpdf(zs_2_1);
}
generated quantities {
  //log_likelihoods for use with the loo package
  vector[N_obs_bin] log_lik_bin;
  vector[N_obs_pos] log_lik_pos;
  vector[N_obs_bin] theta_bin_pred;
  vector[N_obs_pos] D_pos_pred;
    
    theta_bin_pred =   Intercept_bin + 
                    Xs_bin * bs_bin + // linear effects of smoothes
                    Zs_1_1_bin * s_1_1_bin + Zs_1_2_bin * s_1_2_bin + Zs_1_3_bin * s_1_3_bin +
                    Zs_2_1_bin * s_2_1_bin;
  
    D_pos_pred     =    Intercept_pos + 
                    Xs_pos * bs_pos + // linear effects of smoothes
                    Zs_1_1_pos * s_1_1_pos + Zs_1_2_pos * s_1_2_pos + Zs_1_3_pos * s_1_3_pos +
                    Zs_2_1_pos * s_2_1_pos;

  for(i in 1:N_obs_bin){
    log_lik_bin[i]  = bernoulli_logit_lpmf(bin_weight_dens[i] | theta_bin_pred[i]) ;
  }
for(i in 1:N_obs_pos){
    log_lik_pos[i]  =  lognormal_lpdf(pos_weight_dens[i] | D_pos_pred[i] -
                                        0.5*sigma^2,
                                        sigma);
  }

}