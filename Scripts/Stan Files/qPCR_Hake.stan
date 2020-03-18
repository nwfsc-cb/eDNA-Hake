// This is a model file for calculating the qPCR results for Skagit eDNA sampling in 2017.

data { /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Number of observations in various categories
    int N_sample; // Number of unique bottles used combinations observed.
    int N_control_sample; // Number of samples used for negative controls.
    int N_station_depth; // Number of station-depth combinations observed.
    int N_pcr ;    // Number of PCR plates
    
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
    vector[N_obs_bin] bin_vol_obs ;
    vector[N_obs_pos] pos_vol_obs ;
    vector[N_control_bin] bin_vol_control ;
    vector[N_control_pos] pos_vol_control ;

    // Standard Indices
    int pcr_stand_bin_idx[N_stand_bin] ;
    int pcr_stand_pos_idx[N_stand_pos] ;
    int pcr_obs_bin_idx[N_obs_bin] ;
    int pcr_obs_pos_idx[N_obs_pos] ;
    int pcr_control_bin_idx[N_control_bin] ;
    int pcr_control_pos_idx[N_control_pos] ;

    // Station-depth and sample indices
    int sample_idx[N_sample]    ;
    int sample_control_idx[N_control_sample]    ;
    int station_depth_idx[N_station_depth] ;

    // Sample related indices
    int sample_bin_idx[N_obs_bin]      ;
    int sample_pos_idx[N_obs_pos]  ;
    int sample_control_bin_idx[N_control_bin]      ;
    int sample_control_pos_idx[N_control_pos]  ;
    
    int station_depth_bin_idx[N_obs_bin]      ;
    int station_depth_pos_idx[N_obs_pos]  ;
  
    real singleton_bin_idx[N_obs_bin] ;
    real singleton_pos_idx[N_obs_pos] ;
    
    //Offset
    real OFFSET ;
}
transformed data{
}
parameters { /////////////////////////////////////////////////////////////////////////////////////////////
    // Standards regression and logit coeffs

      real beta_0[N_pcr] ;
      real beta_1[N_pcr] ;
      
      real phi_0[N_pcr] ;
      real phi_1[N_pcr];

    // Variance parameters for observation and random effects
      real sigma_stand_int ;
      //real sigma_stand_slope ;
     // real sigma_stand_slope2 ;
      
      // real sigma_stand_int_bar ;
      // real<lower=0> sigma_stand_int_sd ;
      // real sigma_stand_slope_bar ;
      // real<lower=0> sigma_stand_slope_sd ;
      // 
      //vector<lower=0>[N_pcr] sigma_stand_int ;
      
       real<lower=0> tau_sample ;
       real<lower=0> sigma_pcr ; 
    // 
    // // Effects of station-depths and samples
       real D[N_station_depth] ;
       real D_control[N_control_sample] ;
       real delta[N_sample] ;
}
transformed parameters { ////////////////////////////////////////////////////////////////////////////////
    // Latent variables for each log Density of 
      // real<lower=0> beta_0_sd ;
      // real<lower=0> beta_1_sd ;
      
     // vector[N_station_depth] D ;
      vector[N_stand_bin] theta_stand ;
      vector[N_obs_bin] theta_obs ;
      vector[N_control_bin] theta_control ;

      vector[N_stand_pos] kappa_stand ;
      vector[N_obs_pos] kappa_obs ;
      vector[N_control_pos] kappa_control ;
//      vector[N_count_stand] sigma_all_stand ;
  //    vector[N_count_samp] sigma_all_samp ;
      real sigma_all_stand ;
      real sigma_all_samp ;
      
      // vector[N_pcr] beta_1 ;
      // vector[N_pcr] beta_0 ;

      //MATT TRICK FOR BETA_1
      // beta_0 = beta_0_bar + beta_0_sd * beta_0_raw ;
      // beta_1 = beta_1_bar + beta_1_sd * beta_1_raw ;

      // beta_0_sd = exp(log_beta_0_sd) ;
      // beta_1_sd = exp(log_beta_1_sd) ;
  
    // Presence-Absence component of model.
    for(i in 1:N_stand_bin){
       theta_stand[i] = phi_0[pcr_stand_bin_idx[i]] + phi_1[pcr_stand_bin_idx[i]] *  (D_bin_stand[i] - OFFSET) ;
      //theta_stand[i] = phi_0 + phi_1 * (D_bin_stand[i] - OFFSET) ;

    }
    for(i in 1:N_obs_bin){
       theta_obs[i]  = phi_0[pcr_obs_bin_idx[i]] + 
                      phi_1[pcr_obs_bin_idx[i]] * (D[station_depth_bin_idx[i]] + bin_vol_obs[i] + singleton_bin_idx[i] * delta[sample_bin_idx[i]] - OFFSET) ;
      //theta_samp[i]  = phi_0 + phi_1 * (D[bottle_bin_idx[i]] - OFFSET) ;
    }

    for(i in 1:N_control_bin){
       theta_control[i]  = phi_0[pcr_control_bin_idx[i]] + 
                      phi_1[pcr_control_bin_idx[i]] * (D_control[sample_control_bin_idx[i]] + bin_vol_control[i] - OFFSET) ;
      //theta_samp[i]  = phi_0 + phi_1 * (D[bottle_bin_idx[i]] - OFFSET) ;
    }


    // Positive Comonent of the model
    for(i in 1:N_stand_pos){
       kappa_stand[i] = beta_0[pcr_stand_pos_idx[i]] + beta_1[pcr_stand_pos_idx[i]] * (D_pos_stand[i] - OFFSET) ;
      //kappa_stand[i] = beta_0 + beta_1 * D_count_stand[i] ;
    }
    for(i in 1:N_obs_pos){
      kappa_obs[i]  = beta_0[pcr_obs_pos_idx[i]] + 
                        beta_1[pcr_obs_pos_idx[i]] * (D[station_depth_pos_idx[i]] + pos_vol_obs[i] + singleton_pos_idx[i] * delta[sample_pos_idx[i]] - OFFSET) ;
      //kappa_samp[i]  = beta_0 + beta_1 * D[bottle_count_idx[i]] ;
    }
    for(i in 1:N_control_pos){
      kappa_control[i]  = beta_0[pcr_control_pos_idx[i]] + 
                        beta_1[pcr_control_pos_idx[i]] * (D_control[sample_control_pos_idx[i]] + pos_vol_control[i] - OFFSET) ;
      //kappa_samp[i]  = beta_0 + beta_1 * D[bottle_count_idx[i]] ;
    }
    
    // 
    //for( i in 1:N_count_stand){
      sigma_all_stand = sigma_stand_int;
                                //+ sigma_stand_slope * (D_count_stand[i]- OFFSET)),-2)   ;
                                //+ sigma_stand_slope2 * D_count_stand[i]^2),-2) ;

    //}  
    //for( i in 1:N_count_samp){
      sigma_all_samp = //sigma_stand_int;
                        pow(sigma_stand_int^2 +
                            //sigma_stand_slope * (D[bottle_count_idx[i]]- OFFSET)) +
                            //+ sigma_stand_slope2 * D[bottle_count_idx[i]]^2) + 
                            sigma_pcr^2,-2) ;
    //} 
    
}
model {////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Likelihood components
    bin_stand  ~ bernoulli( inv_logit(theta_stand) ) ;
    bin_obs   ~ bernoulli( inv_logit(theta_obs) ) ;
    bin_control   ~ bernoulli( inv_logit(theta_control) ) ;
    
    pos_stand ~ normal(kappa_stand, sigma_all_stand) ;
    pos_obs ~ normal(kappa_obs, sigma_all_samp) ;
    pos_control ~ normal(kappa_control, sigma_all_samp) ;

    // Parameters for the gBLOCKS:
     // Random effects
      D ~ normal(0,4);
      D_control ~ normal(0,4);
      // Priors
      sigma_stand_int ~ gamma(1,1) ;
      sigma_pcr ~ gamma(1,1) ;

      beta_0 ~ normal(40,5) ;
      beta_1 ~ normal(-3.32,0.1) ;

      phi_0 ~ normal(2, 2) ;
      phi_1 ~ normal(4, 2) ;

      delta ~ normal(0,tau_sample) ;
      tau_sample ~ gamma(1,1) ;

}
