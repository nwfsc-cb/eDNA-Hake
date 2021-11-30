functions {
  // method-of-moments function
  vector mm_beta(real Mean, real SD){
    real Var = SD^2 ;
    real alpha = Mean * ((Mean*(1-Mean) /Var) -1) ; 
    real beta = alpha * (1-Mean) / Mean ;
    vector[2] out ;
      out[1] = alpha;
      out[2] = beta; 
  return(out);
  }

}

data {
  int<lower=0> I_master;
  int<lower=0> I_count;
  int<lower=0> I_mifish;
  int<lower=0> I_sebastes;

  int<lower=0> K_master ;
  int<lower=0> K_mifish ;
  int<lower=0> K_sebastes ;
  int<lower=0> K_count ;

  int<lower=0> J_mifish;
  int<lower=0> J_mifish_vec[K_mifish];

  int<lower=0> J_sebastes;
  int<lower=0> J_sebastes_vec[K_sebastes];

  // Number of species-station combinations to estimate
  int<lower=0> N_station_species_master;
  
  // Observations
  int<lower=0> N_count_obs;
  int<lower=0> N_mifish_obs;
  int<lower=0> N_sebastes_obs;
  
  int<lower=0> D_mifish_obs[N_mifish_obs];
  int<lower=0> D_sebastes_obs[N_sebastes_obs];
  int<lower=0> D_count_obs[N_count_obs];
  
  //  covariates / offsets 
  int<lower=0> N_mifish_station_rep_idx;
  int<lower=0> N_sebastes_station_rep_idx;
  vector[N_mifish_station_rep_idx] log_r_mifish;
  vector[N_sebastes_station_rep_idx] log_r_sebastes;
  real N_pcr; // Number of PCR cycles (assumed constant across all PCRs)
  
  //Count covariates / offsets
  vector[N_count_obs] log_V; // Volume of water sampled relative to the reference volume.
  vector[N_count_obs] log_P; // Proportion of jar identified

  // Indexes
  int<lower=0> count_station_idx[N_count_obs];
  int<lower=0> count_sp_idx[N_count_obs];

  int<lower=0> mifish_station_idx[N_mifish_obs];
  int<lower=0> mifish_sp_idx[N_mifish_obs]; // this is the index to the list of species only found by mifish primer
  int<lower=0> mifish_station_rep_idx[N_mifish_obs];

  int<lower=0> sebastes_station_idx[N_sebastes_obs];
  int<lower=0> sebastes_sp_idx[N_sebastes_obs]; // this is the index to the list of species only found by sebastes primer
  int<lower=0> sebastes_station_rep_idx[N_sebastes_obs];
  
  int<lower=0> master_station_idx[N_station_species_master] ;
  int<lower=0> master_sp_idx[N_station_species_master] ;

  // Mapping matrices
  matrix[I_count,I_master] M_to_count ;
  matrix[I_mifish,I_master] M_to_mifish ; 
  matrix[I_sebastes,I_master] M_to_sebastes ;
  
  real log_eta_prior[2]; 
  //real log_eta_mifish_fixed; 

  //real log_eta_sebastes_prior[2]; 
  
  //real a_sebastes; //fix sebastes amp efficiencies, for now
  
  real beta_prior_mifish[2];
  //real beta_prior_sebastes[2];
  
  vector[N_station_species_master] b_master_prior_log_mean;
  
  real<lower=0> a_sebastes_sd;
}

transformed data{
}

parameters {
  // MiFish specific parameters first
   real<lower = 0> tau ;
   //real<lower = 0> tau_sebastes ;
   //real<lower=0> phi;
  
   real<upper=0> log_eta_mifish ; // This is the fraction of total amplicons that is read by the sequencer (a scalar)
   real<upper=0> log_eta_sebastes ;
   //real<upper=0> log_eta_sebastes ; // This is the fraction of total amplicons that is read by the sequencer (a scalar)
   //real<lower = 0> tau; //standard deviation for normal distribution, representing error on the pipetting fraction between step 1 and step 2 PCR
   real<lower=0,upper=1> a_sebastes_all;
   
   vector<lower=0,upper=1>[I_mifish] a_mifish ; // amp efficiency of primer x species
   //vector<lower=0,upper=1>[I_sebastes] a_sebastes ; // amp efficiency of primer x species
   //vector<upper = 0>[N_station_rep_idx] log_r_mifish_real; //library prep; fraction of reaction from PCR2 into library; log scale

  // Master list log-intercept parameters (includes all species, all communities with non-zero observations)
  vector[N_station_species_master] b_master ;// This is a vector possible species / species groups.
  ///real<upper = 0> epsilon[N_mifish_obs]; //realized pipetting fraction PCR1 into PCR2; on log scale
  
  
  
}

transformed parameters{
   vector[I_master] b_master_grid[K_master];
   vector[I_mifish] b_mifish[K_mifish] ;
   vector[I_sebastes] b_sebastes[K_sebastes] ;
   vector[I_count] b_count[K_count] ;
   //vector<lower=0>[2] a_sebastes_beta_params; 

  //transform mean and SD for a_mifish to beta distribution
  //  a_sebastes_beta_params = mm_beta(a_sebastes_mean, a_sebastes_sd);

  for(k in 1:K_count){
    b_master_grid[k] = rep_vector(0,I_master); 
  }
  
  for(j in 1:N_station_species_master){
      b_master_grid[master_station_idx[j],master_sp_idx[j]] = exp(b_master[j]);
  }

  for(k in 1:K_count){
    b_count[k] = M_to_count * b_master_grid[k] ;
  }
  
  for(k in 1:K_mifish){
    b_mifish[k] = (M_to_mifish * b_master_grid[k]) / 
                  sum((M_to_mifish * b_master_grid[k])) ;
  }
  
  // print("master",b_master_grid[1,31]);
  // print("mifish",b_mifish[1,28]);
  
  for(k in 1:K_sebastes){
    b_sebastes[k] = (M_to_sebastes * b_master_grid[k]) / 
                    sum(M_to_sebastes * b_master_grid[k]) ;
  }
}

model {
  { // Local variables declaration for making the Stan program less Memory hungry
    vector[N_mifish_obs] log_lambda_mifish ;
    vector[N_sebastes_obs] log_lambda_sebastes ;
    vector[N_count_obs] log_theta ;

    for(q in 1:N_mifish_obs){
      log_lambda_mifish[q] =  -1.609438 + //log(0.2); the fraction of each sample subsampled. 
                              log(b_mifish[mifish_station_idx[q],mifish_sp_idx[q]]) + 
                              N_pcr*log(1.0 + a_mifish[mifish_sp_idx[q]]) +
                              log_r_mifish[mifish_station_rep_idx[q]] + // add a vector of known sample fractions
                              log_eta_mifish 
                              //epsilon[q] //pipetting fraction of PCR1 into PCR2, normal centered at log(.2), w sd = tau
                              ; 
    }
    
    
    for(q in 1:N_sebastes_obs){
      log_lambda_sebastes[q] =  -1.609438 + //log(0.2); the fraction of each sample subsampled. 
                              log(b_sebastes[sebastes_station_idx[q],sebastes_sp_idx[q]]) + 
                              N_pcr*log(1.0 + a_sebastes_all) +
                              log_r_sebastes[sebastes_station_rep_idx[q]] + // add a vector of known sample fractions
                              log_eta_sebastes 
                              //epsilon[q] //pipetting fraction of PCR1 into PCR2, normal centered at log(.2), w sd = tau
                              ; 
    }
    
    for(q in 1:N_count_obs){
        log_theta[q] = log(b_count[count_station_idx[q],count_sp_idx[q]]);
    }
        log_theta = log_theta + log_V + log_P;
        
    // Likelihoods.
    D_mifish_obs ~ neg_binomial_2(exp(log_lambda_mifish),tau);
    D_sebastes_obs ~ neg_binomial_2(exp(log_lambda_sebastes),tau);
    //D_sebastes_obs ~ poisson_log(log_lambda_sebastes);
    D_count_obs ~ poisson_log(log_theta);
    //epsilon ~ normal(0, tau); // THIS HAS
    } // end local variable declaration
    
    // Priors
    for (i in 1:N_station_species_master){
      b_master[i] ~ normal(b_master_prior_log_mean[i],4);  //not sure if the normal sampling function is vectorized. Making sure w a loop.
    }
    
    
    //log_r_mifish_real ~ normal(log_r_mifish,0.005); 
    target += normal_lpdf(log_eta_mifish | log_eta_prior[1], log_eta_prior[2]) - 1 * normal_lcdf(1 | log_eta_prior[1], log_eta_prior[2]);
    target += normal_lpdf(log_eta_sebastes | log_eta_prior[1], log_eta_prior[2]) - 1 * normal_lcdf(1 | log_eta_prior[1], log_eta_prior[2]);
        
    a_mifish ~ beta(beta_prior_mifish[1],beta_prior_mifish[2]);
    a_sebastes_all ~ beta(beta_prior_mifish[1],beta_prior_mifish[2]);
    //a_sebastes_all ~ beta(a_sebastes_beta_params[1], a_sebastes_beta_params[2]);   //beta(beta_prior_sebastes[1],beta_prior_sebastes[2]);
    //target += normal_lpdf(phi | 0, 1) - 1 * normal_lccdf(0 | 0, 1);
    tau ~ gamma(2, 5); 
    //tau_sebastes ~ gamma(10,10);  //changing this to 10,10 has no discernable effect
}
generated quantities{
  /// HERE MAKE DERIVED QUANTITIES FOR EACH PRIMER / MANUAL COUNT and 
   // vector[I_master] pi_master_grid[K_master];
   // vector[I_mifish] pi_mifish_grid[K_mifish];
   // vector[I_sebastes] pi_sebastes_grid[K_sebastes];

// Need to un-log the below, after moving b_master_grid out of log scale:

   // for(k in 1:K_master){
   //   pi_master_grid[k] = b_master_grid[k] - log_sum_exp(b_master_grid[k])  ;
   // }
   //  pi_master_grid = exp(pi_master_grid);

  
   // for(k in 1:K_mifish){
   //   pi_mifish_grid[k] = M_to_mifish * b_master_grid[k] + N_pcr*log(1+a_mifish) - 
   //                        log_sum_exp(M_to_mifish * b_master_grid[k] + N_pcr*log(1+a_mifish)) ;
   // }
   //  pi_mifish_grid = exp(pi_mifish_grid);
    
   //  for(k in 1:K_sebastes){
   //   pi_sebastes_grid[k] = M_to_sebastes * b_master_grid[k] + N_pcr*log(1+a_sebastes) - 
   //                        log_sum_exp(M_to_sebastes * b_master_grid[k] + N_pcr*log(1+a_sebastes)) ;
   // }
   //  pi_sebastes_grid = exp(pi_sebastes_grid);
    
}

