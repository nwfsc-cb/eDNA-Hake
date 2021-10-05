//
// This Stan program defines a simple model, of
// the estimated proportions of eDNA from certain 
// groups of fish
//
// The
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'. 
// We only have one 
data {
  // species
  //int<lower=0> I_master;
  int<lower=0> I_fish;//n_distinct_morphogroups
  //sites
 // int<lower=0> K_master ;
  int<lower=0> K_fish ;//n_distinct_stations(minus tech)
  //replicates
  int<lower=0> J_fish; // max replication level
  int<lower=0> J_fish_vec[K_fish];
  
  // Number of species-station combinations to estimate
  int<lower=0> N_station_species;
  // Observations
  int<lower=0> N_fish_obs;
  int<lower=0> D_fish_obs[N_fish_obs];
  
  //  covariates / offsets 
  
   int<lower=0> N_fish_station_rep_idx;
   vector[N_fish_station_rep_idx] log_r_mifish;
   real N_pcr_mifish;
   
   // Indexes
#  int<lower=0> count_station_idx[N_count_obs];
  int<lower=0> count_sp_idx[N_count_obs];
  
  int<lower=0> mifish_station_idx[N_mifish_obs];
  int<lower=0> mifish_sp_idx[N_mifish_obs]; // this is the index to the list of species only found by mifish primer
  int<lower=0> mifish_station_rep_idx[N_mifish_obs];
  
  int<lower=0> master_station_idx[N_station_species_master] ;
  int<lower=0> master_sp_idx[N_station_species_master] ;
  
  
  // Mapping matrices
  #matrix[I_mifish,I_master] M_to_mifish ; 
  real log_eta_prior[2]; 
  real beta_prior_mifish[2];
  
}

transformed data{
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
   real<lower = 0> tau_mifish;
 
   real<upper = 0>log_eta_mifish_mean;
 
   vector<upper=0>[N_mifish_station_rep_idx] log_eta_mifish ; // This is the fraction of total amplicons that is read by the sequencer (a scalar)

   vector<lower=0,upper=1>[I_mifish] a_mifish ; // amp efficiency of primer x species for the mifish primers

  // Master list log-intercept parameters (includes all species, all communities with non-zero observations)
  vector[N_station_species_master] b_master ;// This is a vector possible species / species groups.
}

transformed parameters{
#   vector[I_master] b_master_grid[K_master];
   vector[I_mifish] b_mifish[K_mifish] ;
   

  
  for(j in 1:N_station_species_master){
      b_master_grid[master_station_idx[j],master_sp_idx[j]] = exp(b_master[j]);
  }

  for(k in 1:K_mifish){
    b_mifish[k] = (M_to_mifish * b_master_grid[k]) / 
                  sum((M_to_mifish * b_master_grid[k])) ;
  }
  
 
}
// The model to be estimated. We model the output


model {
  
  { // Local variables declaration for making the Stan program less Memory hungry
    vector[N_mifish_obs] log_lambda_mifish ;

    for(q in 1:N_mifish_obs){
      log_lambda_mifish[q] =  -1.609438 + //log(0.2); the fraction of each sample subsampled, pipetting. 
                              log(b_mifish[mifish_station_idx[q],mifish_sp_idx[q]]) + 
                              N_pcr_mifish*log(1.0 + a_mifish[mifish_sp_idx[q]]) +
                              log_r_mifish[mifish_station_rep_idx[q]] + // add a vector of known sample fractions
                              log_eta_mifish[mifish_station_rep_idx[q]] // fraction of amplicons getting sequenced into reads
                              ; 
    }
    
        
    // Likelihoods.
    D_mifish_obs ~ neg_binomial_2(exp(log_lambda_mifish), tau_mifish);  //tau

    log_eta_mifish ~ normal(log_eta_mifish_mean, 0.75);
    
    } // end local variable declaration
    
       
    b_master ~ normal(0,4) ;
    
    target += normal_lpdf(log_eta_mifish_mean | log_eta_prior[1], log_eta_prior[2]) - 1 * normal_lcdf(1 | log_eta_prior[1], log_eta_prior[2]);
 
    a_mifish ~ beta(beta_prior_mifish[1],beta_prior_mifish[2]);
    tau_mifish ~ gamma(5, 5);
   
   
}

