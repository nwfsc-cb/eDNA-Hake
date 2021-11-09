//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  // species
  int<lower=0> I_taxa;
  // stations
  int<lower=0> J_stations;
  //replicates
  int<lower=0> K_reps; //max replication level 
  int<lower=0> K_reps_vector[J_stations]; // vector with the number of replicates for each stations
  //Combos
  int<lower=0> Number_station_taxa; # length of the vector of Combinations to estimate
  //Observations
  int<lower=0> Number_taxa_obs; # length of the vector of input values
  int<lower=0> D_taxa_obs[Number_taxa_obs]; # vector of input values
  // Covariates?
  
  // Indices
  // We need indices for species, stations and replicates
  
   int<lower=0> station_idx[Number_taxa_obs];
   int<lower=0> taxa_idx[Number_taxa_obs];
   int<lower=0> replicate_idx[Number_taxa_obs];
   
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
   real<upper = 0>log_eta_mifish_mean;
   real<lower = 0> phi;
   vector<upper=0>[N_mifish_station_rep_idx] log_eta_mifish ; // This is the fraction of total amplicons that is read by the sequencer (a scalar)
   vector<lower=0,upper=1>[taxa_idx] a_mifish ;
}

transformed parameters{
vector[]  
  
  
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
   { // Local variables declaration for making the Stan program less Memory hungry
    vector[Number_taxa_obs] log_lambda ;

    for(q in 1:Number_taxa_obs){
      log_lambda[q] =  -1.609438 + //log(0.2); the fraction of each sample subsampled, pipetting. 
                              log(b_mifish[station_idx[q],taxa_idx[q]]) + 
                              N_pcr_mifish*log(1.0 + a_mifish[taxa_idx[q]]) +
                              log_r_mifish[replicate_idx[q]] + // add a vector of known sample fractions
                              log_eta_mifish[replicate_idx[q]] // fraction of amplicons getting sequenced into reads
                              ; 
    }
    
    // Likelihoods
    
     D_taxa_obs ~ neg_binomial_2(exp(log_lambda), phi);  //phi

    log_eta_mifish ~ normal(log_eta_mifish_mean, 0.75);
    
    }
    //models 
    b_master ~ normal(0,4) ;
    
    a_mifish ~ beta(beta_prior_mifish[1],beta_prior_mifish[2]);
    tau_mifish ~ gamma(5, 5);
}

