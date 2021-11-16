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

// The input data is a series of vectors - the first one is the observed abundances of each species
// in each replicate. We want to estimate the relative abundace of each species on each point in space
// so we need to map each replicate to their point of origin

// We need to know how many are there I species, J stations, N replicates, IJ combinations
// We need to 


data {
  // species
  int<lower=0> I_taxa; // Number of taxa to estimate
  int<lower=0> N_mifish_station_rep; // Number of unique units sequenced (samples times rep)
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
 
  real beta_prior_mifish[2];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
   real<upper = 0>log_eta_mifish_mean;
   real<lower = 0> phi;
   vector<upper=0>[N_mifish_station_rep] log_eta_mifish ; // This is the fraction of total amplicons that is read by the sequencer (a scalar)
   vector<lower=0,upper=1>[taxa_idx] a_mifish ;
}

transformed parameters{
vector[I_taxa][J_stations]b_mifish  
  
   for(k in 1:J_stations){
    b_mifish[k] = (M_to_mifish * b_master_grid[k]) /
                  sum((M_to_mifish * b_master_grid[k])) ; // everything here is in proportion space, I will probably want to declare b_mifish[] as a simplex for each - YOu don;t have all species on all samples
                  
                  // Move b to the parameters instead of the Bmaster. B points to the data not estimated - True Zeroes  or non-observation. It is more complicated. Assuming that there are there but rare and not detected. 
                  
  } 
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
    
    phi ~ normal(0, 2);
    a_mifish ~ beta(beta_prior_mifish[1],beta_prior_mifish[2]);
   
}

