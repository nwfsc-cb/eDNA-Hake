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
  int<lower=0> D_taxa_obs; # vector of input values
  // Covariates?
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(mu, sigma);
}

