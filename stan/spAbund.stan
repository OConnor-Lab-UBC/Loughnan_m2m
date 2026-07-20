data {
  int<lower=0> N;
  int<lower=0> Nspp;
array[N] int<lower=1, upper=Nspp> species;

int<lower=0> Nsite;
array[N] int<lower=1, upper=Nsite> site;
vector[N] year;
 array[N] int<lower=0> y;
}

parameters {
  // Global (Population) Effects
  real alpha_global;
  real beta_global;

  // Between-group variation (Standard Deviations)
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_alphaS;

  real<lower=0> sigma_beta;

  // Group-level parameters defined DIRECTLY
  vector[Nspp] alpha_group;
    vector[Nsite] alpha_groupS;

  vector[Nspp] beta_group;

  // Negative Binomial 
  real<lower=0> phi;
}

model {
  // 1. Priors for Global parameters
  alpha_global ~ normal(0, 5);
  beta_global ~ normal(0, 2);
  sigma_alpha ~ exponential(1);
  sigma_alphaS ~ exponential(1);
  sigma_beta ~ exponential(1);
  phi ~ exponential(1);

  // 2. Centered Distribution for Group Effects
  // These parameters are "centered" on the global mean and sigma
  alpha_group ~ normal(alpha_global, sigma_alpha);
  alpha_groupS ~ normal(alpha_global, sigma_alphaS);
  beta_group ~ normal(beta_global, sigma_beta);

  // 3. Likelihood
  // Using vectorization for efficiency
  y ~ neg_binomial_2_log(alpha_group[species] +alpha_groupS[site] + beta_group[species] .* year, phi);
}

generated quantities {
 array[N] int y_rep;
  for (n in 1:N) {
    real mu = exp(alpha_group[species[n]] + alpha_groupS[site[n]] + beta_group[species[n]] * year[n]);
    y_rep[n] = neg_binomial_2_rng(mu, phi);
  }
}

 
