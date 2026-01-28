data {
int<lower=1> N_obs;
vector[N_obs] bray; //bayes
// vector[N_obs] lai; // lai
// vector[N_obs] sumTemp;// temp
// int<lower=1> N_fetch;
// vector[N_obs] fetch;
int N_regions;
array[N_obs] int<lower=1, upper=N_regions> region_idx;
int N_samples;
array[N_obs] int<lower=1, upper=N_samples> sample_idx1;
array[N_obs] int<lower=1, upper=N_samples> sample_idx2;}

parameters {
real a_baseline;
vector[N_regions] a_region;
real<lower=0> sigma_region;

vector[N_samples] a_sample_ncp;
real<lower=0> sigma_sample;

// vector[N_regions] bfetch_ncp;
// real mu_bfetch;
// real<lower=0> sigma_bfetch;
// vector[N_regions] blai;
// real mu_blai;
// real<lower=0> sigma_blai;

// vector[N_regions] bsumTemp;
// real mu_bsumTemp;
// real<lower=0> sigma_bsumTemp;

real<lower=0> kappa;
}

transformed parameters {
 vector[N_samples] a_samples = sigma_sample * a_sample_ncp;
// 
// vector[N_regions] bfetch = mu_bfetch+sigma_bfetch*bfetch_ncp;

vector[N_obs] mu;

mu = a_baseline
+ a_samples[sample_idx1] + a_samples[sample_idx2]
+ a_region[region_idx];
// + bfetch[region_idx].* fetch;
// + blai[region_idx] .* lai
// + bsumTemp[region_idx] .* sumTemp;
}

model {
a_baseline ~ normal(0, 1);

a_sample_ncp ~ normal(0, 1);
sigma_sample ~ exponential(4);

// a_sample ~ normal(0, sigma_sample);
// sigma_sample ~ exponential(4);

a_region ~ normal(0, sigma_region);
sigma_region ~ exponential(4);
// blai~normal( mu_blai, sigma_blai );
// mu_blai~normal(0,0.3);
// sigma_blai ~ exponential(4);
// bfetch_ncp~normal(0,1);
// mu_bfetch~normal(0,0.3);
// sigma_bfetch ~ exponential(4);
// bsumTemp~normal(mu_bsumTemp,sigma_bsumTemp);
// mu_bsumTemp~normal(0,0.3);
// sigma_bsumTemp ~ exponential(4);
kappa ~ normal(0, 1);
bray ~ beta_proportion(inv_logit(mu), kappa);
}

generated quantities{
array[N_obs] real bray_pred;
array[N_obs] real res;
for(i in 1:N_obs){
bray_pred[i]= beta_proportion_rng(inv_logit(mu[i]), kappa);
}

for(i in 1:N_obs){
res[i]= bray_pred[i]-bray[i];
}
}
