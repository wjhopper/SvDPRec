data {
  int Nsubs;
  int yeses[3, 3, Nsubs];
  int trials[3, 3, Nsubs];
}

parameters {
  vector[2] mu_pop; // population d-primes
  vector<lower=0>[2] mu_pop_SD; // variablility of population d-primes
  
  real log_sigma_pop; // evidence noise in pop.
  real<lower=0> log_sigma_pop_SD; // variability of evidence noise in pop.
  
  ordered[3] crit_pop; //criterion locations in pop.
  vector<lower=0>[3] crit_pop_SD; //variability of criterion locations in pop.
  
  matrix[Nsubs, 2] mu_sub;
  vector[Nsubs] log_sigma_sub;
  ordered[3] crit_sub[Nsubs];
}

transformed parameters {
  vector<lower=0>[Nsubs] sigma_sub;
  for (sub in 1:Nsubs) {
    sigma_sub[sub] = exp(log_sigma_sub[sub]);
  }
}

model {
  vector[9] p[Nsubs]; // Yes-rates, arranged first by category, then by bias level

  // Population level priors
  mu_pop ~ normal(0, 3); // prior for population d-primes
  mu_pop_SD ~ normal(0, 3); // prior for variabilty in population d-primes
  log_sigma_pop ~ normal(0, 3); // prior for population level varaibility of the evidence strength distributions
  log_sigma_pop_SD ~ normal(0, 3); // prior for uncertainty in level varaibility of the evidence strength distributions
  crit_pop ~ normal(0, 3); // prior for population level criterion placements
  crit_pop_SD  ~ normal(0, 3); // prior for variability population level criterion placements
  
// aside: We have priors on both means and SD's here because we are stating the belief that the individual
// subjects have parameters that are drawn from one population-level distribution, but we are somewhat uncertain
// about the nature of that population level distribution. So, we put a prior on the population mean (a range of
// likely values for the mean) and the population SD (a range of likely values for the SD), and allow the
// individual subjects to have parameter values from this "uncertain" population level distribution.


  for (sub in 1:Nsubs) {

    // Subject level priors
    mu_sub[sub, 1:2] ~ normal(mu_pop, mu_pop_SD);
    log_sigma_sub[sub] ~ normal(log_sigma_pop, log_sigma_pop_SD);
    crit_sub[sub, 1:3] ~ normal(crit_pop, crit_pop_SD);
    
    p[sub, 1:3] = 1 - Phi(crit_sub[sub, 1:3]); // False Alarm Rates
    p[sub, 4:6] = 1 - Phi((crit_sub[sub, 1:3] - mu_sub[sub, 1])/sigma_sub[sub]); // Weak Hits Rates
    p[sub, 7:9] = 1 - Phi((crit_sub[sub, 1:3] - mu_sub[sub, 2])/sigma_sub[sub]); // String Hit Rates

    yeses[1:, 1, sub] ~ binomial(trials[1:, 1, sub], p[sub, 1:3]);
    yeses[1:, 2, sub] ~ binomial(trials[1:, 2, sub], p[sub, 4:6]);
    yeses[1:, 3, sub] ~ binomial(trials[1:, 3, sub], p[sub, 7:9]);
  }

}
