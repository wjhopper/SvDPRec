data {
  int Nsubs;
  int yeses[3, 3, Nsubs];
  int trials[3, 3, Nsubs];
}

parameters {
  vector[2] mu;
  real log_sigma;
  ordered[3] crit;
}

transformed parameters {
  real<lower=0> sigma;
  sigma = exp(log_sigma);
}

model {
  vector[9] p; // Yes-rates, arranged first by category, then by bias level
  p = rep_vector(0.0, 9);

  mu ~ normal(0, 3);
  log_sigma ~ normal(0, 1);
  crit ~ normal(0, 3);

  for (sub in 1:Nsubs) {
    p[1:3] = 1 - Phi(crit); // False Alarm Rates
    p[4:6] = 1 - Phi((crit-mu[1])/sigma); // Weak Hits Rates
    p[7:9] = 1 - Phi((crit-mu[2])/sigma); // String Hit Rates
    
    yeses[1:, 1, sub] ~ binomial(trials[1:, 1, sub], p[1:3]);
    yeses[1:, 2, sub] ~ binomial(trials[1:, 2, sub], p[4:6]);
    yeses[1:, 3, sub] ~ binomial(trials[1:, 3, sub], p[7:9]);
  }
  print("Criterion: ", crit)
  print("Yeses: ", yeses[1:, 3, 14])
}
