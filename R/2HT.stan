data {
  int Nsubs;
  int corrects[3, 3, Nsubs];
  int trials[3, 3, Nsubs];
}

parameters {
  vector<lower=0, upper=1>[6] beta_mean;
  vector<lower=0>[6] beta_pres;
  vector<lower=0, upper=1>[Nsubs] DN;
  vector<lower=0, upper=1>[Nsubs] DO_W;
  vector<lower=0, upper=1>[Nsubs] DO_S;
  vector<lower=0, upper=1>[Nsubs] GO_Lib;
  vector<lower=0, upper=1>[Nsubs] GO_Neut;
  vector<lower=0, upper=1>[Nsubs] GO_Cons;
}

model {
  vector[6] a;
  vector[6] b;
  vector[9] p_correct[Nsubs];
  vector[3] guess_old;
  
  beta_mean ~ beta(5,5);
  beta_pres ~ gamma(1, 1);
  a = beta_mean .* beta_pres;
  b = (1-beta_mean) .* beta_pres;
  
  DN ~ beta(a[1], b[1]);
  DO_W ~ beta(a[2], b[2]);
  DO_S ~ beta(a[3], b[3]);
  GO_Lib ~ beta(a[4], b[4]);
  GO_Neut ~ beta(a[5], b[5]);
  GO_Cons ~ beta(a[6], b[6]);

  for (sub in 1:Nsubs) {

    guess_old = [GO_Lib[sub], GO_Neut[sub], GO_Cons[sub]]';
    p_correct[sub, 1:3] = DN[sub] + ((1-DN[sub]) * (1-guess_old));
    p_correct[sub, 4:6] = DO_W[sub] + (1-DO_W[sub]) * guess_old;
    p_correct[sub, 7:9] = DO_S[sub] + (1-DO_S[sub]) * guess_old;

    corrects[1:, 1, sub] ~ binomial(trials[1:, 1, sub], p_correct[sub, 1:3]);
    corrects[1:, 2, sub] ~ binomial(trials[1:, 2, sub], p_correct[sub, 4:6]);
    corrects[1:, 3, sub] ~ binomial(trials[1:, 3, sub], p_correct[sub, 7:9]);
  }
  
  
}
