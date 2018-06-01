data{
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 0, upper = 1> y[N, J];
  int<lower = 0, upper = 1> X1[N];
  int<lower = 0, upper = 1> X2[N];
  // real<lower=1> nu;
}

parameters{
  vector[N] theta;
  vector<lower = 0>[J] alpha; //discrimination
  vector[J] betas; // item difficulty
  real X1_indirect_eff; 
  real X2_indirect_eff; 
}

transformed parameters{
  vector[J] X1_direct_eff; 
  vector[J] X2_direct_eff;
  
  // global and local variance parameters, and the input weights
  real<lower=0> tau;
  vector<lower=0>[2] lambda1;
  vector<lower=0>[2] lambda2;
  
  tau = .5;
  X1_direct_eff[J] = 0;
  X2_direct_eff[J] = 0;
}

model{
  for (i in 1:J){
    betas[i] ~ normal(0, 4);
    alpha[i] ~ lognormal(0, 2.5);
  }
  
  lambda1 ~ cauchy(0, 1);
  lambda2 ~ cauchy(0, 1);
  
  X1_indirect_eff ~ normal(0, lambda1[1]*tau);
  X2_indirect_eff ~ normal(0, lambda2[1]*tau);
  
  X1_direct_eff[1:J-1] ~ normal(0, lambda1[2]*tau);
  X2_direct_eff[1:J-1] ~ normal(0, lambda2[2]*tau);

  for (i in 1:N){
    theta[i] ~ normal(X1[i]*X1_indirect_eff + X2[i]*X2_indirect_eff, 1);
  }
  
  for (i in 1:N){
    for (j in 1:J){
      y[i, j] ~ bernoulli_logit(alpha[j]*(theta[i] - betas[j] - X1_direct_eff[j]*X1[i] - X2_direct_eff[j]*X2[i]));
    }
  }
}