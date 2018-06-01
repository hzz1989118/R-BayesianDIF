data{
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 0, upper = 1> y[N, J];
  int<lower = 0, upper = 1> X1[N];
}

parameters{
  vector[N] theta;
  real<lower=0> alpha_1;
  real beta_1;
  vector<lower = 0>[J-1] alphas_ref;
  vector[J-1] betas_ref;
  matrix[J-1, 2] eff;
  //vector<lower = 0, upper = 1>[J] guesses;
  real eff_mu;
  real eff_sigma;
  matrix<lower = 0>[J-1, 2] tau2; // auxilliary parameter;
  real<lower = 0> lambda2;
}

transformed parameters{
  vector[J-1] eff_alpha;
  vector[J-1] eff_beta;
  eff_alpha = eff[, 1];
  eff_beta = eff[, 2];
}

model{
  alpha_1 ~ lognormal(0, 2.5);
  beta_1 ~ normal(0, 4);
  eff_sigma ~ normal(0, 4);
  lambda2 ~ gamma(1, 2);
  tau2[,1] ~ gamma(1, lambda2/2);
  tau2[,2] ~ gamma(1, lambda2/2);
  
  for (i in 1:(J-1)){
    betas_ref[i] ~ normal(0, 4);
    alphas_ref[i] ~ lognormal(0, 2.5);
    eff[i,1] ~ normal(0, sqrt(tau2[i, 1]));
    eff[i,2] ~ normal(0, sqrt(tau2[i, 2]));
  }

  for (i in 1:N){
    theta[i] ~ normal(0 + X1[i]*eff_mu, exp(log(1) + X1[i]*eff_sigma));
  }
  
  for (i in 1:N){
    y[i, 1] ~ bernoulli_logit(alpha_1*(theta[i] - beta_1));
    for (j in 1:(J-1)){
      y[i, j+1] ~ bernoulli_logit(exp(log(alphas_ref[j]) + X1[i]*eff_alpha[j])*(theta[i] - betas_ref[j] - X1[i]*eff_beta[j]));
    }
  }
}