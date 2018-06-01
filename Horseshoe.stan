
functions {
  // square root of a vector (elementwise)
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] out;
    for (m in 1:dims(x)[1]){
      out[m] = sqrt(x[m]);
    }
    return out;
  }
}

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
  vector[J-1] eff_alpha;
  vector[J-1] eff_beta;
  //vector<lower = 0, upper = 1>[J] guesses;
  real eff_mu;
  real eff_sigma;
  vector<lower=0>[2*J-2] lambdas;
}

transformed parameters{
  // vector<lower = 0>[J-1] alphas_foc;
  // vector<lower = 0>[J-1] ua;
  real<lower=0> tau;
  
  // for(i in 1:(J-1)){
  //   alphas_foc[i] = exp(log(alphas_ref[i]) + eff_alpha[i]);
  //   ua[i] = fabs(2*(alphas_foc[i] - alphas_ref[i])/(alphas_ref[i]*alphas_foc[i])*log(1 + exp(alphas_ref[i]*alphas_foc[i]*eff_beta[i]/(alphas_foc[i] - alphas_ref[i]))) - eff_beta[i]);
  // }
  // 
  tau = .5;
 
  // for(i in 1: (2*J-2)){
  //   sigma_[i] = tau*lambdas[i];
  // }
  
}

model{
  alpha_1 ~ lognormal(0, 2.5);
  beta_1 ~ normal(0, 4);
  for (i in 1:(J-1)){
    betas_ref[i] ~ normal(0, 4);
    eff_beta[i] ~ normal(0, tau*lambdas[i]);
    alphas_ref[i] ~ lognormal(0, 2.5);
    eff_alpha[i] ~ normal(0, tau*lambdas[(J-1)+i]);
  }

  for (i in 1:N){
    theta[i] ~ normal(0 + X1[i]*eff_mu, exp(log(1) + X1[i]*eff_sigma));
  }
  
  for (i in 1:(J-1)){
    lambdas[i] ~ cauchy(0,1);
  }
  
  for (i in 1:N){
    y[i, 1] ~ bernoulli_logit(alpha_1*(theta[i] - beta_1));
    for (j in 1:(J-1)){
      y[i, j+1] ~ bernoulli_logit(exp(log(alphas_ref[j]) + X1[i]*eff_alpha[j])*(theta[i] - betas_ref[j] - eff_beta[j]*X1[i]));
    }
  }
}

