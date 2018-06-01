data{
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 0, upper = 1> y[N, J];
  int<lower = 0, upper = 1> X1[N];
}

parameters{
  vector[N] theta;
  vector<lower = 0>[J] alpha; //discrimination
  vector[J] betas; // item difficulty
  real X1_indirect_eff; // indirect effect from X1 to items through latent trait
  vector[J-1] X1_direct_eff; // direct effect from X1 to items
}

transformed parameters{
  vector[J] dEff;
  for(i in 1:J){
    if (i == J){
      dEff[i] = 0 - sum(X1_direct_eff);
    } else {
      dEff[i] = X1_direct_eff[i];
    }
  }
}

model{
  for (i in 1:J){
    betas[i] ~ normal(0, 2);
    alpha[i] ~ lognormal(0, 1);
    // X1_direct_eff[i] ~ normal(0, 25);
  }
  
  X1_direct_eff ~ normal(0, 1);
  
 /* X1_indirect_eff ~ normal(0, .01);*/
  
  for (i in 1:N){
    theta[i] ~ normal(X1[i]*X1_indirect_eff, 1);
  }
  
  for (i in 1:N){
    for (j in 1:J){
      y[i, j] ~ bernoulli(inv_logit(alpha[j]*(theta[i] - betas[j] - dEff[j]*X1[i])));
    }
  }
}

