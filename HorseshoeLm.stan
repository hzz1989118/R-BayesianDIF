
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
  int<lower = 1> P;
  int<lower = 0, upper = 1> y[N]; // outputs
  matrix[N,P] X; // inputs
  real<lower=1> nu;
}

parameters {
  // intercept and noise std
  real w0;
  // real<lower=0> sigma;
  // auxiliary variables for the variance parameters
  vector[P] z;
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  vector<lower=0>[P] r1_local;
  vector<lower=0>[P] r2_local;
}

transformed parameters {

  // global and local variance parameters, and the input weights
  real<lower=0> tau;
  vector<lower=0>[P] lambda;
  vector[P] w;
  tau = r1_global * sqrt(r2_global);
  lambda = r1_local .* sqrt_vec(r2_local);
  w = z .* lambda*tau;
}

model {
  // observation model
  y ~ bernoulli_logit(w0 + X*w);
  // half t-priors for lambdas (nu = 1 corresponds to horseshoe)
  z ~ normal(0, 1);
  r1_local ~ normal(0.0, 1.0);
  r2_local ~ inv_gamma(0.5*nu, 0.5*nu);
  // half cauchy for tau
  r1_global ~ normal(0.0, 1.0);
  r2_global ~ inv_gamma(0.5, 0.5);
  // weakly informative prior for the intercept
  w0 ~ normal(0,5);
  // using uniform prior on the noise variance
}