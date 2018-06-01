######## Data Generation ##########
set.seed(1108)
N <- 3000
J <- 10
X <- sample(c(0,1), N, replace = T)
eff_mu <- -.5
eff_sigma <- .5
eff_alpha <- c(0, 1.5, 0, rep(0, 7))
eff_beta <- c(0, .8, .8, rep(0, 7))
# XindEff <- 0
# XdEff <- rep(0, 10)
theta <- NULL
for (i in 1:N){
  theta[i] <- rnorm(n = 1, mean = X[i]*eff_mu, sd = exp(log(1) + X[i]*eff_sigma))
}
beta <- rnorm(J, .5, 0.3)
alpha <- runif(J, 0.5, 1.5)

#### Functions are Needed ###
getIRTres <- function(theta, a, b, eff_a, eff_b, X){
  prob <- (1 + exp(exp((log(a) + outer(eff_a, X)))*(outer(b, theta, "-") + outer(eff_b, X))))^(-1)
  threshold <- runif(length(a)*length(theta), 0, 1)
  out <- (prob > threshold)*1
  return(t(out))
}

y <- getIRTres(theta = theta, a = alpha, b = beta, eff_a = eff_alpha, eff_b = eff_beta, X = X)


sim_dat <- list("N" = N, "J" = J,
                "y" = y, "X1" = X)

stan_model(file = "/Users/hzz1989118/Documents/Study/R/src/bayesianDIF/plainBDIF.stan")
######## Fit #######

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

fit1 <- stan(file = "/Users/zhuangzhuanghan/Study/R/src/bayesianDIF/Horseshoe.stan", 
             data = sim_dat, chains = 1, thin = 5,
             iter = 5000)

fit2 <- stan(file = "/Users/hzz1989118/Documents/Study/R/src/bayesianDIF/plainBDIF.stan", 
             data = sim_dat, chains = 1, thin = 5,
             iter = 5000)

fit3 <- stan(file = "/Users/zhuangzhuanghan/Study/R/src/bayesianDIF/BGroupLasso.stan", 
             data = sim_dat, chains = 1, thin = 5,
             iter = 5000)

fit4 <- stan(file = "/Users/zhuangzhuanghan/Study/R/src/bayesianDIF/BLasso.stan", 
             data = sim_dat, chains = 1, thin = 5,
             iter = 5000)