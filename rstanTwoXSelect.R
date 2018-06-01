N <- 1000
J <- 10
X1 <- sample(c(0,1), N, replace = T)
X2 <- sample(c(0,1), N, replace = T)

X1indEff <- -.5
X1dEff <- rep(0, 10)

X2indEff <- 0
X2dEff <- c(rep(0.8, 8), .5, 0)

theta <- NULL
for (i in 1:N){
  theta[i] <- rnorm(n = 1, mean = X1[i]*X1indEff + X2[i]*X2indEff, sd = 1)
}
beta <- rnorm(J, .5, 0.3)
alpha <- runif(J, 0.5, 1.5)

Prob <- (1 + exp(alpha*(outer(beta, theta, "-") + outer(X1dEff, X1) + outer(X2dEff, X2))))^(-1)
threshold <- runif(length(alpha)*length(theta), 0, 1)

y <- t((Prob > threshold)*1)

sim_dat <- list("N" = N, "J" = J,
                "y" = y, "X1" = X1, "X2" = X2)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

fitTwo <- stan(file = "/Users/zhuangzhuanghan/Study/R/src/bayesianDIF/HorseshoeTwoX.stan", 
             data = sim_dat, chains = 1,
             iter = 2000)