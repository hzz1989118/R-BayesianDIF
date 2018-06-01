N <- 1000
P <-  10

X <- NULL
for ( i in 1:P){
  temp <- rnorm(N, 0, 1)
  X <-  cbind(X, temp)
}

nu <- 1

betas <- runif(5, 1, 2)

# y <- apply(X[,1:5], 1, FUN = function(x) rnorm(n = 1, mean = x%*%betas, sd = .5))

y <- apply(X[,1:5], 1, FUN = function(x) rbinom(n = 1, size = 1, prob = 1/(1+exp(-x%*%betas))))

sim_dat <- list("N" = N, "P" = P,
                "y" = y, "X" = X, "nu" = nu)

fitlmSelect <- stan(file = "/Users/zhuangzhuanghan/Study/R/src/bayesianDIF/HorseshoeLm.stan",
                    data = sim_dat, chains = 1,
                    iter = 2000)