N <- 1000
J <- 10
alpha <- runif(J, .5, 1.5)
beta <- runif(J, -1.5, 1.5)
theta <- rnorm(N, 0, 1)

generate2NOres <- function(theta, alpha, beta){
  core <- (-1)*alpha*outer(beta, theta, "-")
  out <- (runif(length(theta)*length(alpha)) < pnorm(core))*1  
  return(t(out))
}

Y <- generate2NOres(theta, alpha, beta)

init <- list("Z" = matrix(nrow = N, ncol = J),
             "theta" = theta + rnorm(N, 0, .5),
             "xi" = cbind(alpha + rnorm(J, 0.1, .01), alpha*beta + runif(J, -.05, .05)),
             "mu" = 0,
             "sigma" = 1)

#########################################################
runMCMC <- function(Y, init, nkeep, nburnin){
  curr <- init
  out <- NULL

  for (i in 1:nburnin) {
    curr$Z <- gb_sampler_Z(curr, Y) # augment variable
    curr$theta <- gb_sampler_theta(curr)
    curr$xi <- gb_sampler_xi(curr)
  }
  
  curr$Z <- gb_sampler_Z(curr, Y)
  curr$theta <- gb_sampler_theta(curr)
  curr$xi <- gb_sampler_xi(curr)
  
  out$theta <- curr$theta
  out$alpha <- curr$xi[, 1]
  out$gamma <- curr$xi[, 2]
  
  for(i in 2:nkeep){
    curr$Z - gb_sampler_Z(curr, Y)
    curr$theta <- gb_sampler_theta(curr)
    curr$xi <- gb_sampler_xi(curr)
    out$theta <- cbind(out$theta, curr$theta)
    out$alpha <- cbind(out$alpha, curr$xi[,1])
    out$gamma<- cbind(out$gamma, curr$xi[,2])
  }
  
  return(out)
}

library(msm)

gb_sampler_Z <- function(curr, Y){
  I <- dim(Y)[1]
  J <- dim(Y)[2]
  Z <- matrix(nrow = I, ncol = J)
  theta <- curr$theta
  alpha <- curr$xi[, 1]
  gamma <- curr$xi[, 2]
  
  for (i in 1 : I){
    for (j in 1: J){
      eta <- alpha[j]*theta[i] - gamma[j]
      if (Y[i, j] == 1){
        Z[i, j] <- rtnorm(1, eta, sd = 1, lower = 0)
      } else {
        Z[i, j] <- rtnorm(1, eta, sd = 1, upper = 0)
      }
    }
  }
  
  return(Z)
}

library(MASS)

gb_sampler_theta <- function(curr){
  I <- dim(Y)[1]
  J <- dim(Y)[2]
  mu <- curr$mu
  sigma <- curr$sigma
  Z <- curr$Z
  theta <- matrix(nrow = I, ncol = 1)
  alpha <- curr$xi[, 1]
  gamma <- curr$xi[, 2]
  
  for (i in 1 : I){
        theta[i] <- rnorm(n = 1, mean = (sum((Z[i,] + gamma)*alpha) + mu/(sigma^2))/(1/(sigma^2) + sum(alpha^2)), sd = sqrt(1/(1/(sigma^2) + sum(alpha^2))))
  }
  
  return(theta)
}

gb_sampler_xi <- function(curr){
  I <- dim(Y)[1]
  J <- dim(Y)[2]
  xi <- matrix(nrow = J, ncol = 2)
  alpha <- curr$xi[, 1]
  Z <- curr$Z
  theta <- curr$theta
  X <- cbind(theta, rep(-1, I))
  
  for (j in 1 : J){
    xi[j, ] <- mvrnorm(mu = solve(t(X)%*%X)%*%t(X)%*%Z[,j], Sigma = solve(t(X)%*%X))
    if (xi[j, 1] < 0){
      xi[j, 1] <- 0
    } 
  }
  
  return(xi)
}

###############################
system.time(gibbmcmc1 <- runMCMC(Y = Y, init = init, nkeep = 30, nburnin = 50)) 
