set.seed(1108)
N <- 3000
J <- 10
Xs <- sample(c(0,1), N, replace = T)
a <- runif(J, .5, 1.5)
b <- runif(J, -1.5, 1.5)
a_eff <- c(0, .3, rep(0, times = J-2))
b_eff <- c(0, 0, 0.5, -0.5, rep(0, J-4))
th <- rnorm(N, 0 + Xs*.5, 1*exp(Xs*log(0.75)))

########################################################
generate2NOres <- function(theta, alpha, beta, alpha_eff, beta_eff, X){
  I <- length(theta)
  J <- length(alpha)
  out <- matrix(nrow = I, ncol = J)
  for(i in 1 :I){
    core <- (alpha + X[i]*alpha_eff)*(theta[i] -(beta + X[i]*beta_eff))
    out[i, ] <- (runif(J) < pnorm(core))*1  
  }
  return(out)
}

Y <- generate2NOres(th, a, b, a_eff, b_eff, Xs)

init <- list(# "Z" = matrix(nrow = N, ncol = J),
             "Z" = gb_sampler_Z(init, Y, X = Xs),
             "theta" = th + rnorm(N, 0, .5),
             "xi" = cbind(a + rnorm(J, 0.1, .01), a*b + runif(J, -.05, .05)),
             "mu_eff" = 0,
             "sigma_eff" = 0,
             "dif" = matrix(0.1, nrow = J-1, ncol = 2),
             "lambda2" = 1,
             "inv_tau2" = rep(1, J-1))

#########################################################
runMCMC <- function(Y, X, init, nkeep, nburnin){
  curr <- init
  out <- NULL
  
  for (i in 1:nburnin) {
    curr$Z <- gb_sampler_Z(curr, Y, X) # augment variable
    curr$theta <- gb_sampler_theta(curr, X)
    curr$xi <- gb_sampler_xi(curr, X)
    curr$dif <- gb_sampler_dif(curr , X)
    curr$inv_tau2 <- gb_sampler_inv_tau2(curr)
    curr$lambda2 <- gb_sampler_lambda2(curr)
    curr$mu_eff <- gb_sampler_mu_eff(curr, X)
    curr$sigma_eff <- gb_sampler_sigma_eff(curr, X)
  }
  
  curr$Z <- gb_sampler_Z(curr, Y, X)
  curr$theta <- gb_sampler_theta(curr, X)
  curr$xi <- gb_sampler_xi(curr, X)
  curr$dif <- gb_sampler_dif(curr , X)
  curr$inv_tau2 <- gb_sampler_inv_tau2(curr)
  curr$lambda2 <- gb_sampler_lambda2(curr)
  curr$mu_eff <- gb_sampler_mu_eff(curr, X)
  curr$sigma_eff <- gb_sampler_sigma_eff(curr, X)
  
  out$theta <- curr$theta
  out$alpha <- curr$xi[, 1]
  out$gamma <- curr$xi[, 2]
  out$alpha_dif <- curr$dif[, 1]
  out$gamma_dif <- curr$dif[, 2]
  out$inv_tau2 <- curr$inv_tau2
  out$lambda2 <- curr$lambda2
  out$mu_eff <- curr$mu_eff
  out$sigma_eff <- curr$sigma_eff
  
  for(i in 2:nkeep){
    curr$Z <- gb_sampler_Z(curr, Y, X)
    curr$theta <- gb_sampler_theta(curr, X)
    curr$xi <- gb_sampler_xi(curr, X)
    curr$dif <- gb_sampler_dif(curr , X)
    curr$inv_tau2 <- gb_sampler_inv_tau2(curr)
    curr$lambda2 <- gb_sampler_lambda2(curr)
    curr$mu_eff <- gb_sampler_mu_eff(curr, X)
    curr$sigma_eff <- gb_sampler_sigma_eff(curr, X)
    
    out$theta <- cbind(out$theta, curr$theta)
    out$alpha <- cbind(out$alpha, curr$xi[, 1])
    out$gamma <- cbind(out$gamma, curr$xi[, 2])
    out$alpha_dif <- cbind(out$alpha_dif, curr$dif[, 1])
    out$gamma_dif <- cbind(out$gamma_dif, curr$dif[, 2])
    out$inv_tau2 <- cbind(out$inv_tau2, curr$inv_tau2)
    out$lambda2 <- cbind(out$lambda2, curr$lambda2)
    out$mu_eff <- cbind(out$mu_eff, curr$mu_eff)
    out$sigma_eff <- cbind(out$sigma_eff, curr$sigma_eff)
  }
  return(out)
}

############################################

library(msm) # truncated normal

gb_sampler_Z <- function(curr, Y, X){
  I <- dim(curr$Z)[1]
  J <- dim(curr$Z)[2]
  Z <- matrix(nrow = I, ncol = J)
  theta <- curr$theta
  alpha <- curr$xi[, 1]
  gamma <- curr$xi[, 2]
  dif_alpha <- curr$dif[, 1]
  dif_gamma <- curr$dif[, 2]
  mu_eff <- curr$mu_eff
  sigma_eff <- curr$sigma_eff
  
  for (i in 1 : I){
    eta <- alpha[1]*theta[i] - gamma[1]
    if (Y[i, 1] == 1){
      Z[i, 1] <- rtnorm(1, eta, sd = 1, lower = 0)
    } else {
      Z[i, 1] <- rtnorm(1, eta, sd = 1, upper = 0)
    }
    for (j in 2: J){
      eta <- (alpha[j] + X[i]*dif_alpha[j-1])*theta[i] - (gamma[j]+ X[i]*dif_gamma[j-1])
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

gb_sampler_theta <- function(curr, X){
  I <- dim(curr$Z)[1]
  J <- dim(curr$Z)[2]
  mu_eff <- curr$mu_eff
  sigma_eff <- curr$sigma_eff
  Z <- curr$Z
  theta <- matrix(nrow = I, ncol = 1)
  alpha <- curr$xi[, 1]
  gamma <- curr$xi[, 2]
  dif_alpha <- curr$dif[, 1]
  dif_gamma <- curr$dif[, 2]
  
  for (i in 1 : I){
    repeat{
      theta[i] <- rnorm(n = 1, mean = (sum((Z[i,] + gamma + X[i]*c(0, dif_gamma))*(alpha + X[i]*c(0, dif_alpha))) + (0+X[i]*mu_eff)/((1*exp(X[i]*sigma_eff))^2))/(1/((1*exp(X[i]*sigma_eff))^2) + sum((alpha + X[i]*c(0, dif_alpha))^2)), 
                      sd = sqrt(1/(1/((1*exp(X[i]*sigma_eff))^2) + sum((alpha + X[i]*c(0, dif_alpha))^2))))
      if (theta[i] <= 4 && theta >= -4){
        break
      }
    }
  }
  
  return(theta)
}

gb_sampler_xi <- function(curr, X){
  I <- dim(curr$Z)[1]
  J <- dim(curr$Z)[2]
  xi <- matrix(nrow = J, ncol = 2)
  ind <- which(X == 0)
  Z <- curr$Z[ind,] # indicate out augment realizations in reference group
  theta <- curr$theta
  X_1 <- cbind(theta, rep(-1, I))[ind, ]
  
  for (j in 1 : J){
    repeat{
      xi[j, ] <- mvrnorm(mu = solve(t(X_1)%*%X_1)%*%t(X_1)%*%Z[,j], Sigma = solve(t(X_1)%*%X_1))
      if(xi[j, 1] > 0){
        break
      }
    }
  }
  return(xi)
}

gb_sampler_dif <- function(curr, X){
  I <- dim(curr$Z)[1]
  J <- dim(curr$Z)[2]
  dif <- matrix(nrow = J-1, ncol = 2)
  xi <- curr$xi
  ind <- which(X == 1)
  Z <- curr$Z[ind,] # indicate out augment realizations in focal group
  theta <- curr$theta
  X_1 <- cbind(theta, rep(-1, I))[ind, ]
  inv_tau2 <- curr$inv_tau2
  
  for (j in 1 : (J-1)){
    repeat{
      dif[j, ] <- mvrnorm(mu = solve(t(X_1)%*%X_1 + diag(inv_tau2[j], 2, 2))%*%t(X_1)%*%(Z[,(j+1)] - X_1%*%xi[j+1, ]), Sigma = solve(t(X_1)%*%X_1+ diag(inv_tau2[j],2,2)))
      if(xi[j, 1] + dif[j, 1] > 0){
        break
      }
    }
  }
  return(dif)
}

library(statmod)

gb_sampler_inv_tau2 <- function(curr){
  J <- dim(curr$Z)[2]
  lambda2 <- curr$lambda2
  dif <- curr$dif
  inv_tau2 <- rep(NA, J-1)
  for (j in 1:(J-1)){
    inv_tau2[j] <- rinvgauss(n = 1, mean = sqrt(lambda2/(t(dif[j, ])%*%dif[j, ])), shape = lambda2)
  }
  return(inv_tau2)
}

gb_sampler_lambda2 <- function(curr){
  J <- dim(curr$Z)[2]
  p <- dim(curr$dif)[1]*dim(curr$dif)[2]
  k <- J-1
  inv_tau2 <- curr$inv_tau2
  lambda2 <- rgamma(n = 1, shape = (p+k)/2 + 1, rate = (1/2)*sum(1/inv_tau2) + 1.78)
  return(lambda2)
}

gb_sampler_mu_eff <- function(curr, X){
  ind <- which(X == 1)
  theta <- curr$theta[ind] 
  mu_eff <- mean(theta)
  return(mu_eff)
}

gb_sampler_sigma_eff <- function(curr, X){
  ind <- which(X == 1)
  theta <- curr$theta[ind] 
  sigma_eff <- sd(theta)
  return(log(sigma_eff))
}

###############################
system.time(gibbs_mcmc2 <- runMCMC(Y = Y, X = Xs, init = init, 5000, 5000))

system.time(gibbs_mcmc2_1 <- runMCMC(Y = Y, X = Xs, init = init, nkeep = 2000, 3000))

system.time(gibbs_mcmc2_2 <- runMCMC(Y = Y, X = Xs, init = init, nkeep = 5000, 10000))


