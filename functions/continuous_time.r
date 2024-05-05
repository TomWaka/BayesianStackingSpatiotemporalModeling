library(LaplacesDemon)
library(Matrix)
library(MCMCpack)
library(MASS)
library(statmod)  
library(mvtnorm)
library(quadprog)
library(parallel)
library(nloptr)

# ----------- continuous ---------------#

# Spatiotemporal correlation kernel
K_phi <- function(t1, t2, gamma1, gamma2, phi1, phi2) {
  1 / (phi1 * abs(t1 - t2)^2 + 1) * exp(-phi2 * sqrt(sum(gamma1 - gamma2)^2) / sqrt(1 + phi1 * abs(t1 - t2)^2))
}

# Temporal correlation kernel for beta
C_xi <- function(t1, t2, xi) {
  exp(- xi^2 * abs(t1 - t2))
}


# Function to simulate data
simulate_data <- function(n, p, sigma, phi1, phi2, xi, delta_beta, delta_z) {
  # Time and location points
  times <- sort(runif(n))
  locations <- matrix(0, n, 2)  # matrix(runif(n * 2), ncol = 2)
  locations[1,] <- runif(2)
  for (i in 2:n) {
    lll <- locations[i-1,] + runif(2, -1, 1) / 2
    lll[lll > 1] <- 2 - lll[lll > 1]
    lll[lll <0]  <- - lll[lll < 0]
    locations[i,] <- lll
  }
  
  # Generate covariance matrix C for beta using C_xi
  C_beta <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      C_beta[i, j] <- C_xi(times[i], times[j], xi)
    }
  }
  
  # Generate covariance matrix S for z using K_phi
  KK <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      KK[i, j] <- K_phi(times[i], times[j], locations[i,], locations[j,], phi1, phi2)
    }
  }
  
  # Generate z, beta, and eta
  z <- mvrnorm(1, rep(1, n), sigma^2 * delta_z^2 * KK)
  beta <- mvrnorm(p, rep(1, n), sigma^2 * delta_beta^2 * C_beta)
  eta <- rnorm(n, 0, sigma^2)
  
  # Generate X and Y
  x <- matrix(rnorm(n * p, 0, 4), ncol = p)
  trend <- numeric(n)  
  for (i in 1:n) {
    trend[i] <- sum(beta[, i] * x[i, ])
  }
  
  y <- trend + z + eta
  
  list(x = x, y = y, true = trend + z, z=z, beta=beta, times = times, locations = locations)
}

# Compute the posterior
compute_posterior <- function(y, x, times, locations, n, p, a_sigma, b_sigma, phi1, phi2, xi, delta_beta, delta_z) {
  
  # Generate covariance matrix C for beta using C_xi
  C_beta <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      C_beta[i, j] <- C_xi(times[i], times[j], xi)
    }
  }
  
  # Generate covariance matrix S for z using K_phi
  KK <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      KK[i, j] <- K_phi(times[i], times[j], locations[i,], locations[j,], phi1, phi2)
    }
  }
  
  
  CC <- kronecker(diag(p), C_beta) 
  S <- bdiag(diag(n), delta_beta^2 * CC, delta_z^2 * KK) 
  S <- S + diag(1e-10 , nrow(S))
  S_inv <- solve(S)
  
  xx <- matrix(0, n, n * p)
  for (i in 1:p) {
    xx[,((i-1)*n + 1):(i*n)] <- diag(x[, i])
  }
  
  X <- rbind(cbind(xx, diag(n)),
             cbind(diag(n*p), matrix(0, p*n, n)),
             cbind(matrix(0, n, p*n), diag(n))   )
  Y <- c(y, array(0, (p+1)*n) )
  
  SX <- solve(S, X)
  Sigma_inv <- t(X) %*% SX
  Sigma_inv <- (t(as.matrix(Sigma_inv ))+ as.matrix(Sigma_inv ))/2
  Sigma <- solve(Sigma_inv)
  m <- Sigma %*% t(SX) %*% Y
  
  a_sigma_star <- a_sigma + n / 2
  b_sigma_star <- b_sigma + t(Y - X %*% m) %*% S_inv %*% (Y - X %*% m)/2
  sigma2_post <- rinvgamma(1, shape = a_sigma_star, scale = as.matrix(b_sigma_star))
  
  #Sigma <- (t(as.matrix(Sigma)) + as.matrix(Sigma))/2 + diag(1e-3 , nrow(Sigma))
  Sigma <- (t(as.matrix(Sigma ))+ as.matrix(Sigma ))/2
  theta_post <- rmvnorm(1, as.matrix(m), sigma2_post * Sigma)
  
  list(Sigma = Sigma, theta = c(as.matrix(m)),
       a_sigma_star=a_sigma_star, b_sigma_star=as.vector(b_sigma_star),
       y_filt = (X %*% m)[1:n,], y_fil_cov = cbind(xx, diag(n)) %*% Sigma %*% t(cbind(xx, diag(n))),
       beta = c(as.matrix(m))[1:(p*n)], z = c(as.matrix(m))[(p*n+1):(n*(p+1))])
}

# predict z
predict_z0_distribution <- function(T, T0, G, G0, z_sampled, Sigma, sigma2, delta_z, phi1, phi2) {
  # T: Known space-time points (existing time series data)
  # T0: Unknown space-time points (time series data you want to predict)
  # z_sampled: Sampled z values at known space-time points
  # sigma2: Value of sigma^2
  # delta_z: Value of delta_z
  
  n <- length(T)
  n0 <- length(T0)  # This is the number of unknown space-time points
  
  # Initialize covariance matrices
  C_z <- matrix(0, n, n)
  C_z0 <- matrix(0, n, n0)
  C_z00 <- matrix(0, n0, n0)
  
  # Calculate C_z
  for (i in 1:n) {
    for (j in 1:n) {
      C_z[i, j] <- sigma2 * delta_z^2 * K_phi(T[i], T[j], G[i], G[j], phi1, phi2)
    }
  }
  
  # Calculate C_z0
  for (i in 1:n) {
    for (j in 1:n0) {
      C_z0[i, j] <- sigma2 * delta_z^2 * K_phi(T[i], T0[j], G[i], G0[j], phi1, phi2)
    }
  }
  
  # Calculate C_z00
  for (i in 1:n0) {
    for (j in 1:n0) {
      C_z00[i, j] <- sigma2 * delta_z^2 * K_phi(T0[i], T0[j], G0[i], G0[j], phi1, phi2)
    }
  }
  
  common_term <- solve(C_z + diag(1e-10, nrow(C_z)) , C_z0)
  
  # Calculate the predictive mean and variance for z_0
  predictive_mean_z0 <- as.vector(t(common_term) %*% z_sampled)
  predictive_var_z0 <- (C_z00 - t(common_term) %*% C_z0) + t(common_term) %*% Sigma[(p*n+1):((p+1)*n),(p*n+1):((p+1)*n)] %*% (common_term) 
  
  return(list(mean = predictive_mean_z0, variance = predictive_var_z0))
}

# predict beta 
predict_beta_distribution <- function(k, T, T0, beta_sampled, Sigma, sigma2, xi, delta_beta, phi1, phi2, p) {
  # T: Known space-time points (existing time series data)
  # T0: Unknown space-time points (time series data you want to predict)
  # beta_sampled: Sampled beta values at known space-time points
  # sigma2, delta_z, phi1, phi2, p: Model parameters
  
  n <- length(T)
  n0 <- length(T0)
  
  beta <- beta_sampled[(n*(k-1)+1):(n*k)]
  
  # Initialize covariance matrices
  C_beta <- matrix(0, n, n)
  C_beta0 <- matrix(0, n, n0)
  C_beta00 <- matrix(0, n0, n0)
  
  # Calculate C_beta
  for (i in 1:n) {
    for (j in 1:n) {
      C_beta[i, j] <- sigma2 * delta_beta^2 * C_xi(T[i], T[j], xi)
    }
  }
  
  # Calculate C_beta0
  for (i in 1:n) {
    for (j in 1:n0) {
      C_beta0[i, j] <- sigma2 * delta_beta^2 * C_xi(T[i], T0[j], xi)
    }
  }
  
  # Calculate C_beta00
  for (i in 1:n0) {
    for (j in 1:n0) {
      C_beta00[i, j] <- sigma2 * delta_beta^2 * C_xi(T0[i], T0[j], xi)
    }
  }
  
  common_term <- solve(C_beta, C_beta0)
  
  predictive_mean_beta <- as.vector(t(common_term) %*% beta)
  predictive_var_beta <- C_beta00 - t(common_term) %*% C_beta0 + t(common_term) %*% Sigma[(n*(k-1)+1):(n*k),(n*(k-1)+1):(n*k)] %*% (common_term) 
  
  return(list(mean = predictive_mean_beta, variance = predictive_var_beta))
}

pred_dist_continuous <- function(y, x, x0, T, T0, G, G0, n, p, a_sigma, b_sigma, phi1, phi2, xi, delta_beta, delta_z) {
  
  n0 <- length(T0)
  posterior = compute_posterior(y, x, T, G, n, p, a_sigma, b_sigma, phi1, phi2, xi, delta_beta, delta_z)
  
  result <- predict_z0_distribution(T, T0, G, G0, posterior$theta[(p*n+1):((p+1)*n)], posterior$Sigma, posterior$b_sigma_star/posterior$a_sigma_star, delta_z, phi1, phi2)
  predictive_mean_z0 <- result$mean
  predictive_var_z0 <- result$variance
  
  predictive_mean_beta0 <- array(0, dim=c(p,n0))
  predictive_var_beta0  <- array(0, dim=c(p,n0,n0))
  for (k in 1:p) {
    result <- predict_beta_distribution(k, T, T0, posterior$theta[1:(p*n)], posterior$Sigma, posterior$b_sigma_star/posterior$a_sigma_star, xi, delta_beta, phi1, phi2, p)
    predictive_mean_beta0[k,] <- result$mean
    predictive_var_beta0[k,,]  <- result$variance
  }
  
  trend0 <- numeric(n0)  
  for (i in 1:n0) {
    trend0[i] <- sum(predictive_mean_beta0[, i] * x0[i, ])
  }
  predictive_mean_y0 <- trend0 + predictive_mean_z0
  
  xx0 <- matrix(0, n0, n0 * p)
  beta0_var <- matrix(0, n0 * p, n0 * p)
  for (i in 1:p) {
    xx0[,((i-1)*n0+1):(i*n0)] <- diag(x0[, i])
    beta0_var[(n0*(i-1)+1):(n0*i),(n0*(i-1)+1):(n0*i)] <- predictive_var_beta0[i,,]
  }
  
  predictive_var_y0 <- diag(posterior$b_sigma_star/posterior$a_sigma_star, n0) + 
    predictive_var_z0 + xx0 %*% beta0_var %*% t(xx0)
  
  return(list(mean = predictive_mean_y0, variance = predictive_var_y0, df = a_sigma + n/2, z = c(predictive_mean_z0), z_sampled = posterior$theta[(p*n+1):((p+1)*n)]  ))
}




compute_DIC <- function(y, bp_sample, bp_sigma, result_point, T) {
  
  stack_sample <- c()
  for (i in 1:nrow(param_grid)) {
    if (round(result_point$solution*1000)[i] > 0) {
      if (round(result_point$solution*1000)[i] > 1000){
        stack_sample <- 
          rbind(stack_sample, bp_sample[1:1000, ,i])
      } else {
        stack_sample <- 
          rbind(stack_sample, bp_sample[(1:round(result_point$solution*1000)[i]), ,i])
      }
    }
  }    
  sigma_sample <- c()
  for (i in 1:nrow(param_grid)) {
    if (round(result_point$solution*1000)[i] > 0) {
      if (round(result_point$solution*1000)[i] > 1000){
        sigma_sample <- 
          c(sigma_sample, bp_sigma[1:1000, ,i])
      } else {
        sigma_sample <- 
          c(sigma_sample, bp_sigma[(1:round(result_point$solution*1000)[i]), ,i])
      }
    }
  }
  
  stack_mean <- colMeans(stack_sample)
  sigma_mean <- mean(sigma_sample)
  D_bar <- 0
  for (i in 1:length(sigma_sample)) {
    D_bar <- D_bar - 2 * dmvn(y, stack_sample[i,], sigma_sample[i] * diag(T), log = TRUE) /length(sigma_sample)
  }
  D_theta_bar <- - 2 * dmvn(y, stack_mean, sigma_mean * diag(T), log = TRUE)
  p_D <- D_bar - D_theta_bar
  dic <- p_D + D_bar
  
  return(dic)
}

compute_WAIC <- function(y, bp_sample, bp_sigma, result_point, T) {
  stack_sample <- c()
  for (i in 1:nrow(param_grid)) {
    if (round(result_point$solution*1000)[i] > 0) {
      if (round(result_point$solution*1000)[i] > 1000){
        stack_sample <- 
          rbind(stack_sample, bp_sample[1:1000, ,i])
      } else {
        stack_sample <- 
          rbind(stack_sample, bp_sample[(1:round(result_point$solution*1000)[i]), ,i])
      }
    }
  }    
  sigma_sample <- c()
  for (i in 1:nrow(param_grid)) {
    if (round(result_point$solution*1000)[i] > 0) {
      if (round(result_point$solution*1000)[i] > 1000){
        sigma_sample <- 
          c(sigma_sample, bp_sigma[1:1000, ,i])
      } else {
        sigma_sample <- 
          c(sigma_sample, bp_sigma[(1:round(result_point$solution*1000)[i]), ,i])
      }
    }
  }
  
  x <- matrix(0, nrow = 1, ncol = length(sigma_sample))
  for (i in 1:length(sigma_sample)) {
    x[,i] <- dmvn(y, stack_sample[i,], sigma_sample[i] * diag(T), log = TRUE) 
  }
  lppd <- sum(log(rowMeans(exp(x))))
  pWAIC1 <- 2 * sum(log(rowMeans(exp(x))) - rowMeans(x))
  pWAIC2 <- sum(.rowVars(x))
  WAIC <- -2 * lppd + 2 * pWAIC2
  return(list(WAIC = WAIC, lppd = lppd, pWAIC = pWAIC2, pWAIC1 = pWAIC1))
}
