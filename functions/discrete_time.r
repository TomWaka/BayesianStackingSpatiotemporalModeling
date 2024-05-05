library(LaplacesDemon)       
library(Matrix)
library(MCMCpack)
library(MASS)
library(statmod)  
library(mvtnorm)
library(quadprog)
library(parallel)
library(nloptr)


# Matérn kernel 
matern_kernel <- function(s1, s2, phi, nu) {
  # Modified Bessel function of the second kind
  h <- sqrt(sum((s1 - s2)^2))
  if (h==0) {
    return(1)
  } else {
    k <- besselK( h/phi, nu)
    # Matérn covariance
    cov_value <- (2^(1-nu)/gamma(nu)) * (( h/phi)^nu) * k
    return(cov_value)
  }
}


# Function to simulate data
simulate_data <- function(T, n, p, sigma, phi, nu, delta_beta, delta_z, W, alpha) {
  
  # Time and location points
  times <- 1:T
  locations <- matrix(0, n, 2)  # matrix(runif(n * 2), ncol = 2)
  locations[1,] <- runif(2)
  for (i in 2:n) {
    lll <- locations[i-1,] + runif(2, -1, 1) / 2
    lll[lll > 1] <- 2 - lll[lll > 1]
    lll[lll <0]  <- - lll[lll < 0]
    locations[i,] <- lll
  }
  B <- matrix(0, T, n*T)
  for (t in 1:T) {
    j <- t %% n + 1   #sample(c(1:n), siz=1)  # 
    B[t,j+(t-1)*n] <- 1
  }
  
  # Construct kernel matrix
  K_phi <- matrix(0, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      K_phi[i, j] <- matern_kernel(locations[i, ], locations[j, ], phi, nu)
    }
  }
  
  # Generate z, beta, and eta
  z <- matrix(0, T*n)
  z[1:n] <- mvrnorm(1, rep(0, n), sigma^2 * delta_z^2 * K_phi)
  for (t in 1:(T-1)) {
    z[((n*(t)+1):(n*(t+1)))] <- alpha * z[((n*(t-1)+1):(n*t))] + mvrnorm(1, rep(0, n), sigma^2 * delta_z^2 * K_phi)
  }
  w = as.vector(B %*% z)
  
  beta <- matrix(0, T*p)
  beta[1:p] <- mvrnorm(1, rep(0, p), sigma^2 * delta_beta^2 * W)
  for (t in 1:(T-1)) {
    beta[((p*(t)+1):(p*(t+1)))] <- alpha * beta[((p*(t-1)+1):(p*(t)))] + mvrnorm(1, rep(0, p), sigma^2 * delta_beta^2 * W)
  }
  
  eta <- rnorm(T, 0, sigma^2)
  
  # Generate X and Y
  x <- matrix(0, T, p*T)
  for (t in 1:T) {
    x[t,(((t-1)*p+1):(t*p))] <- rnorm(p, 0, 4)
  } 
  trend <- x %*% beta
  
  y <- trend + w + eta
  
  list(x = x, y = y, z=z, B=B, w=w, beta=beta, times = times, locations = locations, true = trend + w)
}


# Compute the posterior
compute_posterior <- function(y, x, x0, B, T, T0, W, locations, p, a_sigma, b_sigma, phi, nu, delta_beta, delta_z, alpha) {
  
  #locations_ <- unique(locations)
  n <- nrow(locations)
  #n <- ncol(B) / T
  
  # Generate covariance matrix S for z using K_phi
  KK <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      KK[i, j] <- matern_kernel(locations[i, ], locations[j, ], phi, nu)
    }
  }
  
  zero_vector <- matrix(0, nrow = T-1, ncol = 1)  
  zero_vector_transpose <- t(zero_vector)  
  identity_matrix <- diag(T-1)  
  single_zero <- matrix(0, nrow = 1, ncol = 1)  
  A <- rbind(
    cbind(zero_vector_transpose, single_zero),
    cbind(alpha * identity_matrix, zero_vector)
  )
  SnA  <- solve(diag(T) - A)
  SIP  <- Matrix(kronecker(SnA, diag(p)))
  SIN  <- Matrix(kronecker(SnA, diag(n)))
  
  IP <- Matrix(kronecker(diag(T) - A, diag(p)))
  IN <- Matrix(kronecker(diag(T) - A, diag(n)))
  
  #S <- bdiag(diag(T), delta_beta^2 * SIP %*% kronecker(diag(T), W) %*% t(SIP), delta_z^2 * SIN %*% kronecker(diag(T), KK) %*% t(SIN))
  #solve(delta_beta^2 * SIP %*% kronecker(diag(T), W) %*% t(SIP))
  #solve(delta_z^2 * SIN %*% kronecker(diag(T), KK) %*% t(SIN))
  KTW <- Matrix(kronecker(diag(T), solve(W)))
  KTK <- Matrix(kronecker(diag(T), solve(KK)))
  S_inv <- Matrix(bdiag(diag(T), delta_beta^{-2} * t(IP) %*% KTW %*% (IP), delta_z^{-2} * t(IN) %*% KTK %*% (IN)))
  
  
  X <- Matrix(rbind(cbind(x, B),
                    cbind(diag(T*p), matrix(0, T*p, ncol(B) )),
                    cbind(matrix(0, n*T, p*T), diag(n*T))
  ))
  Y <- c(y, array(0, (p+n)*T ))
  
  SX <- S_inv %*% X
  Sigma_inv <- t(X) %*% SX
  Sigma <- as.matrix(solve(Sigma_inv))
  m <- Sigma %*% t(SX) %*% Y
  a_sigma_star <- a_sigma + T / 2
  b_sigma_star <- b_sigma + t(Y - X %*% m) %*% S_inv %*% (Y - X %*% m)/2
  sigma2_post  <- rinvgamma(1, shape = a_sigma_star, scale = as.matrix(b_sigma_star))
  
  Sigma <- as.matrix(Sigma)
  
  list(Sigma = Sigma,  theta = c(as.matrix(m)), filt = c((X %*% m)[1:T,]), fil_cov = cbind(x, B) %*% Sigma %*% t(cbind(x, B)), fil_z = c(B %*% m[(p*T+1):(n*T+p*T)]), fil_b = c(m[1:(p*T)]), a_sigma_star=a_sigma_star, b_sigma_star=as.vector(b_sigma_star), df = T / 2)
}




# predict z
predict_z0_distribution <- function(y, x, x0, B, T, T0, locations, new_locations, theta_post, Sigma, sigma2, delta_z, phi, nu) {
  # T: Known space-time points (existing time series data)
  # T0: Unknown space-time points (time series data you want to predict)
  # z_sampled: Sampled z values at known space-time points
  # sigma2: Value of sigma^2
  # delta_z: Value of delta_z
  n <- nrow(locations)
  n0 <- nrow(new_locations)
  
  C_z <- matrix(0, n, n)
  C_z0 <- matrix(0, n, n0)
  C_z00 <- matrix(0, n0, n0)
  
  # Calculate C_z
  for (i in 1:n) {
    for (j in 1:n) {
      C_z[i, j] <- sigma2 * delta_z^2 * matern_kernel(locations[i,], locations[j,], phi, nu)
    }
  }
  
  # Calculate C_z0
  for (i in 1:n) {
    for (j in 1:n0) {
      C_z0[i, j] <- sigma2 * delta_z^2 * matern_kernel(locations[i,], new_locations[j,], phi, nu)
    }
  }
  
  # Calculate C_z00
  for (i in 1:n0) {
    for (j in 1:n0) {
      C_z00[i, j] <- sigma2 * delta_z^2 * matern_kernel(new_locations[i,], new_locations[j,], phi, nu)
    }
  }
  
  common_term <- solve(C_z, C_z0)
  
  z <- theta_post[(length(theta_post)-n+1):length(theta_post)]
  predictive_mean_z0 <- t(common_term) %*% matrix(z)
  predictive_var_z0  <- C_z00 - t(common_term) %*% C_z0 + t(common_term) %*% Sigma[(n*T+p*T-n+1):(n*T+p*T),(n*T+p*T-n+1):(n*T+p*T)] %*% (common_term)     
  
  return(list(mean = predictive_mean_z0, variance = predictive_var_z0))
}


# predict beta 
predict_beta_distribution <- function(y, x, x0, B, T, T0, theta_post, W, Sigma, sigma2, delta_beta, phi, nu, p) {
  # T: Known space-time points (existing time series data)
  # T0: Unknown space-time points (time series data you want to predict)
  # beta_sampled: Sampled beta values at known space-time points
  # sigma2, delta_z, phi1, phi2, p: Model parameters
  
  beta <- theta_post[(p*(T-1)+1):(T*p)]
  predictive_mean_beta <- beta
  predictive_var_beta <- sigma2 * delta_beta^2 * W + sigma2 * Sigma[((T-1)*p+1):(T*p),((T-1)*p+1):(T*p)]
  
  return(list(mean = predictive_mean_beta, variance = predictive_var_beta))
}


pred_dist_dis <- function(y, x, x0, B, T, T0, locations, new_locations, W, p, a_sigma, b_sigma, phi, nu, delta_beta, delta_z) {
  
  new_locations <- matrix(new_locations, ncol=2)
  posterior = compute_posterior(y, x, x0, B, T, T0, W, locations, p, a_sigma, b_sigma, phi, nu, delta_beta, delta_z)
  
  result <- predict_z0_distribution(y, x, x0, B, T, T0, locations, new_locations, posterior$theta, posterior$Sigma, posterior$b_sigma_star/posterior$a_sigma_star, delta_z, phi, nu)
  predictive_mean_z0 <- result$mean
  predictive_var_z0 <- result$variance
  
  result <- predict_beta_distribution(y, x, x0, B, T, T0, posterior$theta, W, posterior$Sigma, posterior$b_sigma_star/posterior$a_sigma_star, delta_beta, phi, nu, p)
  predictive_mean_beta0 <- result$mean
  predictive_var_beta0  <- result$variance
  
  beta0 <- matrix(0, T0*p)
  for (t in 1:T0) {
    beta0[(p*(t-1)+1):(p*t)] <- predictive_mean_beta0 
  }
  
  trend0 <- x0 %*% beta0 
  predictive_mean_y0 <- trend0 + predictive_mean_z0
  predictive_var_y0 <- posterior$b_sigma_star/posterior$a_sigma_star + 
    predictive_var_z0 + (x0) %*% kronecker(diag(T0),predictive_var_beta0) %*% t(x0)
  
  return(list(mean = predictive_mean_y0, variance = predictive_var_y0, df = a_sigma + T/2, z = c(predictive_mean_z0), z_sampled = posterior$theta[(p*n+1):((p+1)*n)]  ))
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

