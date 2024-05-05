# import functions
source("functions/discrete_time.r")

# setting
T_values <- c(50, 70) + 1
p <- 2

# Discrete time
results_dis <- list()
for (current_T in T_values) {
  T <- current_T
  n <- T - 1

  # true parameters
  sigma <- 1
  phi <- 1 / 7
  nu <- 1
  delta_beta <- 1
  delta_z <- 1
  W <- diag(p)
  alpha <- 0.9

  data <- simulate_data(T, n, p, sigma, phi, nu, delta_beta, delta_z, W, alpha)

  test_id <- T
  train_id <- setdiff(c(1:T), test_id)
  y <- data$y[1:(test_id - 1)]
  y0 <- data$y[test_id]
  true <- data$true[1:(test_id - 1)]
  true0 <- data$true[test_id]
  x <- data$x[1:(test_id - 1), (1:((test_id - 1) * 2))]
  x0 <- matrix(data$x[test_id, (((test_id - 1) * 2) + 1):(test_id * 2)], nrow = 1)
  w <- data$w[1:((test_id - 1))]
  w0 <- data$w[test_id]
  z <- data$z[1:((test_id - 1) * n)]
  z0 <- data$z[((test_id - 1) * n + 1):(test_id * n)]
  beta <- data$beta[1:((test_id - 1) * p), ]
  beta0 <- data$beta[((test_id - 1) * p + 1):(test_id * p), ]
  T <- T - 1
  T0 <- 1
  B <- data$B[-test_id, ]
  B0 <- data$B[test_id, ]
  locations <- array(0, dim = c(T, 2))
  new_locations <- array(0, dim = c(T0, 2))
  loc_id <- array(0, T)
  for (t in 1:T) {
    locations[t, ] <- data$locations[which(B[t, ((1 + (t - 1) * n):(t * n))] == 1), ]
    loc_id[t] <- which(B[t, ((1 + (t - 1) * n):(t * n))] == 1)
  }
  locations <- unique(locations)
  new_locations[1, ] <- locations[which(B0 == 1) %% n, ]

  n_new <- length(unique(loc_id))
  B_new <- matrix(0, T, n_new * T)
  for (t in 1:T) {
    j <- loc_id[t]
    B_new[t, j + (t - 1) * n_new] <- 1
  }

  # Define k-fold cross-validation
  K <- 20
  folds <- lapply(1:K, function(k) c(1:(k * (T - 10) / K + 9)))

  a_sigma <- 1 / 2
  b_sigma <- 1 / 2
  # Parameter candidates
  phi_values <- c(1, 1 / 10)
  nu_values <- c(3, 1 / 3)
  delta_beta_values <- c(5, 1 / 5)
  delta_z_values <- c(5, 1 / 5)

  # Define parameter grid
  param_grid <- expand.grid(
    delta_beta = c(delta_beta_values),
    delta_z = c(delta_z_values),
    phi = c(phi_values),
    nu = c(nu_values)
  )

  # Placeholder for Bayes predictors
  bayes_preds <- array(0, dim = c(10, K, nrow(param_grid)))
  log_probs <- array(0, dim = c(K, nrow(param_grid)))
  y_valid <- array(0, dim = c(10, K))

  # Loop over each fold
  for (f in 1:K) {
    # Split data
    {
      train_index <- folds[[f]]
      valid_index <- max(train_index) + 1

      T_train_cv <- length(train_index)
      T_valid_cv <- 1

      x_train_cv <- x[train_index, (1:((length(train_index)) * 2))]
      x_valid_cv <- matrix(x[valid_index, (((length(train_index)) * 2) + 1):((length(train_index) + 1) * 2)], ncol = 2)

      y_train_cv <- matrix(y[train_index], nrow = 1)
      y_valid_cv <- matrix(y[(T_train_cv - 9):T_train_cv])
      y_valid[, f] <- y_valid_cv

      true_train_cv <- matrix(true[train_index], nrow = 1)
      true_valid_cv <- matrix(true[valid_index])

      locations_train_cv <- array(0, dim = c(T_train_cv, 2))
      locations_valid_cv <- array(0, dim = c(T_valid_cv, 2))
      loc_id <- array(0, T_train_cv)
      for (t in 1:T_train_cv) {
        locations_train_cv[t, ] <- data$locations[which(B[t, ((1 + (t - 1) * n):(t * n))] == 1), ]
        loc_id[t] <- which(B[t, ((1 + (t - 1) * n):(t * n))] == 1)
      }
      locations_train_cv <- unique(locations_train_cv)
      locations_valid_cv[1, ] <- data$locations[which(B[T_train_cv + 1, ((1 + (T_train_cv) * n):((T_train_cv + 1) * n))] == 1), ]

      n_train <- length(unique(locations_train_cv[, 1]))
      B_train <- matrix(0, T_train_cv, n_train * T_train_cv)
      for (t in 1:T_train_cv) {
        j <- sample(c(1:n_train), siz = 1)
        B_train[t, j + (t - 1) * n_train] <- 1
      }
    }

    # Loop over each parameter setting
    num_cores <- detectCores() - 1
    results <- do.call(rbind, mclapply(1:nrow(param_grid), function(i) {
      res <- compute_posterior(
        y_train_cv, x_train_cv, x_valid_cv, B_train, T_train_cv, T_valid_cv, W, locations_train_cv, p, a_sigma, b_sigma,
        param_grid[i, ]$phi, param_grid[i, ]$nu, param_grid[i, ]$delta_beta, param_grid[i, ]$delta_z, alpha
      )
      log_prob <- dmvt(y_train_cv, as.vector(res$filt), (as(res$fil_cov, "matrix") + t(as(res$fil_cov, "matrix"))) / 2, df = res$df, log = TRUE) / length(y_train_cv)
      list(prediction = as.vector(res$filt)[(T_train_cv - 9):T_train_cv], log_prob = log_prob)
    }, mc.cores = num_cores))
    bayes_preds[, f, ] <- t(do.call(rbind, results[, 1]))
    log_probs[f, ] <- (do.call(rbind, results[, 2]))
  }

  Qmat <- matrix(0, nrow = nrow(param_grid), ncol = nrow(param_grid))
  for (k in 1:K) {
    Qmat <- Qmat + crossprod(bayes_preds[, k, ])
  }
  dvec <- numeric(nrow(param_grid))
  for (k in 1:K) {
    dvec <- dvec + 2 * t(bayes_preds[, k, ]) %*% y_valid[, k]
  }
  Amat <- diag(nrow(param_grid))
  Amat <- rbind(Amat, rep(1, nrow(param_grid)))
  bvec <- c(rep(0, nrow(param_grid)), 1)
  result_point <- solve.QP(2 * (Qmat + diag(1e-9, ncol(Qmat))), dvec, t(Amat), bvec, meq = 1)


  probs <- exp(log_probs)
  objective_function <- function(w, x) {
    ls <- sum(log((x %*% w)))
    return(-ls)
  }
  initial_weights <- rep(1 / ncol(probs), ncol(probs)) * 0.9
  ui <- rbind(rep(1, ncol(probs)), rep(-1, ncol(probs)), diag(ncol(probs)))
  ci <- c(0.98, -1, rep(0, ncol(probs)))
  result_dist <- constrOptim(
    theta = initial_weights,
    f = objective_function,
    grad = NULL,
    ui = ui,
    ci = ci,
    x = probs
  )

  # Unified function to calculate bp_y and bp_ydist
  calculate_values <- function(i) {
    params <- param_grid[i, ]
    res <- compute_posterior(y, x, x0, B_new, T, T0, W, locations, p, a_sigma, b_sigma, params$phi, params$nu, params$delta_beta, params$delta_z, alpha)
    list(
      mean = res$filt,
      z = res$fil_z,
      b = res$fil_b,
      ydist = dmvt(true, c(res$filt), (t(res$fil_cov) + res$fil_cov) / 2, df = res$df, log = TRUE) / length(true),
      sample = rmvn(1000, c(res$filt), res$b_sigma_star / res$a_sigma_star * (t(res$fil_cov) + res$fil_cov) / 2),
      sigma = rinvgamma(1000, shape = res$a_sigma_star, scale = res$b_sigma_star),
      sig = res$b_sigma_star / (res$a_sigma_star - 1)
    )
  }

  results <- mclapply(1:nrow(param_grid), calculate_values, mc.cores = num_cores)

  bp_y <- do.call(cbind, lapply(results, `[[`, "mean"))
  bp_z <- do.call(cbind, lapply(results, `[[`, "z"))
  bp_b <- do.call(cbind, lapply(results, `[[`, "b"))
  bp_sig <- do.call(cbind, lapply(results, `[[`, "sig"))
  bp_ydist <- do.call(cbind, lapply(results, `[[`, "ydist"))
  bp_sample <- array(do.call(c, lapply(results, `[[`, "sample")), dim = c(1000, T, nrow(param_grid)))
  bp_sigma <- array(do.call(c, lapply(results, `[[`, "sigma")), dim = c(1000, 1, nrow(param_grid)))
}


# Continuous time
# Spatiotemporal correlation kernel
K_phi <- function(t1, t2, gamma1, gamma2, phi1, phi2) {
  1 / (phi1 * abs(t1 - t2)^2 + 1) * exp(-phi2 * sqrt(sum(gamma1 - gamma2)^2) / sqrt(1 + phi1 * abs(t1 - t2)^2))
}
# Temporal correlation kernel for beta
C_xi <- function(t1, t2, xi) {
  exp(-xi^2 * abs(t1 - t2))
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
      KK[i, j] <- K_phi(times[i], times[j], locations[i, ], locations[j, ], phi1, phi2)
    }
  }


  CC <- kronecker(diag(p), C_beta)
  S <- bdiag(diag(n), delta_beta^2 * CC, delta_z^2 * KK)
  S <- S + diag(1e-10, nrow(S))
  S_inv <- solve(S)

  xx <- matrix(0, n, n * p)
  for (i in 1:p) {
    xx[, ((i - 1) * n + 1):(i * n)] <- diag(x[, i])
  }

  X <- rbind(
    cbind(xx, diag(n)),
    cbind(diag(n * p), matrix(0, p * n, n)),
    cbind(matrix(0, n, p * n), diag(n))
  )
  Y <- c(y, array(0, (p + 1) * n))

  SX <- solve(S, X)
  Sigma_inv <- t(X) %*% SX
  Sigma_inv <- (t(as.matrix(Sigma_inv)) + as.matrix(Sigma_inv)) / 2
  Sigma <- solve(Sigma_inv)
  m <- Sigma %*% t(SX) %*% Y

  a_sigma_star <- a_sigma + n / 2
  b_sigma_star <- b_sigma + t(Y - X %*% m) %*% S_inv %*% (Y - X %*% m) / 2
  sigma2_post <- rinvgamma(1, shape = a_sigma_star, scale = as.matrix(b_sigma_star))

  Sigma <- (t(as.matrix(Sigma)) + as.matrix(Sigma)) / 2
  theta_post <- rmvnorm(1, as.matrix(m), sigma2_post * Sigma)


  list(
    Sigma = Sigma, theta = c(as.matrix(m)),
    a_sigma_star = a_sigma_star, b_sigma_star = as.vector(b_sigma_star),
    y_filt = (X %*% m)[1:n, ], y_fil_cov = cbind(xx, diag(n)) %*% Sigma %*% t(cbind(xx, diag(n))),
    beta = c(as.matrix(m))[1:(p * n)], z = c(as.matrix(m))[(p * n + 1):(n * (p + 1))]
  )
}

results_con <- list()
for (current_T in T_values) {
  T <- current_T
  n <- T - 1

  # true parameters
  sigma <- 1
  phi <- 1 / 7
  nu <- 1
  delta_beta <- 1
  delta_z <- 1
  W <- diag(p)
  alpha <- 0.9

  data <- simulate_data(T, n, p, sigma, phi, nu, delta_beta, delta_z, W, alpha)

  test_id <- T
  train_id <- setdiff(c(1:T), test_id)
  y <- data$y[1:(test_id - 1)]
  y0 <- data$y[test_id]
  true <- data$true[1:(test_id - 1)]
  true0 <- data$true[test_id]
  xx <- data$x[1:(test_id - 1), (1:((test_id - 1) * 2))]
  xx0 <- matrix(data$x[test_id, (((test_id - 1) * 2) + 1):(test_id * 2)], nrow = 1)
  w <- data$w[1:((test_id - 1))]
  w0 <- data$w[test_id]
  z <- data$z[1:((test_id - 1) * n)]
  z0 <- data$z[((test_id - 1) * n + 1):(test_id * n)]
  beta <- data$beta[1:((test_id - 1) * p), ]
  beta0 <- data$beta[((test_id - 1) * p + 1):(test_id * p), ]
  T <- T - 1
  T0 <- 1
  B <- data$B[-test_id, ]
  B0 <- data$B[test_id, ]
  locations <- data$locations


  x <- matrix(0, ncol = p, nrow = T)
  for (i in 1:n) {
    x[i, ] <- xx[i, ((i * p - 1):(i * p))]
  }

  # Define k-fold cross-validation
  K <- 20
  folds <- lapply(1:K, function(k) c(1:(k * (T - 10) / K + 9)))

  # Parameter candidates
  phi1_values <- c(3, 1 / 10)
  phi2_values <- c(3, 1 / 10)
  xi_values <- c(3, 1 / 10)
  delta_beta_values <- c(3, 1 / 3)
  delta_z_values <- c(3, 1 / 3)

  # Define parameter grid
  param_grid <- expand.grid(
    phi1 = c(phi1_values),
    phi2 = c(phi2_values),
    xi = c(xi_values),
    delta_beta = c(delta_beta_values),
    delta_z = c(delta_z_values)
  )

  # Placeholder for Bayes predictors
  bayes_preds <- array(0, dim = c(10, K, nrow(param_grid)))
  log_probs <- array(0, dim = c(K, nrow(param_grid)))
  y_valid <- array(0, dim = c(10, K))

  # Loop over each fold
  for (f in 1:K) {
    # Split data
    {
      train_index <- folds[[f]]
      valid_index <- max(train_index) + 1

      T_train_cv <- length(train_index)
      T_valid_cv <- 1

      x_train_cv <- x[train_index, ]
      x_valid_cv <- x[valid_index, ]

      y_train_cv <- matrix(y[train_index], nrow = 1)
      y_valid_cv <- matrix(y[(T_train_cv - 9):T_train_cv])
      y_valid[, f] <- y_valid_cv

      true_train_cv <- matrix(true[train_index], nrow = 1)
      true_valid_cv <- matrix(true[valid_index])

      locations_train_cv <- locations[train_index, ]
    }

    num_cores <- detectCores() - 1
    results <- do.call(rbind, mclapply(1:nrow(param_grid), function(i) {
      res <- compute_posterior(
        y_train_cv, x_train_cv, 1:T_train_cv, locations_train_cv, T_train_cv, p, a_sigma, b_sigma,
        param_grid[i, ]$phi1, param_grid[i, ]$phi2, param_grid[i, ]$xi, param_grid[i, ]$delta_beta,
        param_grid[i, ]$delta_z
      )
      log_prob <- dmvt(c(y_valid_cv), as.vector(res$y_filt)[(T_train_cv - 9):T_train_cv], (res$y_fil_cov[(T_train_cv - 9):T_train_cv, (T_train_cv - 9):T_train_cv]), df = res$a_sigma_star, log = TRUE) / length(y_valid_cv)
      list(prediction = res$y_filt[(T_train_cv - 9):T_train_cv], log_prob = log_prob)
    }, mc.cores = num_cores))
    bayes_preds[, f, ] <- t(do.call(rbind, results[, 1]))
    log_probs[f, ] <- (do.call(rbind, results[, 2]))
  }


  Qmat <- matrix(0, nrow = nrow(param_grid), ncol = nrow(param_grid))
  for (k in 1:K) {
    Qmat <- Qmat + crossprod(bayes_preds[, k, ])
  }
  dvec <- numeric(nrow(param_grid))
  for (k in 1:K) {
    dvec <- dvec + 2 * t(bayes_preds[, k, ]) %*% y_valid[, k]
  }
  Amat <- diag(nrow(param_grid))
  Amat <- rbind(Amat, rep(1, nrow(param_grid)))
  bvec <- c(rep(0, nrow(param_grid)), 1)
  result_point <- solve.QP(2 * (Qmat + diag(1e-9, ncol(Qmat))), dvec, t(Amat), bvec, meq = 1)


  probs <- exp(log_probs)
  objective_function <- function(w, x) {
    ls <- sum(log((x %*% w)))
    return(-ls)
  }
  initial_weights <- rep(1 / ncol(probs), ncol(probs)) * 0.9
  ui <- rbind(rep(1, ncol(probs)), rep(-1, ncol(probs)), diag(ncol(probs)))
  ci <- c(0.98, -1, rep(0, ncol(probs)))
  result_dist <- constrOptim(
    theta = initial_weights,
    f = objective_function,
    grad = NULL,
    ui = ui,
    ci = ci,
    x = probs
  )

  # Unified function to calculate bp_y and bp_ydist
  calculate_values <- function(i) {
    res <- compute_posterior(
      y, x, 1:T, locations, T, p, a_sigma, b_sigma,
      param_grid[i, ]$phi1, param_grid[i, ]$phi2, param_grid[i, ]$xi, param_grid[i, ]$delta_beta,
      param_grid[i, ]$delta_z
    )
    log_prob <- dmvt(c(true), as.vector(res$y_filt), (res$y_fil_cov), df = res$a_sigma_star, log = TRUE) / length(true)

    list(
      mean = res$y_filt,
      z = res$z,
      beta = as.vector(t(matrix(res$beta, ncol = 2))),
      ydist = log_prob,
      sample = rmvn(1000, c(res$y_filt), res$b_sigma_star / res$a_sigma_star * (t(res$y_fil_cov) + (res$y_fil_cov)) / 2),
      sigma = rinvgamma(1000, res$a_sigma_star, res$b_sigma_star),
      sig = res$b_sigma_star / (res$a_sigma_star - 1)
    )
  }

  results <- mclapply(1:nrow(param_grid), calculate_values, mc.cores = num_cores)

  bp_y <- do.call(cbind, lapply(results, `[[`, "mean"))
  bp_z <- do.call(cbind, lapply(results, `[[`, "z"))
  bp_sig <- do.call(cbind, lapply(results, `[[`, "sig"))
  bp_beta <- do.call(cbind, lapply(results, `[[`, "beta"))
  bp_ydist <- do.call(cbind, lapply(results, `[[`, "ydist"))
  bp_sample <- array(do.call(cbind, lapply(results, `[[`, "sample")), dim = c(1000, T, nrow(param_grid)))
  bp_sigma <- array(do.call(cbind, lapply(results, `[[`, "sigma")), dim = c(1000, 1, nrow(param_grid)))
}


### dlm
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
model_file <- "~/dlm.stan"
stan_model <- stan_model(model_file)

results_linear <- list()
for (current_T in T_values) {
  T <- current_T
  n <- T - 1

  # true parameters
  sigma <- 1
  phi <- 1 / 7
  nu <- 1
  delta_beta <- 1
  delta_z <- 1
  W <- diag(p)
  alpha <- 0.9

  data <- simulate_data(T, n, p, sigma, phi, nu, delta_beta, delta_z, W, alpha)

  test_id <- T
  train_id <- setdiff(c(1:T), test_id)
  y <- data$y[1:(test_id - 1)]
  y0 <- data$y[test_id]
  true <- data$true[1:(test_id - 1)]
  true0 <- data$true[test_id]
  xx <- data$x[1:(test_id - 1), (1:((test_id - 1) * 2))]
  xx0 <- matrix(data$x[test_id, (((test_id - 1) * 2) + 1):(test_id * 2)], nrow = 1)
  w <- data$w[1:((test_id - 1))]
  w0 <- data$w[test_id]
  z <- data$z[1:((test_id - 1) * n)]
  z0 <- data$z[((test_id - 1) * n + 1):(test_id * n)]
  beta <- data$beta[1:((test_id - 1) * p), ]
  beta0 <- data$beta[((test_id - 1) * p + 1):(test_id * p), ]
  T <- T - 1
  T0 <- 1
  n <- T
  B <- data$B[-test_id, ]
  B0 <- data$B[test_id, ]

  x <- matrix(0, ncol = p, nrow = T)
  for (i in 1:n) {
    x[i, ] <- xx[i, ((i * p - 1):(i * p))]
  }


  # Prepare the data list for Stan
  stan_data <- list(T = T, K = p, y = y, x = x)

  fit <- sampling(stan_model, data = stan_data, iter = 3000, chains = 4)
  beta_est <- apply(rstan::extract(fit)$state, c(2, 3), mean)
  y_est <- rep(NA, T)
  for (t in 1:T) {
    y_est[t] <- beta_est[t, ] %*% x[t, ]
  }

  stack_sample <- rmvn(1000, y_est, diag(mean(rstan::extract(fit)$sigma[2001:3000])^2, T))
  sigma_sample <- rstan::extract(fit)$sigma[2001:3000]
  stack_mean <- colMeans(stack_sample)
  sigma_mean <- mean(sigma_sample)
  D_bar <- 0
  for (i in 1:length(sigma_sample)) {
    D_bar <- D_bar - 2 * dmvn(y, stack_sample[i, ], sigma_sample[i] * diag(T), log = TRUE) / length(sigma_sample)
  }
  D_theta_bar <- -2 * dmvn(y, stack_mean, sigma_mean * diag(T), log = TRUE)
  p_D <- D_bar - D_theta_bar

  x <- matrix(0, nrow = 1, ncol = length(sigma_sample))
  for (i in 1:length(sigma_sample)) {
    x[, i] <- dmvn(y, stack_sample[i, ], sigma_sample[i] * diag(T), log = TRUE)
  }
  lppd <- sum(log(rowMeans(exp(x))))
  pWAIC1 <- 2 * sum(log(rowMeans(exp(x))) - rowMeans(x))
  pWAIC2 <- sum(.rowVars(x))
}
