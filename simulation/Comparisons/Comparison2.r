# import functions
source("functions/continuous_time.r")


# setting
T_values <- c(50, 70) + 1
p <- 2
W <- diag(p)
a_sigma <- 1 / 2
b_sigma <- 1 / 2


# continuous time
results_con <- list()
for (current_T in T_values) {


  T <- current_T
  n <- T - 1

  # true parameters
  sigma <- 1
  phi1 <- 1 / 2
  phi2 <- 1 / 2
  xi <- 1 / 2
  delta_beta <- 1
  delta_z <- 1

  data <- simulate_data(n, p, sigma, phi1, phi2, xi, delta_beta, delta_z)


  test_id <- T
  train_id <- setdiff(c(1:T), test_id)
  y <- data$y[-test_id]
  true <- data$true[-test_id]
  y0 <- data$y[test_id]
  x <- data$x[-test_id, ]
  z <- data$z[-test_id]
  beta <- data$beta[, -test_id]
  T <- data$times[-test_id]
  T0 <- data$times[test_id]
  G <- data$locations[-test_id, ]
  locations <- data$locations[-test_id, ]
  n <- length(T)
  n0 <- length(T0)

  # Define k-fold cross-validation
  K <- 20
  folds <- lapply(1:K, function(k) c(1:(k * (n - 10) / K + 9)))

  # Parameter candidates
  phi1_values <- c(1, 1 / 4)
  phi2_values <- c(1, 1 / 4)
  xi_values <- c(1, 1 / 4)
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

      locations_train_cv <- locations[train_index, ]
    }


    # Loop over each parameter setting
    num_cores <- detectCores() - 1
    results <- do.call(rbind, mclapply(1:nrow(param_grid), function(i) {
      res <- compute_posterior(
        y_train_cv, x_train_cv, 1:T_train_cv, locations_train_cv, T_train_cv, p, a_sigma, b_sigma,
        param_grid[i, ]$phi1, param_grid[i, ]$phi2, param_grid[i, ]$xi, param_grid[i, ]$delta_beta,
        param_grid[i, ]$delta_z
      )
      log_prob <- dmvt(c(y_valid_cv), as.vector(res$y_filt)[(T_train_cv - 9):T_train_cv],
        (res$y_fil_cov[(T_train_cv - 9):T_train_cv, (T_train_cv - 9):T_train_cv]),
        df = res$a_sigma_star, log = TRUE
      ) / length(y_valid_cv)
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
  result_point <- solve.QP(2 * Qmat + diag(2e-9, ncol(Qmat)), dvec, t(Amat), bvec, meq = 1)


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
      y, x, T, locations, length(T), p, a_sigma, b_sigma,
      param_grid[i, ]$phi1, param_grid[i, ]$phi2, param_grid[i, ]$xi, param_grid[i, ]$delta_beta,
      param_grid[i, ]$delta_z
    )
    log_prob <- dmvt(c(true), as.vector(res$y_filt), (res$y_fil_cov), df = res$a_sigma_star, log = FALSE) / length(true)

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
  bp_sample <- array(do.call(c, lapply(results, `[[`, "sample")), dim = c(1000, n, nrow(param_grid)))
  bp_sigma <- array(do.call(c, lapply(results, `[[`, "sigma")), dim = c(1000, 1, nrow(param_grid)))

}


# discretre time
compute_posterior <- function(y, x, x0, B, T, T0, W, locations, p, a_sigma, b_sigma, phi, nu, delta_beta, delta_z, alpha) {
  # locations_ <- unique(locations)
  n <- nrow(locations)

  # Generate covariance matrix S for z using K_phi
  KK <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      KK[i, j] <- matern_kernel(locations[i, ], locations[j, ], phi, nu)
    }
  }

  zero_vector <- matrix(0, nrow = T - 1, ncol = 1)
  zero_vector_transpose <- t(zero_vector)
  identity_matrix <- diag(T - 1)
  single_zero <- matrix(0, nrow = 1, ncol = 1)
  A <- rbind(
    cbind(zero_vector_transpose, single_zero),
    cbind(alpha * identity_matrix, zero_vector)
  )
  SnA <- solve(diag(T) - A)
  SIP <- Matrix(kronecker(SnA, diag(p)))
  SIN <- Matrix(kronecker(SnA, diag(n)))

  IP <- Matrix(kronecker(diag(T) - A, diag(p)))
  IN <- Matrix(kronecker(diag(T) - A, diag(n)))

  KTW <- Matrix(kronecker(diag(T), solve(W)))
  KTK <- Matrix(kronecker(diag(T), solve(KK)))
  S_inv <- Matrix(bdiag(diag(T), delta_beta^
    {
      -2
    } * t(IP) %*% KTW %*% (IP), delta_z^
    {
      -2
    } * t(IN) %*% KTK %*% (IN)))

  X <- Matrix(rbind(
    cbind(x, B),
    cbind(diag(T * p), matrix(0, T * p, ncol(B))),
    cbind(matrix(0, n * T, p * T), diag(n * T))
  ))
  Y <- c(y, array(0, (p + n) * T))

  SX <- S_inv %*% X
  Sigma_inv <- t(X) %*% SX
  Sigma <- as.matrix(solve(Sigma_inv))
  m <- Sigma %*% t(SX) %*% Y
  a_sigma_star <- a_sigma + T / 2
  b_sigma_star <- b_sigma + t(Y - X %*% m) %*% S_inv %*% (Y - X %*% m) / 2

  Sigma <- as.matrix(Sigma)

  list(
    Sigma = Sigma, theta = c(as.matrix(m)), filt = c((X %*% m)[1:T, ]),
    fil_cov = cbind(x, B) %*% Sigma %*% t(cbind(x, B)), fil_z = c(B %*% m[(p * T + 1):(n * T + p * T)]), fil_b = c(m[1:(p * T)]),
    a_sigma_star = a_sigma_star, b_sigma_star = as.vector(b_sigma_star), df = T / 2
  )
}

results_dis <- list()
for (current_T in T_values) {


  T <- current_T
  n <- T - 1

  # true parameters
  sigma <- 1
  phi1 <- 1 / 2
  phi2 <- 1 / 2
  xi <- 1 / 2
  delta_beta <- 1
  delta_z <- 1

  data <- simulate_data(T, p, sigma, phi1, phi2, xi, delta_beta, delta_z)

  test_id <- T
  train_id <- setdiff(c(1:T), test_id)
  y <- data$y[1:(test_id - 1)]
  true <- data$true[1:(test_id - 1)]
  w <- data$w[1:((test_id - 1))]
  z <- data$z[1:((test_id - 1))]
  beta <- data$beta[1:((test_id - 1) * p)]
  T <- T - 1

  B <- matrix(0, n, n * n)
  for (t in 1:n) {
    B[t, n * (t - 1) + t] <- 1
  }

  locations <- data$locations[1:T, ]

  xx <- data$x[1:(test_id - 1), ]
  x <- matrix(0, n, n * p)
  for (t in 1:n) {
    x[t, (((t - 1) * p + 1):(t * p))] <- xx[t, ]
  }


  # Define k-fold cross-validation
  K <- 20
  folds <- lapply(1:K, function(k) c(1:(k * (T - 10) / K + 9)))

  # Parameter candidates
  phi_values <- c(2, 1 / 2)
  nu_values <- c(2, 1 / 2)
  delta_beta_values <- c(1 / 2, 1 / 10)
  delta_z_values <- c(1 / 2, 1 / 10)

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

      locations_train_cv <- locations[train_index, ]

      n_train <- length(unique(locations_train_cv[, 1]))
      B_train <- matrix(0, T_train_cv, n_train * T_train_cv)
      for (t in 1:T_train_cv) {
        B_train[t, t + (t - 1) * n_train] <- 1
      }
    }

    # Loop over each parameter setting
    num_cores <- detectCores() - 1
    results <- do.call(rbind, mclapply(1:nrow(param_grid), function(i) {
      res <- compute_posterior(
        y_train_cv, x_train_cv, x_valid_cv, B_train, T_train_cv, T_valid_cv, W, locations_train_cv, p, a_sigma, b_sigma,
        param_grid[i, ]$phi, param_grid[i, ]$nu, param_grid[i, ]$delta_beta, param_grid[i, ]$delta_z, 1
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
  result_point <- solve.QP(2 * Qmat + diag(1e-7, ncol(Qmat)), dvec, t(Amat), bvec, meq = 1)

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
    res <- compute_posterior(y, x, x0, B, T, T0, W, locations, p, a_sigma, b_sigma, params$phi, params$nu, params$delta_beta, params$delta_z, 1)
    list(
      mean = res$filt,
      z = res$fil_z,
      b = res$fil_b,
      ydist = dmvt(true, c(res$filt), (t(res$fil_cov) + res$fil_cov) / 2, df = res$df, log = TRUE) / length(true),
      sample = rmvn(1000, c(res$filt), res$b_sigma_star / res$a_sigma_star * (t(res$fil_cov) + res$fil_cov) / 2),
      sigma = rinvgamma(1000, res$a_sigma_star, res$b_sigma_star),
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
  phi1 <- 1 / 2
  phi2 <- 1 / 2
  xi <- 1 / 2
  delta_beta <- 1
  delta_z <- 1

  data <- simulate_data(n, p, sigma, phi1, phi2, xi, delta_beta, delta_z)

  test_id <- T
  train_id <- setdiff(c(1:T), test_id)
  y <- data$y[-test_id]
  true <- data$true[-test_id]
  y0 <- data$y[test_id]
  x <- data$x[-test_id, ]
  z <- data$z[-test_id]
  beta <- data$beta[, -test_id]
  T <- data$times[-test_id]
  T0 <- data$times[test_id]
  G <- data$locations[-test_id, ]
  locations <- data$locations[-test_id, ]
  n <- length(T)
  n0 <- length(T0)


  # Prepare the data list for Stan
  stan_data <- list(T = n, K = p, y = y, x = x)

  fit <- sampling(stan_model, data = stan_data, iter = 3000, chains = 3)
  beta_est <- apply(rstan::extract(fit)$state, c(2, 3), mean)
  y_est <- rep(NA, n)
  for (t in 1:n) {
    y_est[t] <- beta_est[t, ] %*% x[t, ]
  }

  stack_sample <- rmvn(1000, y_est, diag(mean(rstan::extract(fit)$sigma[2001:3000])^2, n))
  sigma_sample <- rstan::extract(fit)$sigma[2001:3000]
  stack_mean <- colMeans(stack_sample)
  sigma_mean <- mean(sigma_sample)
  D_bar <- 0
  for (i in 1:length(sigma_sample)) {
    D_bar <- D_bar - 2 * dmvn(y, stack_sample[i, ], sigma_sample[i] * diag(n), log = TRUE) / length(sigma_sample)
  }
  D_theta_bar <- -2 * dmvn(y, stack_mean, sigma_mean * diag(n), log = TRUE)
  p_D <- D_bar - D_theta_bar

  x <- matrix(0, nrow = 1, ncol = length(sigma_sample))
  for (i in 1:length(sigma_sample)) {
    x[, i] <- dmvn(y, stack_sample[i, ], sigma_sample[i] * diag(n), log = TRUE)
  }
  lppd <- sum(log(rowMeans(exp(x))))
  pWAIC1 <- 2 * sum(log(rowMeans(exp(x))) - rowMeans(x))
  pWAIC2 <- sum(.rowVars(x))

}
