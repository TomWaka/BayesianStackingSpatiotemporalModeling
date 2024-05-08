# Load the necessary script from the functions directory
source("functions/continuous_time.r")
library(dplyr)


# Function to preprocess data: loads data, filters, and scales time
preprocess_data <- function(data_path, pid, start_time_offset, end_time_offset) {
  data <- read.csv(data_path)
  filtered_data <- dplyr::filter(data, PID == pid) %>%
    dplyr::select(Latitude, Longitude, time, MAG, DistFromHome, Slope, DTP, NDVI)
  
  filtered_data$time <- as.POSIXct(filtered_data$time, format="%Y-%m-%d %H:%M:%S")
  
  start_time <- min(filtered_data$time) + as.difftime(start_time_offset, units="days")
  end_time <- start_time + as.difftime(end_time_offset, units="hours")
  
  dataset <- subset(filtered_data, time >= start_time & time <= end_time)
  dataset <- dataset[1:650,]
  
  dataset$time <- as.POSIXct(dataset$time, format="%Y-%m-%d %H:%M:%S PDT", tz="America/Los_Angeles")
  time_seconds <- as.numeric(difftime(dataset$time, min(dataset$time), units="secs"))
  time_scaled <- (time_seconds - min(time_seconds)) / (max(time_seconds) - min(time_seconds))
  
  return(list(dataset = dataset, time_scaled = time_scaled))
}

# Function to perform grid search over a set of parameters
run_grid_search <- function(dataset, time_scaled, K, param_grid) {
  test_id <- c(seq(2, nrow(dataset), by=5))
  train_id <- setdiff(c(1:nrow(dataset)), test_id)
  
  y  <- log(dataset$MAG[-test_id])
  y0 <- log(dataset$MAG[test_id])
  x  <- cbind(dataset$Slope[-test_id], dataset$NDVI[-test_id])
  x0 <- cbind(dataset$Slope[test_id], dataset$NDVI[test_id])
  T  <- time_scaled[-test_id]
  T0 <- time_scaled[test_id]
  G  <- cbind(dataset$Latitude[-test_id], dataset$Longitude[-test_id])
  G0 <- cbind(dataset$Latitude[test_id], dataset$Longitude[test_id])
  n  <- length(T)
  n0 <- length(T0)
  p  <- ncol(x)
  
  a_sigma <- 1/2
  b_sigma <- 1/2
  
  folds <- lapply(1:K, function(start) c(1:n)[seq(start, length(train_id), by= K)])
  
  bayes_preds <- array(0, dim=c(n/K, K, nrow(param_grid)))
  log_probs <- array(0, dim=c(K, nrow(param_grid)))
  y_valid   <- array(0, dim=c(n/K, K))
  
  num_cores <- detectCores() - 1
  
  for(f in 1:K) {
    print(f)
    
    valid_index <- folds[[f]]
    train_index <- setdiff(1:length(train_id), valid_index)
    
    G_train_cv <- G[train_index,]
    G_valid_cv <- G[valid_index,]
    
    T_train_cv <- T[train_index]
    T_valid_cv <- T[valid_index]
    
    x_train_cv <- x[train_index,]
    x_valid_cv <- x[valid_index,]
    
    y_train_cv  <- y[train_index]
    y_valid_cv  <- y[valid_index]
    y_valid[,f] <- y_valid_cv
    
    results <- do.call(rbind, mclapply(1:nrow(param_grid), function(i) {
      res <- pred_dist_continuous(y_train_cv, x_train_cv, x_valid_cv, T_train_cv, T_valid_cv, G_train_cv, G_valid_cv, length(y_train_cv), p, a_sigma, b_sigma, 
                                  param_grid[i,]$phi1, param_grid[i,]$phi2, param_grid[i,]$xi, param_grid[i,]$delta_beta, 
                                  param_grid[i,]$delta_z)
      log_prob <- dmvt(y_valid_cv, as.vector(res$mean), (as(res$variance,"matrix")+t(as(res$variance,"matrix")))/2 , df=res$df , log=TRUE) / length(y_valid_cv)
      list(prediction = res$mean, log_prob = log_prob)
    }, mc.cores=num_cores))
    
    bayes_preds[, f, ] <- t(do.call(rbind, results[,1]))
    log_probs[f, ] <- (do.call(rbind, results[,2]))
  }
  
  return(list(bayes_preds = bayes_preds, log_probs = log_probs, y_valid = y_valid, y = y, y0 = y0, x = x, x0 = x0, T = T, T0 = T0, G = G, G0 = G0, n = n, n0 = n0, p = p))
}

# Function to determine the point weights
determine_point_weights <- function(bayes_preds, param_grid) {
  Qmat <- matrix(0, nrow = nrow(param_grid), ncol = nrow(param_grid))
  for (k in 1:K) {
    Qmat <- Qmat + crossprod(bayes_preds[, k, ])
  }
  
  dvec <- numeric(nrow(param_grid))
  for (k in 1:K) {
    dvec <- dvec + 2 * t(bayes_preds[, k, ]) %*% y_valid[,k]
  }
  
  Amat <- diag(nrow(param_grid))
  Amat <- rbind(Amat, rep(1, nrow(param_grid)))
  bvec <- c(rep(0, nrow(param_grid)), 1)
  
  result_point <- solve.QP(2*Qmat+diag(1e-10, ncol(Qmat)), dvec, t(Amat), bvec, meq = 1)
  return(result_point$solution)
}



# Main processing
pid <- 285
start_time_offset <- 1.03
end_time_offset <- 17

# Preprocess data and set up the parameter grid
processed_data <- preprocess_data(data_path, pid, start_time_offset, end_time_offset)
dataset <- processed_data$dataset
time_scaled <- processed_data$time_scaled

# Define parameter values and run grid search
K <- 10
phi1_values <- c(100, 1)
phi2_values <- c(100, 1)
xi_values  <- c(100, 1)
delta_beta_values <- c(20, 1/5)
delta_z_values    <- c(20, 1/5)

param_grid <- expand.grid(phi1 = c(phi1_values),
                          phi2 = c(phi2_values),
                          xi = c(xi_values),
                          delta_beta = c(delta_beta_values),
                          delta_z = c(delta_z_values))

grid_search_results <- run_grid_search(dataset, time_scaled, K, param_grid)
bayes_preds <- grid_search_results$bayes_preds
log_probs <- grid_search_results$log_probs
y_valid <- grid_search_results$y_valid
y <- grid_search_results$y
y0 <- grid_search_results$y0
x <- grid_search_results$x
x0 <- grid_search_results$x0
T <- grid_search_results$T
T0 <- grid_search_results$T0
G <- grid_search_results$G 
G0 <- grid_search_results$G0
n <- grid_search_results$n
n0 <- grid_search_results$n0
p <- grid_search_results$p

point_weights <- determine_point_weights(bayes_preds, param_grid)

# computation and output
calculate_values <- function(i) {
  params <- param_grid[i,]
  res <- pred_dist_continuous(y, x, x0, T, T0, G, G0, n, p, a_sigma, b_sigma, params$phi1, params$phi2, params$xi, params$delta_beta, params$delta_z)
  list(mean = res$mean, 
       z = res$z,
       z_sampled = res$z_sampled,
       yvar = (t(res$variance) + res$variance) / 2,
       ydist = dmvt(y0, c(res$mean), (t(res$variance) + res$variance) / 2, df = n / 2, log = TRUE) / n0)
}

num_cores <- detectCores() - 1
results <- mclapply(1:nrow(param_grid), calculate_values, mc.cores = num_cores)

bp_y <- do.call(cbind, lapply(results, `[[`, "mean"))
bp_z <- do.call(cbind, lapply(results, `[[`, "z"))
bp_zsampled <- do.call(cbind, lapply(results, `[[`, "z_sampled"))
bp_ydist <- do.call(cbind, lapply(results, `[[`, "ydist"))
bp_yvar <- do.call(cbind, lapply(results, `[[`, "yvar"))

print(mean((point_weights %*% t(bp_y) - y0)^2))
print(point_weights %*% t(bp_ydist))

