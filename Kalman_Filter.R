####### Kalman Filter:  Forward Pass Filtering #########

require(Matrix)

# Function to create the time variable vector H
get.H <- function(omega.j, T) {
  H <- list()
  for(t in 1:T) {
    H[[t]] <- matrix(c(rbind(cos(omega.j * t), sin(omega.j * t))), nrow = 1)
  }
  return(H)
}


# Function to implement Forward Pass Filtering Equations (Kalman)

Kalman.Filter <- function(y, F, H, G, Q, sigma2) {
  
  # Length of the time series
  T <- length(y)
  
  # Number of states  [(Trend + N.Harmonics) * 2]
  n <- length(omega.j) * 2
  
  
  # Checking dimensionality
  if(!( all(dim(F) == c(n, n)) & all(dim(G) == c(n, n)) & 
        all(dim(Q) == c(n, n)) & all(dim(H[[T]]) == c(1, n)))) {
    stop("Matrix dimensions are incorrect")
  }
  
  # States
  states <- list()
  states[[1]] <- matrix(rep(0, n), nrow = n, ncol = 1)
  
  # Covariances matrix
  covariance.states <- list()
  covariance.states[[1]] <- diag(n)
  
  # Q.r: matrix containing (scaled) system noise 
  Q.r <- Q / sigma2
  
  # Intermediate calculation steps
  prediction.x <- list()
  prediction.P <- list()
  
  
  ### Forward Pass Filtering Equations ###
  
  for(t in 2:T) {
    # Prediction
    prediction.x[[t]] <- F %*% states[[t - 1]]
    
    prediction.P[[t]] <- F %*% covariance.states[[t - 1]] %*% t(F) +
      G %*% Q.r %*% t(G)
    
    # Correction
    states[[t]] <- prediction.x[[t]] + prediction.P[[t]] %*% t(H[[t]]) %*% 
      solve(1 + H[[t]] %*% prediction.P[[t]] %*% t(H[[t]]),
            y[t] - H[[t]] %*% prediction.x[[t]])
    
    covariance.states[[t]] <- prediction.P[[t]] - prediction.P[[t]] %*% t(H[[t]]) %*% 
      solve(1 + H[[t]] %*% prediction.P[[t]] %*% t(H[[t]]),
            H[[t]] %*% prediction.P[[t]])
  }
  
  filtered.process <- list(states = states, 
                           covariance.states = covariance.states,
                           prediction.x = prediction.x, 
                           prediction.P = prediction.P,
                           y = y, F = F, H = H, 
                           G = G, Q = Q, sigma2 = sigma2)
  
  class(filtered.process) <- "Kalman.Filter"
  
  return(filtered.process)
}


# Function to compute the fitted value of y, given the states,
# and the model. 

fitted.DHR <- function(H, states) {
  
  T <- length(H)
  fitted <- c()
  
  for(t in 1:T) {
    fitted[t] <- as.numeric(H[[t]] %*% states[[t]])
  }
  return(fitted)
}

