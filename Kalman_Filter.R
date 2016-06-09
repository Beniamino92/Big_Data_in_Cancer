####### Kalman Filter & Smoother #########

require(Matrix)

##### Function: get time variable matrix H ####
get.H <- function(omega.j, T) {
  H <- list()
  for(t in 1:T) {
    H[[t]] <- matrix(c(rbind(cos(2*pi*omega.j * t), sin(2*pi*omega.j * t))), nrow = 1)
  }
  return(H)
}

get.H.2 <- function(omega.j, time) {
  H <- list()
  for(t in 1:length(time)) {
    H[[t]] <- matrix(c(rbind(cos(omega.j * time[t]), 
                             sin(omega.j * time[t]))), nrow = 1)
  }
  return(H)
}


#### Function: Kalman Filter #####
Kalman.Filter <- function(y, F, H, G, Q, sigma2) {
  
  # Length of the time series
  T <- length(y)
  
  # Number of states  [(Trend + N.Harmonics) * 2]
  n <- dim(F)[1]
  
  
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


#### Function: Kalman - Smoother ####
Kalman.Smoother <- function(filtered.process, plot = FALSE) {
  
  if(class(filtered.process)!="Kalman.Filter") {
    stop("object \'filtered.process\' is not of class \'Kalman.Filter\'")
  }
  
  ### Getting objects 
  states <- filtered.process$states
  covariance.states <- filtered.process$covariance.states
  prediction.P <- filtered.process$prediction.P
  y <- filtered.process$y 
  F <- filtered.process$F
  H <- filtered.process$H
  G <- filtered.process$G
  sigma2 <- filtered.process$sigma2
  Q <- filtered.process$Q; Q.r <- Q/sigma2
  N <- length(y)
  R <- length(H[[1]])
  
  ### Initializations
  I <- diag(R)
  # L
  L <- list()
  L[[N]] <- rep(0, R)
  # x.smoothed
  smoothed.states <- list()
  smoothed.states[[N]] <- states[[N]]
  # P.smoothed
  smoothed.covariance <- list()
  smoothed.covariance[[N]] <- covariance.states[[N]]
  
  #### Backward Pass Smoothing Equations ####
  for(t in (N - 1):1) {
    
    # Lagrange multiplier
    A <- I - covariance.states[[t + 1]] %*% 
      t(H[[t + 1]]) %*% H[[t + 1]]
    B <- t(F) %*% L[[t + 1]] - (t(H[[t + 1]]) %*% 
                                  (y[t + 1] - H[[t + 1]] %*% states[[t + 1]]))
    L[[t]] <- t(A) %*% B
    
    # Smoothed states
    smoothed.states[[t]] <- solve(F) %*% 
      (smoothed.states[[t + 1]] + G %*% Q.r %*% t(G) %*% L[[t]])
    
    # Smoothed covariance states
    C <- solve(smoothed.covariance[[t + 1]] - prediction.P[[t + 1]],
               prediction.P[[t + 1]])
    smoothed.covariance[[t]] <- covariance.states[[t]] +
      solve(covariance.states[[t]] %*% t(F), prediction.P[[t + 1]]) %*%
      C %*% F %*% covariance.states[[t]]
  }
  
  #### Plotting ####
  if(plot == TRUE) {
    fitted <- fitted.DHR(H, states)
    smooth.fitted <- fitted.DHR(H, smoothed.states) 
    
    par(mfrow = c(1, 1))
    plot(1:N, y, col = "grey", pch = 20)
    points(1:N, y, type = "o", col = "grey", pch = 20)
    lines(1:N, fitted, type = "o", col = "green")
    lines(1:N, smooth.fitted, type = "o", col = "purple", pch = 20)
  }
  
  return(list(smoothed.states = smoothed.states, 
              smoothed.covariance = smoothed.covariance,
              F = F, 
              H = H))
}


# Function: k-step ahead forecasting on filtered.process
Kalman.Forecasting <- function(filtered.process, t, K) {
  
  # Getting objects
  F <- filtered.process$F
  H <- filtered.process$H
  Q <- filtered.process$Q
  G <- filtered.process$G
  sigma2 <- filtered.process$sigma2
  Q.r <- Q/sigma2
  
  a <- list()
  R <- list()
  forecast <- c()
  forecast.variance <- c()
  a[[1]] <- filtered.process$states[[t]]
  R[[1]] <- filtered.process$covariance.states[[t]]
  
  # k-step ahead prediction states
  for(k in 2:(K+1)) {
    a[[k]] <- F %*% a[[k - 1]]
    R[[k]] <- F %*% R[[k - 1]] %*% t(F) + G %*% Q.r %*% t(G)
  }
  a <- a[-1]
  # k-step ahead prediction observation
  for(k in 1:K) {
    forecast[k] <- as.vector(H[[t + k]] %*% a[[k]])
    forecast.variance[k] <- as.vector(H[[t + k]] %*% R[[k]] %*% t(H[[t + k]]) +
                                        sigma2)
  }
  return(list(forecast = forecast, forecast.variance = forecast.variance))
}

# Function: k-step ahead forecasting on smoothed.process
Kalman.Forecasting.Smooth <- function(smoothed.process, t, K) {
  
  # Getting objects
  F <- smoothed.process$F
  H <- smoothed.process$H
  a <- list()
  forecast <- c()
  a[[1]] <- smoothed.process$smoothed.states[[t]]
  
  # k-step ahead prediction states
  for(k in 2:(K + 1)) {
    a[[k]] <- F %*% a[[k - 1]]
  }
  a <- a[-1]
  # k-step ahead prediction observations
  for(k in 1:K) {
    forecast[k] <- as.vector(H[[t + k]] %*% a[[k]])
  }
  
  return(forecast)
}


#### Function: computing the fitted DHR ####
fitted.DHR <- function(H, states) {
  
  T <- length(H)
  fitted <- c()
  
  for(t in 1:T) {
    fitted[t] <- as.numeric(H[[t]] %*% states[[t]])
  }
  return(fitted)
}

