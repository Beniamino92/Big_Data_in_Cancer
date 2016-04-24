# Dynamic Harmonic Regression: Kalman Filter


# We can put H as a list; each component represent the matrix
# at time t. 

# Forward Pass Filtering Equations:

DHR.Filtering <- function(y, H, F, G, Q, sigma, T) {
  
  for(t in 1:T) {
    
    H <- H[[t]]
    y <- y[t]
    filtered.states <- list()
    filtered.variances <- list()
    
    # Prediction:
    x_prediction <- F %*% x_previous_state
    P_prediction <- F %*% P_previous_state %*% t(F) + G %*% Q %*% t(G)
    
    # Correction:
    # Remember to change the computation of the inverse.
    x_correction <- x_prediction + P_prediction %*% t(H) %*% 
      solve(1 + H %*% P_prediction %*% t(H)) %*% (y - H %*% x_prediction)
    
    P_correction <- P_prediction - P_prediction %*% t(H) %*% 
      solve(1 + H %*% P_prediction %*% t(H)) %*% (H %*% P_previous_state)
    
    
    filtered.states[[t]] <- x_correction
    filtered.variances[[t]] <- P_correction
    
    X_previous_state <- x_correction
    P_previous_state <- P_correction
  }
}

H <- H[[t]]
y <- y[t]

filtered.states <- list()

# Prediction:
x_prediction <- F %*% x_previous_state
P_prediction <- F %*% P_previous_state %*% t(F) + G %*% Q %*% t(G)


# Correction:
# Remember to change the computation of the inverse.
x_correction <- x_prediction + P_prediction %*% t(H) %*% 
  solve(1 + H %*% P_prediction %*% t(H)) %*% (y - H %*% x_prediction)

P_correction <- P_prediction - P_prediction %*% t(H) %*% 
  solve(1 + H %*% P_prediction %*% t(H)) %*% (H %*% P_previous_state)


filtered.states[[t]] <- x_correction
filtered.variances[[t]] <- P_correction

X_previous_state <- x_correction
P_previous_state <- P_correction


