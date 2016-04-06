#### THE COMPLETE DHR ESTIMATION ALGORITHM 

require(astsa)

##################### STEP 1 ##################

#### Estimate AR(p) spectrum of the observation process
estimate.AR <- function(data, p, n.freq) {
  
  # Fitting AR(p) to data, getting the coefficients and the residual variance
  fit <- arima(data, order = c(p, 0, 0), method = "ML")
  coef <- fit$coef[1:p]
  sigma2 <- fit$sigma2
  
  # Getting spectrum AR(p)
  ar.fit <- arma.spec(ar = coef, n.freq = n.freq, var.noise = fit$sigma2)
  
  # Not considering frequency 0
  spectrum <- ar.fit$spec[-1]
  freq <- ar.fit$freq[-1]
  
  return(list(spectrum = spectrum, frequency = freq , sigma2 = sigma2))
}



######################### STEP 2 #####################


#### Find Linear Least Square estimate of NVR parameter vector

NVR.LinearLS <- function(empirical.spectrum, sigma2, 
                         frequency, omega.j, structures) {
  
  # Number of frequencies
  T <- length(frequency)
  # Number of harmonic components
  R <- length(omega.j) - 1
  
  # Preparing data for evaluating LLS NVR estimate
  Y <- empirical.spectrum - (sigma2/(2*pi))
  X <- matrix(NA, nrow = T, ncol = (R + 1))
  
  for(r in 1:(R + 1)) {
    X[, r] <- S(frequency, omega.j[r], type = structures[r])
  }
  
  X <- sigma2 * X
  
  # Least square estimate
  NVR.estimate <- solve(t(X) %*% X, t(X) %*% Y)
  
  return(NVR.estimate)
}



######################### STEP 3 #####################


#### Function to evaluate the logarithmic error.

E <- function(NVR, freq, empirical.spectrum, omega.j, structures, sigma2) {
  
  
  # Number of frequencies
  T <- length(freq)
  # Number of harmonic components
  R <- length(omega.j) - 1
  
  # Matrix T * (R + 1) with components X[i, j] <- S(omega_i, omega_j)
  X <- matrix(NA, nrow = T, ncol = (R + 1))
  for(r in 1:(R + 1)) {
    X[, r] <- S(freq, omega.j[r], type = structures[r])
  }
  
  fitted <- sigma2 * ((X %*% NVR) + (1/(2*pi)))
  
  # Evaluating E(sigma)_L; 
  error <- sum((log(empirical.spectrum) - 
                  log(fitted))^2)
  
  return(error)
}

