#### For now we assume that n is even.

# Function: Periodogram
periodogram <- function(data) {
  n <- length(data)
  I <- (abs(fft(data))^2)/(n*2*pi)
  return(I)
}

# Function: Gaussian Kernel
gaussian.kernel <- function(x, b = 1, scaled = FALSE) {
  
  if(scaled == TRUE) {
    out <- (1/b) * (1/sqrt(2*pi)) * exp(-((x/b)^2)/2)
    return(out)
  }
  
  out <- (1/sqrt(2*pi)) * exp(-(x^2)/2)
  return(out)
}

# Function: Smoothed periodogram for single frequency
single.smoothed.periodogram <- function(x, I, b) {
  
  # Auxiliary useful variables
  n <- length(I)
  n.tilda <- n/2
  k <- seq(from = -n.tilda, to = n - 1)
  I.aux <- c(I[n.tilda:1], 0,  I[1:(n - 1)])
  freq.aux <- 2*pi*k/n
  
  # Kernel coefficients
  temp <- gaussian.kernel(freq.aux - x, b, scaled = TRUE)
  
  # Evaluating smoothed periodogram
  numerator <- temp %*% I.aux
  denominator <- sum(temp)
  I.smoothed <- numerator/denominator
  
  return(I.smoothed)
}

# Function: Smoothed periodogram
smoothed.periodogram <- function(freq, I, b) {
  out <- c()
  for(i in 1:length(freq)) {
    out[i] <- single.smoothed.periodogram(freq[i], I, b)
  }
  return(out)
}


### Function: unbiased estimate for Residual Sum of Squares (RSS)
unbiased.RSS <- function(b, I) {
  
  # n, \tilde{n}
  n <- length(I)
  n.tilda <- n/2
  # Frequencies
  k <- 1:n.tilda
  freq.rad <- 2*pi*k/n
  
  # Smoothed.periodogram
  I.smooth <- smoothed.periodogram(freq.rad, I, b = b)
  
  # Auxiliary variables
  k.aux <- seq(from = -n.tilda, to = n - 1)
  freq.aux <- 2*pi*k.aux/n
  
  # Constant for bias
  K.0 <- gaussian.kernel(0, b = b, scaled = T)
  norm.const <- sum(gaussian.kernel(freq.aux, b = b, scaled = T))
  W.b <- K.0/norm.const
  
  # RSS
  RSS <- sum((I[1:n.tilda] - I.smooth)^2)
  # Bias
  bias <- ((1 - 2*W.b)/2) * sum(I[1:n.tilda]^2)
  # Unbiased estimate
  out <- RSS - bias
  
  return(out)
}


# Function: find c that minimise the unbised RSS
find.c <- function(I) {
  
  # Length time series
  n <- length(I)
  # Grid
  c <- (1:1e3)/1e3
  b <- c * (n^(-1/5))
  
  # Risk values
  risk.values <- c()
  for(i in 1:length(c)) {
    risk.values[i] <- unbiased.RSS(b = b[i], I)
  }
  
  # Values of C that minimises the risk value
  out <- c[which(risk.values == min(risk.values))]
  
  return(out)
}








