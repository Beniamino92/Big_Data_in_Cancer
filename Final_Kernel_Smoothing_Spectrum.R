library(ggplot2)
library(pracma)

# Periodogram
periodogram <- function(x) {
  I <- spec.pgram(x, detrend = FALSE, taper = 0, plot = FALSE)$spec
  freq <- spec.pgram(x, detrend = FALSE, taper = 0, plot = FALSE)$freq
  freq.rad <- freq*2*pi
  return(list(I = I, freq = freq, freq.rad = freq.rad))
}

# Function: scaled gaussian kernel
gauss.kernel <- function(omega, b) {
  (1/b)* (1/sqrt(2*pi)) * exp(-((omega/b)^2)/2)
}



# Function: evaluate the risk function
risk.function <- function(freq.rad, I, b) {
  
  # Parameters and auxiliary objects
  n.tilda <- length(I)
  freq.aux <- c(freq.rad[n.tilda:1], 0, freq.rad[1:n.tilda], freq.rad[(n.tilda - 1):1])
  W.b <- gauss.kernel(0, b)/sum(gauss.kernel(freq.aux, b))
  
  # Smoothed periodogram
  I.smooth <- smooth.periodogram(freq.rad, I, b)
  
  # Get estimate, bias, and unbiased estimator
  bias <- ((1 - 2*W.b)/2) * sum(c(0, I[1:(n.tilda - 1)])^2)
  est <- sum((c(0, I[1:(n.tilda - 1)]) - c(0, I.smooth[1:(n.tilda  - 1)]))^2)
  out <- est - bias
  
  return(out)
}

# Function: Smoothed Periodogram 
smooth.periodogram <- function(freq.rad, I, b) {
  
  n.tilda <- length(I)
  n <- n.tilda * 2
  I.smooth <- c()
  I.aux <- c(I[n.tilda:1], 0, I[1:n.tilda], I[(n.tilda - 1):1])
  freq.aux <- c(freq.rad[n.tilda:1], 0, freq.rad[1:n.tilda], freq.rad[(n.tilda - 1):1])
  
  for(k in 1:n.tilda) {
    kern.aux <- gauss.kernel(freq.aux - freq.rad[k], b = b)
    norm.const <- sum(gauss.kernel(freq.aux - freq.rad[k], b = b))
    I.smooth[k] <- (kern.aux %*% I.aux)/norm.const
  }
  return(I.smooth)
}


# Function: find c
get.c <- function(freq.rad, I) {
  
  # Setting parameters
  n.tilda <- length(I)
  n <- n.tilda * 2
  # Grid
  c.val <- (1:1e3)/1e3
  b <- c.val*(n^(-1/5))
  
  # Finding c
  risk.val <- c()
  for(j in 1:length(b)) {
    risk.val[j] <- risk.function(freq.rad, I, b = b[j])
  }
  out <- c.val[which(risk.val == min(risk.val))]
  
  return(out)
}



# Function: boostrap spectral estimate
smooth.periodogram.bootstrap <- function(freq.rad, I.smooth.int, 
                                         sample.residuals, b.final) {
  
  # Parameters and auxiliary objects
  n.tilda <- length(freq.rad)
  
  residuals.aux <- c(sample.residuals[n.tilda:1], 0, sample.residuals[1:n.tilda], 
                     sample.residuals[(n.tilda - 1):1])
  I.smooth.aux <- c(I.smooth.int[n.tilda:1], 0, I.smooth.int[1:n.tilda], 
                    I.smooth.int[(n.tilda - 1):1])
  freq.aux <- c(freq.rad[n.tilda:1], 0, freq.rad[1:n.tilda], 
                freq.rad[(n.tilda - 1):1])
  I.star <- I.smooth.aux * residuals.aux
  
  out <- c()
  for(k in 1:n.tilda) {
    kern.aux <- gauss.kernel(freq.aux - freq.rad[k], b = b.final)
    norm.const <- sum(gauss.kernel(freq.aux - freq.rad[k], b = b.final))
    out[k] <- (I.star %*% kern.aux)/norm.const
  }
  return(out)
}

# Function: Spectral Resampling method
Spectral.Resampling <- function(freq.rad, I, c.val, R) {
  
  n.tilda <- length(I)
  n <- n.tilda * 2
  I.smooth.boot <- matrix(NA, ncol = R, nrow = n.tilda)
  
  b.star <- c.val*(n^(-1/4))
  b.int <- c.val*(n^(-1/6))
  b.final <- c.val*(n^(-1/5))
  
  # Step 1:
  I.smooth <- smooth.periodogram(freq.rad, I, b = b.star)
  residuals <- I/I.smooth
  residuals.scaled <- residuals/(sum(residuals)/n.tilda)
  
  # Step 2:
  for(r in 1:R) {
    sample.residuals <- sample(residuals.scaled, size = n.tilda, replace = TRUE)
    I.smooth.int <- smooth.periodogram(freq.rad, I, b = b.int)
    I.smooth.boot[, r] <- smooth.periodogram.bootstrap(freq.rad, I.smooth.int, 
                                                       sample.residuals, b.final)
  }
  return(I.smooth.boot)
}


### Function: Get Average, and Confindence Intervals of SR
###           and an optional plot.

Spectral.Resampling.CI <- function(freq, SR.sample) {
  
  # Getting dimensions and setting frequencies
  n.tilda <- length(freq)
  n <- n.tilda*2
  
  
  # Creating confidence intervals for spectrum estimate 
  conf.int <- matrix(NA, nrow = n.tilda, ncol = 2)
  for(i in 1:n.tilda) {
    conf.int[i, ] <- quantile(SR.sample[i, ], probs = c(.05, .95))
  }
  # Average
  mean.spectrum.SR <- rowMeans(SR.sample)
  
  # Auxiliary data.frame for plotting
  dat <- data.frame(freq = freq, spec = mean.spectrum.SR, 
                    ci.up = conf.int[, 2], ci.low = conf.int[, 1])
  # Plotting mean and C.I
  p <- ggplot() +
    geom_line(data = dat, 
              aes(x = freq, y = spec, linetype = "c"),
              size = 1.0, colour = "red") +
    geom_line(data = dat, 
              aes(x = freq, y = ci.up), 
              size = 0.6, colour = "black", linetype = "dashed") + 
    geom_line(data = dat, 
              aes(x = freq, y = ci.low), 
              size = 0.6, colour = "black", linetype = "dashed") +
    xlab("Frequency") + ylab("Power Spectrum") +
    theme(legend.position = "none", 
          plot.title = element_text(size = 19)) +
    scale_x_continuous(minor_breaks = seq(min(freq), max(freq), 0.05))
  
  
  return(list(spec = mean.spectrum.SR, 
              ci.up = conf.int[, 2], ci.low = conf.int[, 1], 
              plot = p))
}

