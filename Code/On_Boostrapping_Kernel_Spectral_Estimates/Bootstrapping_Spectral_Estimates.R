library(ggplot2)
library(pracma)
source("Kernel_Smoothing_Spectrum.R")

### Function: single estimate for bootstrap spectrum
single.f.star <- function(x, I.smooth, residuals, b) {
  
  # Setting length parameters
  n <- length(I.smooth) * 2
  n.tilda <- length(I.smooth)
  
  # Parameters
  k <- seq(from = -n.tilda, to = n - 1)
  freq.aux <- 2*pi*k/n
  
  # Auxiliary objects
  residuals.aux <- c(residuals[n.tilda:1], 0, residuals[1:n.tilda],
                     residuals[(n.tilda - 1) : 1])
  I.smooth.aux <- c(I.smooth[n.tilda:1], 0, I.smooth[1:n.tilda],
                    I.smooth[(n.tilda - 1):1])
  
  # Bootstrap periodogram
  I.star <- I.smooth.aux * residuals.aux
  
  # Kernel coefficients
  temp <- gaussian.kernel(freq.aux - x, b, scaled = TRUE)
  
  # Evaluating bootstrap spectral estimates
  numerator <- temp %*% I.star
  denominator <- sum(temp)
  I.smoothed.BR <- numerator/denominator
  
  return(I.smoothed.BR)
}

### Function: vector estimate for bootstrap periodogram
f.star <- function(freq, I.smooth, residuals, b) {
  out <- c()
  for(i in 1:length(freq)) {
    out[i] <- single.f.star(freq[i], I.smooth, residuals, b)
  }
  return(out)
}


### Function: Spectral Resampling algorithm
Spectral.Resampling <- function(data, c, R) {
  
  # Setting length parameters 
  n <- length(data)
  n.tilda <- n/2
  k <- 0:(n.tilda - 1)
  freq.rad <- 2*pi*k/n
  
  # Setting bandwiths
  b.start <- c * (n ^ (-1/4))
  b.intermediate <- c * (n ^ (-1/6))
  b.final <- c * (n ^ (-1/5))
  
  # Output boostrap spectral estimate
  I.smooth.BS <- matrix(NA, nrow = n.tilda, ncol = R)
  
  for(r in 1:R) {
    
    # STEP 1:
    I <- periodogram(data)
    I[1] <- 0
    I.smooth <- smoothed.periodogram(freq.rad, I, b = b.start)
    residuals <- I[1:n.tilda]/I.smooth[1:n.tilda]
    residuals.norm <- residuals/(sum(residuals)/n.tilda)

    
    # STEP 2:
    residuals.BS <- sample(residuals.norm, size = n.tilda,
                           replace = TRUE, prob = rep(1/n.tilda, n.tilda))

    I.smooth.intermediate <- smoothed.periodogram(freq.rad, I, b = b.intermediate)
    
    # Final boostrap estimate: 
    I.smooth.BS[, r] <- f.star(freq.rad, I.smooth.intermediate,
                               residuals = residuals.BS, b = b.final)
  }
  class(I.smooth.BS) <- "SR"
  return(I.smooth.BS)
}




### Function: Get Average, and Confindence Intervals of SR
###           and an optional plot.

Spectral.Resampling.CI <- function(SR.sample) {
  
  # Getting dimensions and setting frequencies
  n.tilda <- nrow(SR.sample)
  n <- n.tilda*2
  k <- 0:(n.tilda - 1)
  freq.rad <- (2*pi*k)/ n 
  
  # Checking input
  if(class(SR.sample) != "SR") 
    stop("input not of class SR")
  
  # Creating confidence intervals for spectrum estimate 
  conf.int <- matrix(NA, nrow = n.tilda, ncol = 2)
  for(i in 1:n.tilda) {
    conf.int[i, ] <- quantile(SR.sample[i, ], probs = c(.05, .95))
  }
  # Average
  mean.spectrum.SR <- rowMeans(SR.sample)
  
  # Auxiliary data.frame for plotting
  dat <- data.frame(freq = freq.rad, spec = mean.spectrum.SR, 
                    ci.up = conf.int[, 2], ci.low = conf.int[, 1])
  # Plotting mean and C.I
  p <- ggplot() +
    geom_line(data = dat, 
              aes(x = freq.rad, y = spec, linetype = "c"),
              size = 1.0, colour = "red") +
    geom_line(data = dat, 
              aes(x = freq.rad, y = ci.up), 
              size = 0.6, colour = "black", linetype = "dashed") + 
    geom_line(data = dat, 
              aes(x = freq.rad, y = ci.low), 
              size = 0.6, colour = "black", linetype = "dashed") +
    xlab("Frequency") + ylab("Power Spectrum") +
    ggtitle("Spectral Resampling (SR)") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 19)) +
    scale_x_continuous(minor_breaks = seq(0, pi, 0.1))

    
  return(list(spec = mean.spectrum.SR, 
              ci.up = conf.int[, 2], ci.low = conf.int[, 1], 
              plot = p))
}




### Function: distribution of the peaks of SR sample.
#             This function also provides an histogram
SR.distribution <- function(SR.sample, n.peaks) {
  
  # Checking input
  if(class(SR.sample) != "SR") 
    stop("input not of class SR")
  
  # This library implement the function findpeaks
  library(pracma)
  
  # Setting dimensions and frequency parameters
  R <- dim(SR.sample)[2]
  n.tilda <- dim(SR.sample)[1]
  n <- n.tilda * 2
  k <- 0:(n.tilda - 1)
  freq.rad <- 2*pi*k/n
  
  # Largest peaks in each bootstrap sample
  counts <- matrix(NA, nrow = n.peaks, ncol = R)
  
  # Getting largest peaks for each bootstrap sample
  for(n in 1:n.peaks) {
    for(r in 1:R) {
      counts[n, r] <- freq.rad[findpeaks(SR.sample[, r])[n, 2]]
    }
  }
  
  # Plotting histogram
  dat <- data.frame(counts = c(counts))
  p <- ggplot(data = dat,
              aes(dat$counts)) +
    geom_histogram(fill = "blue", col = "black", breaks = seq(0, pi, by = .04)) +
    coord_cartesian(xlim = c(0, pi)) +
    scale_x_continuous(minor_breaks = seq(0, pi, 0.1)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm")) +
    xlab("Frequency") + ylab("Power Spectrum")
  
  return(list(counts = counts, plot = p))
}
