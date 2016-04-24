setwd("C:/Users/Beniamino/Desktop/Mini_Project_1/Code/On_Boostrapping_Kernel_Spectral_Estimates")
source("Kernel_Smoothing_Spectrum.R")
source("Bootstrapping_Spectral_Estimates.R")

library(astsa)
data("sunspotz")


########     KERNEL SMOOTHING SPECTRUM    ########



# Testing on Sunspotz Data (considering n even)
x <- as.vector(sunspotz);
x <- x[-1]
n <- length(x)

# Frequency (radiant)
n.tilda <- n/2
k <- 0:(n.tilda - 1)
freq.rad <- 2*pi*k/n

# Periodogram, and setting I(0) = 0
I <- periodogram(x)
I[1] <- 0

# Smoothed periodogram 
# Notice: to evaluate smoothed.periodogram we need to use freq.rad
I.smooth <- smoothed.periodogram(freq.rad, I, b = 0.14)

# Plotting periodgram and smoothed periodogram
plot(freq.rad, I[1:(n.tilda)], type = "l")
lines(freq.rad, I.smooth[1:(n.tilda)], type = "l", 
      col = "red", lwd = 2)

# We want to find the c, for which the unbised RSS(b) is minimised.
c <- find.c(I)
c;

# Therefore the optimal b is given by:
b <- c * (n ^ (-1/5)); b
I.smooth <- smoothed.periodogram(freq.rad, I, b = b)

# The plotted results:
plot(freq.rad, I[1:(n.tilda)], type = "l", ylab = "Power spectrum",
     xlab = "Frequency")
lines(freq.rad, I.smooth[1:(n.tilda)], type = "l", 
      col = "red", lwd = 2)

# Plotting the unbiased.RSS over the grid:
# c <- (1:1e3)/1e3
# b <- c * (n ^ (-1/5))

# Getting the values of unbiased.RSS
# values.RSS <- c()
# for(i in 1:1000) {
#  values.RSS[i] <- unbiased.RSS(b[i], I)
# }

# plot(c, values.RSS, type = "l", main = "Choosing c", 
#     ylab = "R(b)", xlab  = "c", col = "blue", lwd = 2)






########     SPECTRAL RESAMPLING (SR)   ########

# Applying SR methodology
SR.sample <- Spectral.Resampling(data = x, c = c, R = 100)
# Getting average, and confindence intervals for the spectrum
SR <- Spectral.Resampling.CI(SR.sample)
# Getting distributions of the SR peaks (n.peaks = 3)
SR.distr <- SR.distribution(SR.sample, n.peaks = 3)

# Plotting results
plot.SR.CI <- SR$plot; plot.SR.CI
plot.SR.distr <- SR.distr$plot; plot.SR.distr


# We need to think about the histogram for the period
hist(2*pi/c(SR.distr$counts), breaks = 50, col = "pink")

mean(SR.distr$counts[1, ])

length(SR$spec)
period <- 0:(n.tilda -1)
length(grid)
plot(period, SR$spec, type = "l")

2*pi/quantile(SR.distr$counts[3, ], c(.95, .05))
