setwd("C:/Users/Beniamino/Desktop/Mini_Project_1/Code/Dynamic_Harmonic _Regression")
source("Spectrums.R")


#### Testing Spectrums

# Frequencies:
freq <- seq(from = 0, to = 0.5, len = 1000)
freq <- freq[-1]; freq # Not considering frequency 0.


# Creating spectrums for different value of NVR (sigma2.eta/sigma.2)

spectrums.IRW <- list()
spectrums.RW <- list()

for(i in 1:6) {
  spectrums.IRW[[i]] <- normalise(log(pseudo.spectrum.GRW(freq, sigma2.eta = 5 * (10 ^ (-i)),
                                        sigma2 = 5, type = "IRW")))
}

for(i in 1:6) {
  spectrums.RW[[i]] <- normalise(log(pseudo.spectrum.GRW(freq, sigma2.eta = 5 * (10 ^ (-i)),
                                            sigma2 = 5, type = "RW")))
}


# Plotting normalised (log) spectrums RW for different value of NVR
par(mfrow = c(1, 2))

plot(freq, spectrums.RW[[1]], xlim = c(0, 0.2),
     type = "l", col = "red", lwd = 2, 
     xlab = "Frequency", ylab = "Power",
     main = "RW")
for(i in 2:5) {
  lines(freq, spectrums.RW[[i]], type = "l", lty = 2)
}
lines(freq, spectrums.RW[[6]], type = "l", col = "red")

# Plotting normalised (log)spectrums IRW for different value of NVR

plot(freq, spectrums.IRW[[1]], xlim = c(0, 0.2),
     type = "l", col = "red", lwd = 2,
     xlab = "Frequency", ylab = "Power",
     main = "IRW")
for(i in 2:5) {
  lines(freq, spectrums.IRW[[i]], type = "l", lty = 2)
}
lines(freq, spectrums.IRW[[6]], type = "l", col = "red")






####### Figure.3, pseudo-spectra DHR model (AirPassengers Data)

library(astsa)
data("AirPassengers")

# Variances: Obs variance, Trend variance, harmonic component variances
sigma <- c(1, 1/10, rep(1/10, 5)) 

# Trend frequency, Harmonic components frequencies
omega.j <- c(0, 1/12, 1/6, 1/4, 1/3, 1/(2.4))


# Finding AR(p) to approximate spectrum of AirPassanger
performance.ar.p(AirPassengers, p = 20, plot = TRUE, legend = F)

# Fitting AR(14) to AirPassangers
p <- 14
fit <- arima(AirPassengers, order = c(p, 0, 0), method="ML")
coef <- fit$coef[1:p]

# Spectrum AR(14)
ar.fit <- arma.spec(ar = coef, n.freq = length(freq) + 1,
                    var.noise = fit$sigma2)

# Plotting log(spectrum)
par(mfrow = c(1, 1))
plot(ar.fit$freq, log(ar.fit$spec), type = "l", col = "red", lwd = 3,
     xlab = "Frequency", ylab = "Power Spectrum", ylim = c(0, 30))

# Adding lines of the pseudo spectrum of IRW trend, and RW 
# harmonic components.
lines(freq, log(sigma[2] * S(freq, omega.j[1], type = "IRW")/(2*pi)),
      lty = 2, lwd = 2)
for(i in 2:6) {
  lines(freq, log(sigma[i + 1] * S(freq, omega.j[i], type = "IRW")/(2*pi)),
        lty = 2, col = "blue")
  }


