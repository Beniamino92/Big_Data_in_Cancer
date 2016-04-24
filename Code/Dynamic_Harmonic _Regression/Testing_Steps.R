# Testing steps of Complete DHR Estimation Algorithm 

setwd("C:/Users/Beniamino/Desktop/Mini_Project_1/Code/Dynamic_Harmonic _Regression")


###### Testing Step 1 and Step 2

source("Spectrums.R")
source("Steps_DHR_Estimation_Algorithm.R")



# Getting data
data.AP <- AirPassengers
T <- length(data.AP)

# Harmonics
omega.j <- c(0, 1/12, 1/6, 1/4, 1/3, 1/(2.4))
omega.j.AP <- c(0, 1/12, 1/6, 1/4, 1/3, 1/(2.4))

# IRW Trend, RW harmonic components
structures.AP <- c("IRW", rep("RW", 5))

# Which p for the AR(p)?
performance.ar.p(data.AP, p = 20, plot = TRUE, legend = FALSE)

# Getting spectrum and residuals variance from AR(P) fitted to AirPassengers
AP.AR14 <- estimate.AR(data = data.AP, p = 14, n.freq = length(data.AP))
freq <- AP.AR14$frequency
spec <- AP.AR14$spectrum
sigma2 <- AP.AR14$sigma2

# Estimating NVR with Step2 (LLS)
NVR.LSS.estimate <- NVR.LinearLS(empirical.spectrum = spec, 
                                 sigma2 = sigma2,
                                 frequency = freq, 
                                 omega.j = omega.j.AP, 
                                 structures = structures.AP)

# Sigma estimate
sigma.est <- c(sigma2, sigma2 * NVR.LSS.estimate)

# Plotting AR(14) spectrum and DHR spectrum with estimated sigma (LLS)
plot(freq, log(spec), type = "l", col = "red", lwd = 2)
lines(freq, log(fullDHR.spectrum(freq,  sigma.est, omega.j.AP)),
      lwd = 2, col = "green")




###### Testing Step 3

# Optimising logarithmic error: (THERE ARE WARNINGS)
optimising.NVR <- optim(NVR.LSS.estimate, E, 
                        freq = freq, empirical.spectrum = spec,
                        omega.j = omega.j.AP, structures = structures.AP,
                        sigma2 = sigma2, method = "Nelder-Mead")

# Final NVR estimate
NVR.nonlinear.estimate <- optimising.NVR$par

# Final sigma^2 estimate
final.sigma.est <- c(sigma2, sigma2 * NVR.nonlinear.estimate)

# Plotting AR(14) spectrum and DHR spectrum with estimated sigma, 
# by using log loss function

plot(freq, log(spec), type = "l", col = "red", lwd = 2)
lines(freq, log(fullDHR.spectrum(freq,  final.sigma.est, omega.j.AP)),
      lwd = 2, col = "green")






# Testing Step 4
