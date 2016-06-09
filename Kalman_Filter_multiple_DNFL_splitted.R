

##### Sourcing functions #####
setwd("C:/Users/Beniamino/Desktop/Mini_Project_1/Code/Dynamic_Harmonic_Regression")
source("Kalman_Filter.R")
source("Spectrums.R")
source("Steps_DHR_Estimation_Algorithm.R")
setwd("C:/Users/Beniamino/Desktop/Mini_Project_1/Code/Inference_on_periodicity_of_Circadian_Time_Series")
source("DNFL_model.R")
setwd("C:/Users/Beniamino/Desktop/Mini_Project_1/Code/On_Boostrapping_Kernel_Spectral_Estimates")
source("Final_Kernel_Smoothing_Spectrum.R")





############################### 1 ######################################





##### Getting DNFL model 
model <- get.DNFL2(N = 400, omega = 200, mu = 6, alpha = 10, 
                   mu2 = 2, alpha2 = 15, 
                   SDE = TRUE)$M
N <- length(model)
signal.1 <- model[1:200]
signal.2 <- model[201:400]
model[1:200] <- model[1:200] + rnorm(200, 0, 50)
model[201:400] <- model[201:400] + rnorm(200, 0, 25)
plot(1:N, model, type = "o", col = "red", pch = 20)

dat <- data.frame(mRNA = model)

plot.series <- ggplot() + 
  geom_point(data = dat, aes(x = 1:N, y = mRNA), colour = "red", size = 2) +
  geom_line(data = dat, aes(x = 1:N, y = mRNA, linetype = "o"), colour = "red") +
  theme(panel.background = element_rect(colour = 'black'), 
        axis.title.x = element_blank(), legend.position = "none") 
plot.series

####################### Entire Series ########################

##### Periodogram
I <- periodogram(model)$I
freq <- periodogram(model)$freq
n.tilda <- length(freq)

##### Selecting p
performance <- performance.ar.p(model, p = 60, 
                                plot = TRUE, legend = FALSE)
p <- which(performance$AIC == min(performance$AIC))
p <- 22

##### Fitting AR(p)
fit <- arima(model, order = c(p, 0, 0), method = "ML")
coef <- fit$coef[1:p]
sigma2 <- fit$sigma2


##### Deciding how many and which frequencies to use

# Auxiliary AR(p) object
temp.ar.fit <- arma.spec(ar = coef, var.noise = sigma2, n.freq = n.tilda)
# Auxiliary  Frequencies
temp.ar.freq <- temp.ar.fit$freq; temp.ar.freq <- temp.ar.freq[-1]
# Auxiliary  Spectrum
temp.ar.spec <- temp.ar.fit$spec; temp.ar.spec <- temp.ar.spec[-1]
plot(temp.ar.freq, log(temp.ar.spec), type = "l", col = "red")

# Peaks 
peaks <- findpeaks(as.vector(temp.ar.spec)); peaks
driving.freq <- c()
for(j in 1:nrow(peaks)) {
  driving.freq[j] <- temp.ar.freq[peaks[j, 2]]
}
driving.freq <- c(freq[17], freq[34])
## AR(p) Object 
ar.fit <- arma.spec(ar = coef, var.noise = sigma2, n.freq = 500)
## Frequencies (no omega_0)
ar.freq <- ar.fit$freq; ar.freq <- ar.freq[-1]
## Spectrum (no omega_0)
ar.spec <- ar.fit$spec; ar.spec <- ar.spec[-1]
plot(ar.freq, log(ar.spec), type = "l")


### Setting foundmental frequency and other harmonics
omega.model <- c(0, driving.freq); omega.model
### Modelling the evolution of the coefficient process 
structures.model <- c("IRW", rep("RW", length(driving.freq)))


###### Estimation of the NVR coefficients

# Linear Least Square approximation
NVR.LSS <- NVR.Linear.LS(ar.spec, sigma2, ar.freq, omega.model, 
                         structures.model)
# Relative approximation of the variances
sigma.LSS <- c(sigma2, sigma2 * NVR.LSS)
# Getting full DHR spectrum
middle.DHRspectrum <- fullDHR.spectrum(ar.freq, sigma.LSS, omega.model)
plot(ar.freq, log(middle.DHRspectrum), type = "l", col = "green")
lines(ar.freq, log(ar.spec), type = "l", col = "red") 



# Non linear approximation
NVR <- optim(NVR.LSS, E, 
             freq = ar.freq, empirical.spectrum = ar.spec, 
             omega.j = omega.model, 
             structures = structures.model, 
             sigma2 = sigma2, 
             method = "Nelder-Mead")$par
# Relative approximation of the variances
sigma.est <- c(sigma2, sigma2 * NVR)
# Getting full DHR spectrum
full.DHRspectrum <- fullDHR.spectrum(ar.freq, sigma.est, omega.model)

plot(ar.freq, log(full.DHRspectrum), col = "green", type = "l")
lines(ar.freq, log(ar.spec), type = "l", col = "red")


# Non linear approximation
for(j in 1:1e2) {
  NVR <- optim(NVR, E, 
               freq = ar.freq, empirical.spectrum = ar.spec, 
               omega.j = omega.model, 
               structures = structures.model, 
               sigma2 = sigma2, 
               method = "Nelder-Mead")$par
}

# Relative approximation of the variances
sigma.est <- c(sigma2, sigma2 * NVR)
# Getting full DHR spectrum
full.DHRspectrum <- fullDHR.spectrum(ar.freq, sigma.est, omega.model)

plot(ar.freq, log(full.DHRspectrum), col = "green", type = "l")
lines(ar.freq, log(ar.spec), type = "l", col = "red")


############# Filtering ###########

# Getting F matrix
F.IRW <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
F.RW <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
F <- bdiag(c(list(F.IRW), replicate(length(omega.model) - 1, 
                                    F.RW, simplify = FALSE)))
# Getting G matrix
G.IRW <- matrix(c(0, 0, 0, 1), nrow = 2, ncol = 2)
G.RW <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
G <- bdiag(c(list(G.IRW), replicate(length(omega.model) - 1,
                                    G.RW, simplify = FALSE)))
# Getting Q matrix
Q <- diag(rep(sigma2 * NVR, 1, each = 2))
# Getting H list
H <- get.H(omega.model, N)


### Kalman Filter
filtered.process <- Kalman.Filter(y = model, F = F, H = H, 
                                  G = G, Q = Q, sigma2 = sigma2)
### Kalman Smoothing
smoothed.process <- Kalman.Smoother(filtered.process)



######   Visualisation   #####
fitted <- fitted.DHR(H, filtered.process$states)
smooth.fitted <- fitted.DHR(H, smoothed.process$smoothed.states) 

plot(1:N, model, col = "grey", pch = 20)
points(1:N, model, type = "o", col = "grey", pch = 20)
lines(1:N, fitted, type = "o", col = "green")
lines(1:N, smooth.fitted, type = "o", col = "purple", pch = 20)

lines(1:200, signal.1, lty = 3, col = "red", lwd = 2)
lines(201:400, signal.2, lty = 3, col = "red", lwd = 2)

residuals <- model - smooth.fitted
residuals.norm <- residuals/(sd(residuals))
plot.ts(residuals.norm, type = "l")


########  Forecasting 24 hours
t <- 100
K <- 48

forecasting.process <- Kalman.Forecasting(filtered.process, t, K)
forecast <- forecasting.process$forecast
hwid <- -qnorm((1 - 0.95)/2) * sqrt(forecasting.process$forecast.variance)
plot(1:t, model[1:t], type = "o", col = "grey", pch = 20, 
     xlim = c(1, (t + K + 1)), ylim = c(min(forecast - hwid), max(c(forecast + hwid, model))))
lines(1:t, fitted[1:t], type = "o", col = "green")
lines(1:t, smooth.fitted[1:t], type = "o", col = "purple", pch = 20)
lines((t + 1):(t + K), forecast, type = "o", col = "red", pch = 20)
lines((t + 1):(t + K), forecast + hwid, lty = 3)
lines((t + 1):(t + K), forecast - hwid, lty = 3)
points((t + 1):(t + K), model[(t + 1):(t + K)], col = "grey", type = "o")

dat.forecast <- data.frame(forecast = forecast, upr = forecast + hwid,
                           lwr = forecast - hwid, 
                           grid = (t + 1):(t + K),
                           model = model[(t + 1):(t + K)])
head(dat.forecast)
           
my.plot <- ggplot() + 
  geom_point(data = dat.forecast, aes(x = grid, y = model), colour = "grey", size = 3) +
  geom_line(data = dat.forecast, aes(x = grid, y = model, linetype = "o"), colour = "grey") +
  geom_point(data = dat.forecast, aes(x = grid, y = forecast), colour = "blue", size = 3) +
  geom_line(data = dat.forecast, aes(x = grid, y = forecast, linetype = "dotdashed"), colour = "blue") + 
  geom_ribbon(data = dat.forecast, aes(x = grid, ymin = lwr, ymax = upr), 
              alpha = 0.15, fill = "blue") + 
  ylab("mRNA") + 
  theme(panel.background = element_rect(colour = 'black'), 
        axis.title.x = element_blank(), legend.position = "none") 
  
my.plot

grid.2 <- (t + 1):(t + K) - 200
signal.2[grid.2]




plot((t + 1):(t + K),  model[(t + 1):(t + K)], 
     type = "o", col = "grey", pch = 20)
lines((t + 1):(t + K), forecast, type = "o", col = "red")
(t + 1):(t + K)

lines((t + 1):(t + K), signal.2[(t - 200 + 1):(t + K - 200 + 1)])

(t - 200 + 1)
length((t + 1):(t + K))
length(forecast)
length(forecast + hwid)
length(forecast - hwid)
length(model[(t + 1):(t + K)])
plot.ts(signal)

##### Plotting on ggplot2

library(ggplot2)
dat <- data.frame(model = model, fitted = fitted, smooth.fitted = smooth.fitted)
dat.forecast <- data.frame(forecast = forecast, upr = forecast + hwid,
                           lwr = forecast - hwid)
partial.dat <- data.frame(model = model[1:t], fitted = fitted[1:t])

# Plot: observed, filtered, and smoothed.
p1 <- ggplot() + geom_point(data = dat,
                            aes(x = 1:N, y = model), colour = "grey", size = 2) + 
  geom_line(data = dat, aes(x = 1:N, y = model, linetype = "o"), colour = "grey") +
  geom_point(data = dat, aes(x = 1:N, y = fitted), shape = 1,  colour = "green", size = 2.7) +
  geom_line(data = dat, aes(x = 1:N, y = fitted, linetype = "o"), colour = "green") +
  # geom_point(data = dat, aes(x = 1:N, y = smooth.fitted), colour = "purple", size = 2.7) +
  # geom_line(data = dat, aes(x = 1:N, y = smooth.fitted, linetype = "o"), colour = "purple") +
  theme(panel.background = element_rect(colour = 'black'), 
        legend.position = "none")
p1

# Plot: observed, smoothed and forecasting of 24 hours
p2 <- ggplot() + 
  geom_point(data = partial.dat,
             aes(x = 1:t, y = model), colour = "grey", size = 3) +
  geom_line(data = partial.dat, aes(x = 1:t, y = model, linetype = "o"), colour = "grey") +
  geom_point(data = partial.dat, aes(x = 1:t, y = fitted), colour = "green", size = 2) +
  xlim(c(0, t + K + 1)) +
  geom_line(data = partial.dat, aes(x = 1:t, y = fitted, linetype = "o"), colour = "green") +
  geom_ribbon(data = dat.forecast, aes(x = (t + 1):(t + K), ymin = lwr, ymax = upr), 
              alpha = 0.15, fill = "blue") + 
  geom_line(data = dat.forecast, aes(x = (t + 1):(t + K), y = forecast, linetype = "o"), 
            colour = "blue") +
  geom_point(data = dat.forecast, aes(x = (t + 1):(t + K), y = forecast), colour = "blue", size = 2) +
  ylab("mRNA") + 
  theme(panel.background = element_rect(colour = 'black'),
        axis.title.x = element_blank(), legend.position = "none")
p2






# mRNA before Chemo and during Chemo
preChemo <- model[1:200]
Chemo <- model[201:400]
plot(1:200, preChemo, type = "o", col = "red", pch = 20)
plot(201:400, Chemo, type = "o", col = "red", pch = 20)



###################################
#####      preChemo           #####
###################################



##### Periodogram
preChemo.I <- periodogram(preChemo)$I
preChemo.freq <- periodogram(preChemo)$freq
preChemo.n.tilda <- length(preChemo.freq)

##### Selecting p
preChemo.performance <- performance.ar.p(preChemo, p = 30, 
                                plot = TRUE, legend = TRUE)
preChemo.p <- which(preChemo.performance$AIC == min(preChemo.performance$AIC))
# preChemo.p <- 15


##### Fitting AR(p)
preChemo.fit <- arima(preChemo, order = c(preChemo.p, 0, 0), method = "ML")
preChemo.coef <- preChemo.fit$coef[1:preChemo.p]
preChemo.sigma2 <- preChemo.fit$sigma2


##### Deciding how many and which frequencies to use

# Auxiliary AR(p) object
preChemo.temp.ar.fit <- arma.spec(ar = preChemo.coef, 
                                  var.noise = preChemo.sigma2, 
                                  n.freq = preChemo.n.tilda)
# Auxiliary  Frequencies
preChemo.temp.ar.freq <- preChemo.temp.ar.fit$freq; 
preChemo.temp.ar.freq <- preChemo.temp.ar.freq[-1]
# Auxiliary  Spectrum
preChemo.temp.ar.spec <- preChemo.temp.ar.fit$spec; 
preChemo.temp.ar.spec <- preChemo.temp.ar.spec[-1]
plot(preChemo.temp.ar.freq, log(preChemo.temp.ar.spec),
     type = "l")

# Peaks 
preChemo.peaks <- findpeaks(as.vector(preChemo.temp.ar.spec)); preChemo.peaks
preChemo.driving.freq <- c()
for(j in 1:nrow(preChemo.peaks)) {
  preChemo.driving.freq[j] <- preChemo.temp.ar.freq[preChemo.peaks[j, 2]]
}

#preChemo.driving.freq <- c(preChemo.temp.ar.freq[9], 
#                           preChemo.temp.ar.freq[18],
#                           preChemo.temp.ar.freq[34], 
#                           preChemo.temp.ar.freq[46], 
#                           preChemo.temp.ar.freq[60], 
#                           preChemo.temp.ar.freq[71])

## AR(p) Object 
preChemo.ar.fit <- arma.spec(ar = preChemo.coef, 
                             var.noise = preChemo.sigma2, 
                             n.freq = 500)
## Frequencies (no omega_0)
preChemo.ar.freq <- preChemo.ar.fit$freq; 
preChemo.ar.freq <- preChemo.ar.freq[-1]
## Spectrum (no omega_0)
preChemo.ar.spec <- preChemo.ar.fit$spec; 
preChemo.ar.spec <- preChemo.ar.spec[-1]
plot(preChemo.ar.freq, log(preChemo.ar.spec), type = "l")


### Setting foundmental frequency and other harmonics
preChemo.omega <- c(0, preChemo.driving.freq); preChemo.omega
### Modelling the evolution of the coefficient process 
preChemo.structures <- c("IRW", rep("RW", length(preChemo.driving.freq)))


###### Estimation of the NVR coefficients

# Linear Least Square approximation
preChemo.NVR.LSS <- NVR.Linear.LS(preChemo.ar.spec,
                                  preChemo.sigma2, 
                                  preChemo.ar.freq, preChemo.omega, 
                                  preChemo.structures)
# Relative approximation of the variances
preChemo.sigma.LSS <- c(preChemo.sigma2, preChemo.sigma2 * preChemo.NVR.LSS)
# Getting full DHR spectrum
preChemo.middle.DHRspectrum <- fullDHR.spectrum(preChemo.ar.freq, 
                                                preChemo.sigma.LSS, 
                                                preChemo.omega)
plot(preChemo.ar.freq, log(preChemo.middle.DHRspectrum), type = "l", col = "green")
lines(preChemo.ar.freq, log(preChemo.ar.spec), type = "l", col = "red") 


# Non linear approximation
preChemo.NVR <- optim(preChemo.NVR.LSS, E, 
             freq = preChemo.ar.freq, 
             empirical.spectrum = preChemo.ar.spec, 
             omega.j = preChemo.omega, 
             structures = preChemo.structures, 
             sigma2 = preChemo.sigma2, 
             method = "Nelder-Mead")$par
# Relative approximation of the variances
preChemo.sigma.est <- c(preChemo.sigma2, preChemo.sigma2 * preChemo.NVR)
# Getting full DHR spectrum
preChemo.full.DHRspectrum <- fullDHR.spectrum(preChemo.ar.freq,
                                              preChemo.sigma.est,
                                              preChemo.omega)

plot(preChemo.ar.freq, log(preChemo.full.DHRspectrum), col = "green", type = "l")
lines(preChemo.ar.freq, log(preChemo.ar.spec), type = "l", col = "red")


# Non linear approximation
preChemo.NVR <- optim(preChemo.NVR, E, 
                      freq = preChemo.ar.freq, 
                      empirical.spectrum = preChemo.ar.spec, 
                      omega.j = preChemo.omega, 
                      structures = preChemo.structures, 
                      sigma2 = preChemo.sigma2, 
                      method = "Nelder-Mead")$par
# Relative approximation of the variances
preChemo.sigma.est <- c(preChemo.sigma2, preChemo.sigma2 * preChemo.NVR)
# Getting full DHR spectrum
preChemo.full.DHRspectrum <- fullDHR.spectrum(preChemo.ar.freq,
                                              preChemo.sigma.est,
                                              preChemo.omega)

plot(preChemo.ar.freq, log(preChemo.full.DHRspectrum), col = "green", type = "l")
lines(preChemo.ar.freq, log(preChemo.ar.spec), type = "l", col = "red")






# Getting F matrix
F.IRW <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
F.RW <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
preChemo.F <- bdiag(c(list(F.IRW), replicate(length(preChemo.omega) - 1, 
                                    F.RW, simplify = FALSE)))
# Getting G matrix
G.IRW <- matrix(c(0, 0, 0, 1), nrow = 2, ncol = 2)
G.RW <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
preChemo.G <- bdiag(c(list(G.IRW), replicate(length(preChemo.omega) - 1,
                                    G.RW, simplify = FALSE)))
# Getting Q matrix
preChemo.Q <- diag(rep(preChemo.sigma2 * preChemo.NVR, 1, each = 2))
# Getting H list
preChemo.H <- get.H(preChemo.omega, length(preChemo))



###################################
#####         Chemo           #####
###################################

##### Periodogram
Chemo.I <- periodogram(Chemo)$I
Chemo.freq <- periodogram(Chemo)$freq
Chemo.n.tilda <- length(Chemo.freq)

##### Selecting p
Chemo.performance <- performance.ar.p(Chemo, p = 30, 
                                         plot = TRUE, legend = TRUE)
Chemo.p <- which(Chemo.performance$AIC == min(Chemo.performance$AIC))



##### Fitting AR(p)
Chemo.fit <- arima(Chemo, order = c(Chemo.p, 0, 0), method = "ML")
Chemo.coef <- Chemo.fit$coef[1:Chemo.p]
Chemo.sigma2 <- Chemo.fit$sigma2


##### Deciding how many and which frequencies to use

# Auxiliary AR(p) object
Chemo.temp.ar.fit <- arma.spec(ar = Chemo.coef, 
                                  var.noise = Chemo.sigma2, 
                                  n.freq = Chemo.n.tilda)
# Auxiliary  Frequencies
Chemo.temp.ar.freq <- Chemo.temp.ar.fit$freq; 
Chemo.temp.ar.freq <- Chemo.temp.ar.freq[-1]
# Auxiliary  Spectrum
Chemo.temp.ar.spec <- Chemo.temp.ar.fit$spec; 
Chemo.temp.ar.spec <- Chemo.temp.ar.spec[-1]
plot(Chemo.temp.ar.freq, log(Chemo.temp.ar.spec),
     type = "l")

# Peaks 
Chemo.peaks <- findpeaks(as.vector(Chemo.temp.ar.spec))
Chemo.driving.freq <- c()
for(j in 1:nrow(Chemo.peaks)) {
  Chemo.driving.freq[j] <- Chemo.temp.ar.freq[Chemo.peaks[j, 2]]
}

# Chemo.driving.freq <- c(Chemo.temp.ar.freq[17])
## AR(p) Object 
Chemo.ar.fit <- arma.spec(ar = Chemo.coef, 
                             var.noise = Chemo.sigma2, 
                             n.freq = 500)
## Frequencies (no omega_0)
Chemo.ar.freq <- Chemo.ar.fit$freq; 
Chemo.ar.freq <- Chemo.ar.freq[-1]
## Spectrum (no omega_0)
Chemo.ar.spec <- Chemo.ar.fit$spec; 
Chemo.ar.spec <- Chemo.ar.spec[-1]
plot(Chemo.ar.freq, log(Chemo.ar.spec), type = "l")


### Setting foundmental frequency and other harmonics
Chemo.omega <- c(0, Chemo.driving.freq); Chemo.omega
### Modelling the evolution of the coefficient process 
Chemo.structures <- c("IRW", rep("RW", length(Chemo.driving.freq)))


###### Estimation of the NVR coefficients

# Linear Least Square approximation
Chemo.NVR.LSS <- NVR.Linear.LS(Chemo.ar.spec,
                                  Chemo.sigma2, 
                                  Chemo.ar.freq, Chemo.omega, 
                                  Chemo.structures)
# Relative approximation of the variances
Chemo.sigma.LSS <- c(Chemo.sigma2, Chemo.sigma2 * Chemo.NVR.LSS)
# Getting full DHR spectrum
Chemo.middle.DHRspectrum <- fullDHR.spectrum(Chemo.ar.freq, 
                                                Chemo.sigma.LSS, 
                                                Chemo.omega)
plot(Chemo.ar.freq, log(Chemo.middle.DHRspectrum), type = "l", col = "green")
lines(Chemo.ar.freq, log(Chemo.ar.spec), type = "l", col = "red") 


# Non linear approximation
Chemo.NVR <- optim(Chemo.NVR.LSS, E, 
                      freq = Chemo.ar.freq, 
                      empirical.spectrum = Chemo.ar.spec, 
                      omega.j = Chemo.omega, 
                      structures = Chemo.structures, 
                      sigma2 = Chemo.sigma2, 
                      method = "Nelder-Mead")$par
# Relative approximation of the variances
Chemo.sigma.est <- c(Chemo.sigma2, Chemo.sigma2 * Chemo.NVR)
# Getting full DHR spectrum
Chemo.full.DHRspectrum <- fullDHR.spectrum(Chemo.ar.freq,
                                              Chemo.sigma.est,
                                              Chemo.omega)

plot(Chemo.ar.freq, log(Chemo.full.DHRspectrum), col = "green", type = "l")
lines(Chemo.ar.freq, log(Chemo.ar.spec), type = "l", col = "red")


# Non linear approximation
Chemo.NVR <- optim(Chemo.NVR, E, 
                   freq = Chemo.ar.freq, 
                   empirical.spectrum = Chemo.ar.spec, 
                   omega.j = Chemo.omega, 
                   structures = Chemo.structures, 
                   sigma2 = Chemo.sigma2, 
                   method = "Nelder-Mead")$par
# Relative approximation of the variances
Chemo.sigma.est <- c(Chemo.sigma2, Chemo.sigma2 * Chemo.NVR)
# Getting full DHR spectrum
Chemo.full.DHRspectrum <- fullDHR.spectrum(Chemo.ar.freq,
                                           Chemo.sigma.est,
                                           Chemo.omega)

plot(Chemo.ar.freq, log(Chemo.full.DHRspectrum), col = "green", type = "l")
lines(Chemo.ar.freq, log(Chemo.ar.spec), type = "l", col = "red")



# Getting F matrix
Chemo.F <- bdiag(c(list(F.IRW), replicate(length(Chemo.omega) - 1, 
                                             F.RW, simplify = FALSE)))
# Getting G matrix
Chemo.G <- bdiag(c(list(G.IRW), replicate(length(Chemo.omega) - 1,
                                             G.RW, simplify = FALSE)))
# Getting Q matrix
Chemo.Q <- diag(rep(Chemo.sigma2 * Chemo.NVR, 1, each = 2))
# Getting H list
Chemo.H <- get.H.2(Chemo.omega, 201:400)





###########################################################
############## Filtering preChemo and Chemo ###############
###########################################################

### preChemo: Kalman Filter
preChemo.filtered.process <- Kalman.Filter(y = preChemo, F = preChemo.F, 
                                           H = preChemo.H, 
                                  G = preChemo.G,
                                  Q = preChemo.Q, 
                                  sigma2 = preChemo.sigma2/10)
preChemo.fitted <- fitted.DHR(preChemo.H,
                              preChemo.filtered.process$states)
### preChemo: Kalman Smoothing
# preChemo.smoothed.process <- Kalman.Smoother(preChemo.filtered.process)
# preChemo.smooth.fitted <- fitted.DHR(preChemo.H, 
#                                    preChemo.smoothed.process$smoothed.states) 

### Chemo: Kalman Filter
Chemo.filtered.process <- Kalman.Filter(y = Chemo, F = Chemo.F, 
                                           H = Chemo.H, 
                                           G = Chemo.G,
                                           Q = Chemo.Q, 
                                           sigma2 = Chemo.sigma2/10)
Chemo.fitted <- fitted.DHR(Chemo.H,
                              Chemo.filtered.process$states)
### Chemo: Kalman Smoothing
# Chemo.smoothed.process <- Kalman.Smoother(Chemo.filtered.process)
# Chemo.smooth.fitted <- fitted.DHR(Chemo.H, 
#                                     Chemo.smoothed.process$smoothed.states) 



#### Visualisation ###
plot(1:400, model, col = "grey", pch = 20)
points(1:400, model, col = "grey", type = "o", pch = 20)
lines(1:400, c(preChemo.fitted, preChemo.fitted[200], Chemo.fitted[-1]), 
      type = "l", col = "purple", 
      pch = 20, lty = 3, lwd = 2, cex = 0.9)
lines(1:400, c(preChemo.fitted, preChemo.fitted[200], Chemo.fitted[-1]), type = "o", col = "purple", 
      pch = 20, lty = 3, lwd = 1, cex = 0.9)



###### Trying to  make prediction
prediction.x <- preChemo.filtered.process$prediction.x
prediction.P <- preChemo.filtered.process$prediction.P
F <- preChemo.F
H <- preChemo.H

k <- 100
a <- F %*% prediction.x[[100]] 
prediction <- H[[101]] %*% a
prediction
preChemo[]

prediction.x <- filtered.process$prediction.x
prediction.P <- filtered.process$prediction.P


preChemo.F




############################   2   ################################






##### Getting DNFL model 
model <- get.DNFL(N = 400, omega = 200, mu = 2, 
                  alpha = 10, SDE = TRUE)
mRNA <- model$M; mRNA <- mRNA[201:400]
N <- length(mRNA)
mRNA <- c(((1:100)/20) * mRNA[1:100], ((100:1)/20) * mRNA[101:200])

mRNA <- mRNA + rnorm(N, 0, 20)
plot(1:N, mRNA, type = "o", col = "red", pch = 20)

##### Periodogram
I <- periodogram(mRNA)$I
freq <- periodogram(mRNA)$freq
n.tilda <- length(freq)

##### Selecting p
performance <- performance.ar.p(mRNA, p = 30, 
                                plot = TRUE, legend = TRUE)
p <- which(performance$AIC == min(performance$AIC))
p 

##### Fitting AR(p)
fit <- arima(mRNA, order = c(p, 0, 0), method = "ML")
coef <- fit$coef[1:p]
sigma2 <- fit$sigma2


##### Deciding how many and which frequencies to use

# Auxiliary AR(p) object
temp.ar.fit <- arma.spec(ar = coef, var.noise = sigma2, n.freq = n.tilda)
# Auxiliary  Frequencies
temp.ar.freq <- temp.ar.fit$freq; temp.ar.freq <- temp.ar.freq[-1]
# Auxiliary  Spectrum
temp.ar.spec <- temp.ar.fit$spec; temp.ar.spec <- temp.ar.spec[-1]
plot(temp.ar.freq, log(temp.ar.spec), type = "l")

# Peaks 
peaks <- findpeaks(as.vector(temp.ar.spec)); peaks
driving.freq <- c()
for(j in 1:nrow(peaks)) {
  driving.freq[j] <- temp.ar.freq[peaks[j, 2]]
}


## AR(p) Object 
ar.fit <- arma.spec(ar = coef, var.noise = sigma2, n.freq = 500)
## Frequencies (no omega_0)
ar.freq <- ar.fit$freq; ar.freq <- ar.freq[-1]
## Spectrum (no omega_0)
ar.spec <- ar.fit$spec; ar.spec <- ar.spec[-1]
plot(ar.freq, log(ar.spec), type = "l")


### Setting foundmental frequency and other harmonics
omega.mRNA <- c(0, driving.freq); omega.mRNA
### Modelling the evolution of the coefficient process 
structures.mRNA <- c("IRW", rep("RW", length(driving.freq)))


###### Estimation of the NVR coefficients

# Linear Least Square approximation
NVR.LSS <- NVR.Linear.LS(ar.spec, sigma2, ar.freq, omega.mRNA, 
                         structures.mRNA)
# Relative approximation of the variances
sigma.LSS <- c(sigma2, sigma2 * NVR.LSS)
# Getting full DHR spectrum
middle.DHRspectrum <- fullDHR.spectrum(ar.freq, sigma.LSS, omega.mRNA)
plot(ar.freq, log(middle.DHRspectrum), type = "l", col = "green")
lines(ar.freq, log(ar.spec), type = "l", col = "red") 


# Non linear approximation
NVR <- optim(NVR.LSS, E, 
             freq = ar.freq, empirical.spectrum = ar.spec, 
             omega.j = omega.mRNA, 
             structures = structures.mRNA, 
             sigma2 = sigma2, 
             method = "Nelder-Mead")$par
# Relative approximation of the variances
sigma.est <- c(sigma2, sigma2 * NVR)
# Getting full DHR spectrum
full.DHRspectrum <- fullDHR.spectrum(ar.freq, sigma.est, omega.mRNA)

plot(ar.freq, log(full.DHRspectrum), col = "green", type = "l")
lines(ar.freq, log(ar.spec), type = "l", col = "red")

# Non linear approximation
NVR <- optim(NVR, E, 
             freq = ar.freq, empirical.spectrum = ar.spec, 
             omega.j = omega.mRNA, 
             structures = structures.mRNA, 
             sigma2 = sigma2, 
             method = "Nelder-Mead")$par
# Relative approximation of the variances
sigma.est <- c(sigma2, sigma2 * NVR)
# Getting full DHR spectrum
full.DHRspectrum <- fullDHR.spectrum(ar.freq, sigma.est, omega.mRNA)

plot(ar.freq, log(full.DHRspectrum), col = "green", type = "l")
lines(ar.freq, log(ar.spec), type = "l", col = "red")



############# Filtering ###########

# Getting F matrix
F.IRW <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
F.RW <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
F <- bdiag(c(list(F.IRW), replicate(length(omega.mRNA) - 1, 
                                    F.RW, simplify = FALSE)))
# Getting G matrix
G.IRW <- matrix(c(0, 0, 0, 1), nrow = 2, ncol = 2)
G.RW <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
G <- bdiag(c(list(G.IRW), replicate(length(omega.mRNA) - 1,
                                    G.RW, simplify = FALSE)))
# Getting Q matrix
Q <- diag(rep(sigma2 * NVR, 1, each = 2))
# Getting H list
H <- get.H(omega.mRNA, N)


### Kalman Filter
filtered.process <- Kalman.Filter(y = mRNA, F = F, H = H, 
                                  G = G, Q = Q, sigma2 = sigma2/5)
### Kalman Smoothing
smoothed.process <- Kalman.Smoother(filtered.process)



######   Visualisation   #####
fitted <- fitted.DHR(H, filtered.process$states)
smooth.fitted <- fitted.DHR(H, smoothed.process$smoothed.states) 

plot(1:N, mRNA, col = "grey", pch = 20)
points(1:N, mRNA, type = "o", col = "grey", pch = 20)
lines(1:N, fitted, type = "o", col = "purple", pch = 20)
lines(1:N, smooth.fitted, type = "o", col = "purple", pch = 20)

get.periods(smooth.fitted)

residuals <- fitted - mRNA
residuals.norm <- residuals/(sd(residuals))


plot.ts(residuals.norm, pch = 20, col = "red")
qqnorm(residuals.norm)


## Trying Prediction
# Getting objects
states <- filtered.process$states
covariance.states <- filtered.process$covariance.states
prediction.x <- filtered.process$prediction.x
prediction.P <- filtered.process$prediction.P
F <- filtered.process$F
H <- filtered.process$H
G <- filtered.process$G
sigma2 <- filtered.process$sigma2
Q <- filtered.process$Q; Q.r <- Q/sigma2

# a_t(k) & P_t
t <- 125
K <- 25
a <- list(); a[[1]] <- prediction.x[[t]]
# R <- list(); R[[1]] <- covariance.states[[t]]
f <- c();
Q <- c();

for(k in 1:K) {
  # k-step ahead forecast for states
  a[[k + 1]] <- as.vector(F %*% a[[k]])
}

plot(125:149, f, col = "green", type = "o")

length(f)
as.vector(H[[60 + 1]] %*% a[[1]])

mRNA[60]
mRNA[61]

plot(61:70, f, type = "o")

length(f)

optim(ciao, my.E, 
      freq = ar.freq, empirical.spectrum = ar.spec, 
      omega.j = omega.model, 
      structures = structures.model, 
      sigma2 = sigma2, 
      method = "Nelder-Mead")$par

my.E(log(NVR.LSS), ar.freq, ar.spec, 
     omega.model, structures.model, sigma2)

my.E(c(10, 10, 1), ar.freq, ar.spec, 
     omega.model, structures.model, sigma2)


log(exp(NVR.LSS))

exp(log(NVR.LSS))

ciao <- log(NVR.LSS)
my.E()
