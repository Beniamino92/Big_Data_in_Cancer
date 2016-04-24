data(faithful)
summary(faithful)

Y <- faithful$eruptions
x <- faithful$waiting

X <- cbind(rep(1, 272), x)

plot(y ~ x, pch = 19)
model1 <- lm(y ~ x)
coef <- model1$coefficients

E <- function(beta, Y, X) {
  sum((Y - X%*%beta)^2)
}



Beta <- c(0, 0)
E(Beta, Y, X)
optim(c(0, 0), E, Y = Y, X = X)
my.coef <- optim(c(0, 0), E, y = y, x = x)$par
my.coef
coef

#### Function to evaluate the logarithmic error.

### LA FINE E' SBAGLIATA 

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

E(NVR.LSS.estimate)


NVR.LSS.estimate


for(r in 1:(R + 1)) {
  X[, r] <- S(frequency, omega.j[r], type = structures[r])
}



source("Steps_DHR_Estimation_Algorithm.R")

test <- E(NVR = NVR.LSS.estimate,
          freq = freq , 
          empirical.spectrum = spec, 
          omega.j = omega.j.AP, 
          structures = structures.AP,
          sigma2 = sigma2)

# Getting data
data.AP <- AirPassengers

T <- length(data.AP)

# Harmonics
omega.j.AP <- c(0, 1/12, 1/6, 1/4, 1/3, 1/(2.4))

R <- length(omega.j.AP) - 1

# IRW Trend, RW harmonic components
structures.AP <- c("IRW", rep("RW", 5))
length(structures.AP)
# Getting spectrum and residuals variance from AR(P) fitted to AirPassengers
AP.AR14 <- estimate.AR(data = data.AP, p = 14, n.freq = length(data.AP))
freq <- AP.AR14$frequency
spec <- AP.AR14$spectrum
sigma2 <- AP.AR14$sigma2


R
length(omega.j.AP)



# Let's try to create a function E(sigma) (no logaritmic)
# and optimisit, instead that with least squares, with optim.
R <- length(omega.j) - 1

X <- matrix(NA, nrow = length(spec), ncol = (R + 1))
structures <- c("IRW", rep("RW", 5))
omega.j <- c(0, 1/12, 1/6, 1/4, 1/3, 1/(2.4))
for(r in 1:(R + 1)) {
  X[, r] <- S(freq, omega.j[r], type = structures[r])
}
spec

X %*% NVR.LSS.estimate
sigma2

E <- function(NVR) {
  fitted <- sigma2*(X%*%NVR + (1/(2*pi)))
  error <- sum((spec - fitted)^2)
  return(error)
}

E.log <- function(NVR) {
  fitted <- sigma2*(X%*%NVR + (1/(2*pi)))
  error <- sum((log(spec) - log(fitted))^2)
  return(error)
}

NVR.estimate <- optim(NVR.LSS.estimate, E.log)$par

NVR.estimate


as.matrix(optim(rep(0, 6), E)$par, nrow = 6)
NVR.LSS.estimate

fitted
E(NVR.LSS.estimate)
fitted <- sigma2*(X%*%NVR.LSS.estimate + (1/(2*pi)))
fitted
spec
NVR.LSS.estimate


# Sigma estimate
sigma.est <- c(sigma2, sigma2 * NVR.estimate)

# Plotting AR(14) spectrum and DHR spectrum with estimated sigma (LLS)
plot(freq, log(spec), type = "l", col = "red", lwd = 2)
lines(freq, log(fullDHR.spectrum(freq,  sigma.est, omega.j.AP)),
      lwd = 2, col = "green")


# Optimising logarithmic error: (THERE ARE WARNINGS)
optimising.NVR <- optim(NVR.LSS.estimate, E, 
                        freq = freq, empirical.spectrum = spec,
                        omega.j = omega.j.AP, structures = structures.AP,
                        sigma2 = sigma2, method = "Nelder-Mead")

ciao <- spg(NVR.LSS.estimate, E, 
    freq = freq, empirical.spectrum = spec,
    omega.j = omega.j.AP, structures = structures.AP,
    sigma2 = sigma2)$par
ciao
NVR.estimate
optimising.NVR$par
NVR.nonlinear.estimate
