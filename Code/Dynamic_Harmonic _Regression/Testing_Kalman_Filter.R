######  Testing Kalman Filter #####

source("Spectrums.R")
source("Steps_DHR_Estimation_Algorithm.R")
source("Testing_Steps.R")
source("Kalman_Filter.R")


# Elements defined in "Testing_Steps.R"
data.AP 
T
sigma2
NVR.nonlinear.estimate


####### Preparing system, and observational matrixes #####
R <- 1 + 5 # Trend + Harmonic TVP's
n <- 2*R
omega.j <- c(0, 1/12, 1/6, 1/4, 1/3, 1/(2.4))

# Matrix H
H <- get.H(omega.j, T)

# Matrix F
F.IRW <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
F.RW <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
F <- bdiag(c(list(F.IRW), replicate((R - 1), F.RW, simplify = FALSE)))

# Matrix G
G.IRW <- matrix(c(0, 0, 0, 1), nrow = 2, ncol = 2)
G.RW <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
G <- bdiag(c(list(G.IRW), replicate((R - 1), G.RW, simplify = FALSE)))

# Matrix Q
Q <- diag(rep(sigma2 * NVR.nonlinear.estimate, 1, each = 2))



#### Filtering:
filtered.process <- Kalman.Filter(y = data.AP, F = F, H = H, G = G, 
                                  Q = Q, sigma2 = sigma2)

states <- filtered.process$states
covariance.states <- filtered.process$covariance.states

# Fitted values
fitted <- fitted.DHR(H, states)                               

# Plotting true and fitted value (after Filtering)
plot(1:t, data.AP, type = "o")
lines(1:t, fitted, type = "o", col = "red", pch = 20)


# Smoothing
smoothed.filtered.process <- Backward.Smoothing(filtered.process)

smoothed.states <- smoothed.filtered.process$smoothed.states
smoothed.fitted <- fitted(H, smoothed.states)
plot(1:t, data.AP, type = "o")
lines(1:t, fitted, type = "o", col =" red", pch = 20)
lines(1:t, smoothed.fitted, type = "o", col = "green", pch = 20)
