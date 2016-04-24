# Function to obtain AIC and BIC for the fitting to an
# AR model of order p, for i in 1:p.

performance.ar.p <- function(residuals.data, p, plot = TRUE, legend = TRUE) {
  
  T <- length(residuals.data)
  
  AIC <- rep(0, p) 
  BIC <- rep(0, p)
  
  for(k in 1:p) {
    fit <- ar(residuals.data, order = k, aic = FALSE)
    sigma2 <- var(fit$resid, na.rm = TRUE)
    BIC[k] <- log(sigma2) + (k*log(T)/T)
    AIC[k] <- log(sigma2) + ((T+2*k)/T)
  }
  
  if (plot == TRUE) {
    IC <- cbind(AIC, BIC)
    ts.plot(IC, type = "o", xlab = "p",
            ylab = "AIC and BIC", col = c("red", "blue"),
            lwd = 2)
    points(which(BIC == min(BIC)), BIC[which(BIC == min(BIC))], 
           col = "black", pch = 19, lwd = 4)
    points(which(AIC == min(AIC)), AIC[which(AIC == min(AIC))], 
           col = "black", pch = 19, lwd = 4)
    if (legend == TRUE) {
    legend("bottomright", c("AIC", "BIC"), lty = 1, lwd = 2,
           col = c("red", "blue"))
    }
  }
  
  return(list(AIC = AIC, BIC = BIC))
}


# Normalise
normalise <- function(x) {
  return( (x - min(x))/(max(x) - min(x)))
}


# Pseudo-Spectrum Generalised Random Walk (GRW)
pseudo.spectrum.GRW <- function(omega, sigma2.eta, sigma2,
                                type = "IRW", alpha = 1) {
  if(type == "IRW") {
    alpha <- 1
  }
  else if(type == "RW") {
    alpha <- 0
  }
  
  spectrum <- (1/(2*pi)) * (  (sigma2.eta/((1 + (alpha^2) - (2*alpha*cos(omega)) ) *
                                             (2 - 2*cos(omega)))) + sigma2   )
  return(spectrum)
}



# Function to compute pseudo-spectrum single frequency DHR term

singleDHR.spectrum <- function(omega, omega.j, sigma2, type) {
  
  if(type == "IRW") {
    out <- (1/(2*pi)) * ((sigma2/(4*(1 - cos(omega - omega.j))^2)) +
                           (sigma2/(4*(1 - cos(omega + omega.j))^2)))
    return(out)
  }
  
  else if(type == "RW") {
    out <- (1/(2*pi)) * ((sigma2/(2*(1 - cos(omega - omega.j)))) +
                           (sigma2/(2*(1 - cos(omega + omega.j)))))
    return(out)
  }
}



# Function to compute: S(omega, omega.j) (p.9 DHR)
S <- function(omega, omega.j, type) {
  
  if(type == "IRW") {
    out <- (1/(2*pi)) * ((1/(4*(1 - cos(omega - omega.j))^2)) +
                           (1/(4*(1 - cos(omega + omega.j))^2)))
    return(out)
  }
  
  else if (type == "RW") {
    out <- (1/(2*pi)) * ((1/(2*(1 - cos(omega - omega.j)))) +
                           (1/(2*(1 - cos(omega + omega.j)))))
    return(out)
  }
}



###### Function to compute pseudo-spectrum full DHR (specific example)

# sigma: vector of variances
#        1st component is the observational variance
#        2nd component is the trend variance
#        3rd...Rth components are harmonic variances.

# omega.j: vector of fondamental frequency and harmonics
#          omega.j[1] must be 0 (if the trend is included)

fullDHR.spectrum <- function(omega, sigma, omega.j) {
  
  R <- length(omega.j) - 1
  
  # Trend component:
  spectrum <- sigma[2] * S(omega, omega.j[1], type = "IRW")
  print(spectrum)
  
  # Harmonic components:
  for(j in 1:R) {
    spectrum = spectrum + sigma[j + 2] * S(omega, omega.j[j + 1], type = "RW")
  }
  
  return(spectrum + (sigma[1]/(2*pi)))
  
}











########## Consideration (Problems):

# Ask Mark about normalisation of the spectrum, in this case.

# What is the pseudo spectrum?

# Show my result to Barbel about spectrum Figure.3

# Regarding, the estimation of the hyper parameters via minimisation of an
# error function. Do we optimise it with numerical methods? such as optim?

# Almost got the graph!! Figure 3, I just don't get the picks
# of the spectrums. :( )

# Speak with Mark, about the curious case of the frequencies....,
# and evaluation of p.(8 - 9)

# If I don't have enough frequency, for example let's say if I have 100 freq, it
# doesn't work. Some values are Inf!

# Is sigma2 the variance of the innovations? 

# Different number of frequencies change the estimates of sigma!

# Smoothing by elimination of high frequency component?

# Not sure about (p.4); Why for RW you use the differences? 

# Talk about the structures to put on Kalman Filter. 

# How to structure the time variable matrix, H_t:
# Is that [cos(omega0_t) , sin(omega0 _t), ... , 
#           cos(omegaR_t) , sin(omegaR_t)]


# Question: considering the trend component, as omega_0 = 0; this leads
# to the sine part of alpha_0 cos(w_0 *  t) + beta_0 sin(w_0 * t), alwats to be 
# zero right? So what's the point of including beta_0 ? 


# In general: not sure how many parameters there are in the model? 12? 24?
# I'm a bit lost under this point of view.


# What to put, as first step of the Kalman, Filter. I mean the initial 
# values of the states? at t = 1.

# Do I need to multiply the coefficient I got, to something? Like sigma2.

# What's the difference between cyclical and seasonal?

# Initial condition for backward smoothing? X_{N|N} = \hat{X_N}?

# What about the order of the equations for the smoothing?! L(t)

# Try different G, and think about it.

# Not sure about the section "Estimation of the SSM disturbances"

# Problems inverting F.

# What's the variance of the innovations?

# How can we do confidence intervals?