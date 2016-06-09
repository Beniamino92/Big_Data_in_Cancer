library(pracma)


####### Function: Euler's method to simulate DNLF model #######

# parameters.M: hc, v1, k1, v2, k2, mu, sigma2
# parameters.B: v3, v4

Euler <- function(h, T, parameters.M, parameters.P, SDE = FALSE) {
  
  # Setting Parameters (M)
  hc <- parameters.M[1]
  v1 <- parameters.M[2]
  k1 <- parameters.M[3]
  v2 <- parameters.M[4]
  k2 <- parameters.M[5]
  mu <- parameters.M[6]
  alpha <- parameters.M[7]

  
  # Setting Parameters (P)
  v3 <- parameters.P[1]
  v4 <- parameters.P[2]
  
  # Step size
  n.steps <- T/h; n.steps
  
  # Initial conditions:
  M <- c(); M[1] <- 10
  P <- c(); P[1] <- 10
  
  #### Euler's sytem
  for(n in 1:n.steps) {
    
    ############ mRNA ##########
    
    drift.M <- (v1/(1 + (P[n]/k1)^hc)) - v2*(M[n]/(k2 + M[n]))
    
    if(SDE == FALSE) {
      M[n + 1] <- M[n] + h * drift.M
    }
    else {
      diffusion.M <- sqrt( (v1/(1 + (P[n]/k1)^hc)) + v2*(M[n]/(k2 + M[n])))
      dW.M <- rnorm(1, 0, sd = sqrt(h))
      
      M[n + 1] <- M[n] + h * drift.M + dW.M * diffusion.M
      # Negative values can appear due to discretisaion
      if(M[n + 1] < 0) {
        M[n + 1] <- 0
      }
    }
    
    
    # Delay feedback
    past <- seq(from = 1, to = n)
    M.past <- M[n:1]
    death.prob <- dgamma(past, alpha, (alpha/mu))
    delay.M <- M.past %*% death.prob
    
    
    ############ Protein ##########
    
    drift.P <- v3*delay.M - v4*P[n]
    if(SDE == FALSE) {
      P[n + 1] <- P[n] + h * drift.P
    }
    else {
      diffusion.P <- sqrt(v3*delay.M + v4*P[n])
      dW.P <- rnorm(1, 0, sd = sqrt(h))
      P[n + 1] <- P[n] + h * drift.P  + dW.P * diffusion.P
    }
  }
  
  return(list(M = M, P = P))
}




#### Function: get distribution of the period:
get.periods <- function(data) {
  
  # Positions of peaks
  peaks <- findpeaks(data)
  # Number of peaks
  n.peaks <- nrow(peaks)
  
  periods <- c()
  # Getting periods
  for(i in 1:(n.peaks - 1)) {
    periods[i] <- peaks[(i + 1), 2] - peaks[i, 2]
  }
  return(periods)
}


### Function: MODEL 1
get.DNFL <- function(N, omega, mu, alpha,  SDE = FALSE) {
  
  # Setting Parameters (M)
  hc <- 4
  v1 <- 0.5*omega
  k1 <- 2*omega
  v2 <- 0.3*omega
  k2 <- 0.2*omega
  parameters.M.1 <- c(hc, v1, k1, v2, k2, mu, alpha)
  # Setting Parameters (P)
  v3 <- 2
  v4 <- 0.5
  parameters.P.1 <- c(v3, v4)
  
  # Euler-Maruyama scheme
  if(SDE == TRUE) {
    x <- Euler(h = 1, N, parameters.M.1, parameters.P.1, SDE = TRUE)
  }
  else{
    x <- Euler(h = 1, N, parameters.M.1, parameters.P.1, SDE = FALSE)
  }
  # Not considering M.start, and P.start
  x[[1]] <- x[[1]][-1]
  x[[2]] <- x[[2]][-1]
  
  return(x)
}


# Function: plot mRNA and protein
double.plot.DNFL <- function(model) {
  
  mRNA <- model$M; protein <- model$P
  n <- length(mRNA)
  dat <- data.frame(mRNA = mRNA, protein = protein) 
  
  double.plot <- ggplot() +   
    geom_point(data = dat, aes(x = c(1:n), y = mRNA), colour = "red", size = 2) + 
    geom_line(data = dat, aes(x = c(1:n), y = mRNA, linetype = "o"), colour = "red") + 
    geom_point(data = dat, aes(x = c(1:n), y = protein), colour = "blue", size = 2) + 
    geom_line(data = dat, aes(x = c(1:n), y = protein, linetype = "o"), colour = "blue") + 
    xlab("Time") + ylab("") + theme(legend.position = "none") + 
    ylim(0, 1500)
  
  return(double.plot)
}


# Euler2: Function for series with different periods, oscilaltions, bla bla

Euler2 <- function(h, T, parameters.M, parameters.P, 
                   mu2, alpha2, SDE = FALSE) {
  
  # Setting Parameters (M)
  hc <- parameters.M[1]
  v1 <- parameters.M[2]
  k1 <- parameters.M[3]
  v2 <- parameters.M[4]
  k2 <- parameters.M[5]
  mu <- parameters.M[6]
  alpha <- parameters.M[7]
  
  # Setting Parameters (P)
  v3 <- parameters.P[1]
  v4 <- parameters.P[2]
  
  # Step size
  n.steps <- T/h; n.steps
  
  # Initial conditions:
  M <- c(); M[1] <- 10
  P <- c(); P[1] <- 10
  
  #### Euler's sytem
  for(n in 1:n.steps) {
    
    ############ mRNA ##########
    
    drift.M <- (v1/(1 + (P[n]/k1)^hc)) - v2*(M[n]/(k2 + M[n]))
    
    if(SDE == FALSE) {
      M[n + 1] <- M[n] + h * drift.M
    }
    else {
      diffusion.M <- sqrt( (v1/(1 + (P[n]/k1)^hc)) + v2*(M[n]/(k2 + M[n])))
      dW.M <- rnorm(1, 0, sd = sqrt(h))
      
      M[n + 1] <- M[n] + h * drift.M + dW.M * diffusion.M
      # Negative values can appear due to discretisaion
      if(M[n + 1] < 0) {
        M[n + 1] <- 0
      }
    }
    
    
    # Delay feedback
    if(n  >= (n.steps/2)) {
      mu <- mu2
      alpha <- alpha2
    }
    
    past <- seq(from = 1, to = n)
    M.past <- M[n:1]
    death.prob <- dgamma(past, alpha, (alpha/mu))
    delay.M <- M.past %*% death.prob
    
    
    ############ Protein ##########
    
    drift.P <- v3*delay.M - v4*P[n]
    if(SDE == FALSE) {
      P[n + 1] <- P[n] + h * drift.P
    }
    else {
      diffusion.P <- sqrt(v3*delay.M + v4*P[n])
      dW.P <- rnorm(1, 0, sd = sqrt(h))
      P[n + 1] <- P[n] + h * drift.P  + dW.P * diffusion.P
    }
  }
  
  return(list(M = M, P = P))
}


### Function: get the model using Euler2
get.DNFL2 <- function(N, omega, mu, alpha, mu2, alpha2, SDE = FALSE) {
  
  # Setting Parameters (M)
  hc <- 4
  v1 <- 0.5*omega
  k1 <- 2*omega
  v2 <- 0.3*omega
  k2 <- 0.2*omega
  parameters.M.1 <- c(hc, v1, k1, v2, k2, mu, alpha)
  # Setting Parameters (P)
  v3 <- 2
  v4 <- 0.5
  parameters.P.1 <- c(v3, v4)
  
  # Euler-Maruyama scheme
  if(SDE == TRUE) {
    x <- Euler2(h = 1, N, parameters.M.1, parameters.P.1,
                mu2, alpha2,  SDE = TRUE)
  }
  else{
    x <- Euler2(h = 1, N, parameters.M.1, parameters.P.1,
                mu2, alpha2, SDE = FALSE)
  }
  
  # Not considering M.start, and P.start
  x[[1]] <- x[[1]][-1]
  x[[2]] <- x[[2]][-1]
  
  return(x)
}





