# Backward Smoothing 

# I need this objects:
states
covariance.states
F
H
G
Q.r <- Q / sigma2

#### Filtering:
filtered.process <- Kalman.Filter(y = data.AP, F = F, H = H, G = G, 
                                  Q = Q, sigma2 = sigma2)
class(filtered.process)
states <- filtered.process$states
covariance.states <- filtered.process$covariance.states
prediction.x <- filtered.process$prediction.x
prediction.P <- filtered.process$prediction.P
sigma2 <- filtered.process$sigma2

length(states)
length(covariance.states)
length(prediction.x)
length(prediction.P)

t <- T
L <- list()
smoothed.states <- list()
smoothed.covariance.states <- list()





# Initial Step
smoothed.states[[T]] <- states[[T]]
L[[T]] <- as.matrix(rep(0, n), ncol = 1)
smoothed.covariance.states[[T]] <- covariance.states[[T]]

t <- (T-1)

# Lagrange Multiplier
L[[t]] <- t((diag(n) - covariance.states[[t + 1]] %*% t(H[[t + 1]]) %*% 
             H[[t + 1]])) %*% ( t(F) %*% L[[t + 1]] - t(H[[t + 1]]) %*%
                                  (y[t + 1] - H[[t + 1]] %*% states[[t + 1]]))


# Smoothed states
temp <- (smoothed.states[[t + 1]] + G %*% Q.r %*% t(G) %*% L[[t]])
smoothed.states[[t]] <- solve(F, smoothed.states[[t+1]] +  G %*% Q.r %*% t(G)
                              %*% L[[t]])

smoothed.states[[t]] <- F %*% smoothed.states[[t + 1]] - G %*% Q.r %*% t(G) %*% 
  L[[t]]


# Smoothed covariances
smoothed.covariance.states[[t]] <- covariance.states[[t]] + 
  covariance.states[[t]] %*% t(F) %*% solve(prediction.P[[t + 1]]) %*%
  (smoothed.covariance.states[[t + 1]] - prediction.P[[t + 1]]) %*%
  solve(prediction.P[[t + 1]]) %*% F %*% covariance.states[[t]]


A <- solve(prediction.P[[t + 1]], 
      smoothed.covariance.states[[t + 1]] - prediction.P[[t + 1]])

smoothed.covariance.states[[t]] <- covariance.states[[t]] + 
  covariance.states[[t]] %*% t(F) %*% A %*% 
  solve(prediction.P[[t + 1]], F %*% covariance.states[[t]])

# %*% solve(prediction.P[[t + 1]], F %*% covariance.states[[t]])
      
solve(A, F %*% covariance.states[[t]])


# Warning: I can't evaluate F^(-1), because it's singular.
# Therefore I used the generalised inverse.

require(matrixcalc)
require(MASS)



# DOESN'T WORK AT ALL!!! 

# NEED TO REDO IT.
# Backward Smoothing Function

Backward.Smoothing <- function(filtered.process) {

  
  if(!(class(filtered.process) == "Kalman.Filter")) {
    stop("object \'filtered.process\' has to be of class \'Kalman.Filter\'")
  }
  
  # Required elements for backward smoothing, given by Kalman.Filter
  states <- filtered.process$states
  covariance.states <- filtered.process$covariance.states
  prediction.P <- filtered.process$prediction.P
  y <- filtered.process$y
  F <- as.matrix(filtered.process$F)
  H <- filtered.process$H
  G <- filtered.process$G
  Q.r <- (filtered.process$Q)/(filtered.process$sigma2)
  
  T <- length(y)
  n <- dim(Q.r)[1]
  
  # New elements for backward smoothing
  L <- list()
  smoothed.states <- list()
  smoothed.covariance.states <- list()
  
  # Initial Step
  smoothed.states[[T]] <- states[[T]]
  L[[T]] <- as.matrix(rep(0, n), ncol = 1)
  smoothed.covariance.states[[T]] <- covariance.states[[T]]
  
  # F could be singular, in that case we use the generalised inverse
  if(is.singular.matrix(F)) {
    F.inv <- ginv(F, tol = sqrt(.Machine$double.eps))
  }
  else {
    F.inv <- solve(F)
  }
  
  ### Backward Pass Smoothing Equations:
  for(t in (T-1):1) {
    
    # Lagrange Multiplier
    L[[t]] <- t((diag(n) - covariance.states[[t + 1]] %*% 
                   t(H[[t + 1]]) %*%  H[[t + 1]])) %*% 
      ( t(F) %*% L[[t + 1]] - t(H[[t + 1]]) %*% 
          (y[t + 1] - H[[t + 1]] %*% states[[t + 1]]))
    
    # Smoothed States
    smoothed.states[[t]] <- F.inv %*% (smoothed.states[[t + 1]] + 
                                            G %*% Q.r %*% t(G) %*% L[[t]])
  
    # Smoothed Covariance States
    A <- solve(prediction.P[[t + 1]], 
               smoothed.covariance.states[[t + 1]] - prediction.P[[t + 1]])
    
    smoothed.covariance.states[[t]] <- covariance.states[[t]] + 
      covariance.states[[t]] %*% t(F) %*% A %*% 
      solve(prediction.P[[t + 1]], F %*% covariance.states[[t]])
  }
  
  return(list(smoothed.states = smoothed.states,
              smoothed.covariance.states = smoothed.covariance.states))
}

tryal <- Backward.Smoothing(filtered.process)

states.smooth <- tryal$smoothed.states


str(filtered.process, 1)

class(filtered.process) <- "Kalman.Filter"
class(filtered.process)
class(filtered.process) == "Kalman.Filter"


# Trying stuff with solve

n <- 5
A <- matrix(rnorm(n*n), n, n)
B <- matrix(rnorm(n*n), n, n)
C <- matrix(rnorm(n*n), n, n)
D <- matrix(rnorm(n*n), n, n)

solution <- D + D %*% t(C) %*% solve(A) %*% (B - A) %*% solve(A) %*% C %*% D

D + D %*% t(C) %*% solve(A, (B - A)) %*% solve(A) %*% C %*% D

E <- solve(A, (B - A))

D + D %*% t(C) %*% E %*% solve(A, C %*% D)

solve(A)

for(t in (T-1):1) {
  print(t)
}

covariance.states[[t + 1]]



# Function to implement Forward Pass Filtering Equations (Kalman)

test.F <- as.matrix(F)
class(test.F)
require(matrixcalc)
is.singular.matrix(test.F)
det(test.F)
require(MASS)
ginv(test.F, tol = sqrt(.Machine$double.eps)) %*% test.F

solve(A) 
A <- matrix(rnorm(5 * 5), 5, 5)
solve(A)
ginv(A, tol = sqrt(.Machine$double.eps))


F.mod <- bdiag(replicate(R , F.IRW, simplify = FALSE))
F.mod <- as.matrix(F.mod)

solve(F)
solve(F.mod)
ginv(F.mod)
            