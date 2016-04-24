####### Function: Drift for mRNA molecules 
drift.M <- function(M, P, parameters.M) {
  
  # Getting parameters
  v1 <- parameters.M[1]
  v2 <- parameters.M[2]
  k1 <- parameters.M[3] 
  k2 <- parameters.M[4] 
  hc <- parameters.M[5]
  
  # Drift
  out <- ((v1*k1)/((k1 + P)^hc)) - ((v2*M)/(k2 + M))
  
  return(out)
}


####### Function:  Drift for protein
drift.P <- function(M, P, parameters.P) {
  
  # Getting parameters
  alpha <- parameters.P[1]
  delta.P <- parameters.P[2]
  
  # WARNING: We're using just M(t), not g(M(t))
  
  # Drift
  out <- (alpha * M) - (delta.P * P)
  
  return(out)
}


##### Function: Solving the ODE, with Euler scheme.

Euler <- function(M.start, P.start, parameters.M,
                  parameters.P, h, T, start = 0) {
  
  # NUmber of steps required
  n.steps <- (T - start) / h
  
  # Variables and initial conditions
  M <- c(); M[1] <- M.start
  P <- c(); P[1] <- P.start
  t <- start
  
  # Euler method
  for(n in 1:n.steps) {
    M[n + 1] <- M[n] + h * drift.M(M[n], P[n], parameters.M)
    P[n + 1] <- P[n] + h * drift.P(M[n], P[n], parameters.P)
  }
  
  return(list(M = M, P = P))
}


# Testing Euler Scheme.

parameters.M <- c(1.5, 1.3, (0.2)^4, 0.2, 4)
parameters.P <- c(2, 0.5)

test <- Euler(M.start = 2, P.start = 0, 
              parameters.M, parameters.P, 
              h = 1e-1, T = 100)

plot.ts(test$M) 
plot.ts(test$P[10:100])











# Trying sovling ODE with Euler, simple scenario:
# URL: http://www.sosmath.com/diffeq/system/euler/example1/answer.html


# Drift for M
drift.M <- function(t, M, P) {
  out <- -2*t*M + 3*(P^2)
  return(out)
}

# Drift fot P
drift.P <- function(t, M, P) {
  out <- -3*(M^2)*(1 - P)
}

# Initial condition:
M0 <- -1
P0 <- 2

h <- 1/1000

M <- c()
P <- c()

M[1] <- M0
P[1] <- P0

t <- 0
N <- 100



for(n in 2:N) {
  t <- t + h
  M[n] <- M[n - 1] + h * drift.M(t, M[n - 1], P[n - 1])
  P[n] <- P[n - 1] + h * drift.P(t, M[n - 1], P[n - 1])
}

plot.ts(P)


Euler <- function(M.start, P.start, T, h, start = 0) {
  
  # NUmber of steps required
  n.steps <- (T - start) / h
  
  # Variables and initial conditions
  M <- c(); M[1] <- M.start
  P <- c(); P[1] <- P.start
  t <- start
  
  # Euler method
  for(n in 1:n.steps) {
    t <- t + h
    M[n + 1] <- M[n] + h * drift.M(t, M[n], P[n])
    P[n + 1] <- P[n] + h * drift.P(t, M[n], P[n])
  }
  
  return(list(M = M, P = P))
}

test <- Euler(M.start = -1, P.start = 2, T = 1, 
              h = 1e-1)

length(seq(0, 2, by = 1e-2))
length(test$M)
plot(seq(0, 1, by = 1e-1), test$M, type = "l")
test$M



############ Retry again:

T <- 10
start <- 0
h <- 1/10
tau <- 2

n.steps <- (T - start)/h; n.steps

M.start <- 100
P.start <- 0

M <- c(); M[1] <- M.start
P <- c(); P[1] <- P.start

t <- 0

for(n in 1:n.steps) {
  
  t <- t + h
  M[n + 1] <- M[n] + h * drift.M(M[n], P[n], parameters.M)
  
  if(t < tau) {
    P[n + 1] <- P[n] + h * drift.P(M[n], P[n], parameters.P)
  }
  
  else{
    cat('t: ',t, '\n')
    cat('n: ',n, '\n')
    cat('M[', (t - tau) + 1,'] = ', M[(t - tau) + 1])
    cat('\n')
    
    g.M <- M[(t - tau) + 1]
    P[n + 1] <- P[n] + h * drift.P(g.M, P[n], parameters.P)
  }
}

M[-2]

plot.ts(P[-1])

# NUmber of steps required
n.steps <- (T - start) / h

# Variables and initial conditions
M <- c(); M[1] <- M.start
P <- c(); P[1] <- P.start
t <- start








############ Questions and Considerations:

# Initial value for ODE/SDE
# Not sure about g(M(t)).
