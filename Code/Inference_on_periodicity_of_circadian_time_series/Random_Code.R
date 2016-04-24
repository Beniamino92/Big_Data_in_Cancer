# When we work with frequencies we can either use (T even)
# freq <- 0:(T/2 - 1)/T
# freq <- 2*pi*(0:(T/2 - 1))/T


# Periodogram function:
periodogram <- function(data) {
  T <- length(data)
  I <- (abs(fft(data))^2)/T
  return(I)
}

# Gaussian Kernel:
gaussian.kernel <- function(x, b = 1, scaled = FALSE) {
  if(!(FALSE)) {
    out <- ((1/(b*sqrt(2*pi))) * exp(-((x/b)^2)/2))
    return(out)
  }
  else {
    out <- (1/(sqrt(2*pi))) * exp(-(x^2)/2)
    return(out)
  }
}


# Smoothed.Periodogram
smoothed.periodogram <- function(omega.k, omega, b, I) {
  
}

# Data
x <- as.vector(AirPassengers)
n <- length(AirPassengers)

# Frequency:
n.tilda <- n/2
k <- 0:((n.tilda) - 1)
freq2 <- 2*pi*k/n

# Getting periodogram
I <- periodogram(x, n)
plot(freq, log(I[1:n.tilda]), type = "l")

log(I[1:n.tilda])
length(freq)
freq

freq <- 1:(T/2 - 1)/T

c <- 10
# Bandwith
b <- c * 1/(T^5); b

I <- periodogram(x, (T/2))
plot(freq, log(I[1:(T/2)]), type = "l")


test <- c(freq, 0, freq)
smoothed.periodogram()
length(test)

freq <- seq(from = 0, to = 0.5 - (1/T), by = 1/T)



# Proviamo frequency paper
k.test <- seq(from = - n.tilda, to = ( n - 1), by = 1)

omega.k <- 0.7

freq <- 2*pi*k.test/n

gaussian.kernel(freq - omega.k, b = 1, scaled = TRUE)

# Proviamo frequency solito Mark

# This seems right to me.
freq <- (-T/2):(T - 1)/T
omega.k <- 0.2
ker <- gaussian.kernel(freq - omega.k,
                       b = 0.3701072, scaled = TRUE)
I <- periodogram(freq)
sum(ker * freq)/sum(ker)


# NO 

# Function to evaluate smoothed periodogram (1.2)
smoothed.periodogram <- function(freq, b, I) {
  
  # Length of the series
  T <- length(I)
  
  # Extended frequency and periodogram
  auxiliary.freq <- (-T/2):(T-1)/T
  auxiliary.I <- c(I[(T/2):1], I[2:(T - 1)])
  
  # Scaled kernel evaluation
  ker <- gaussian.kernel(auxiliary.freq - freq,
                         b = b, scaled = TRUE)
  
  smoothed.period <- sum(ker * auxiliary.freq)/sum(ker)
  return(smoothed.period)
}

freq <- 0:(T/2 - 1)/T
freq

omega <- 0:(T/2 - 1)/T

negro <- smoothed.periodogram(omega, b = (T^(-(1/5))), 
                     I = I)

ciao <- c()
for(t in 1:length(omega)) {
  ciao[t] <- smoothed.periodogram(omega[t], b = (T^(-(1/5))), 
                                  I = I)
}

plot(omega, ciao, type = "l")

I <- periodogram(as.vector(AirPassengers))
smoothed.periodogram(freq, b = (T^(-(1/5))), 
                     I = I)
periodogram(0)
T^(-(1/5))

c(1, 2, 3) * c(4, 5, 6)

rave <- function(money, drugs, music) {
  return("smashed")
}

freq2 <- seq(0, to = pi, len = (T/2))

length(freq2)

plot(freq, log(I[1:(T/2)]), type = "l")
plot(freq2, log(I[1:(T/2)]), type = "l")

log(I[1:(T/2)])
(freq[13])^(-1)

((freq2[13])/(2*pi))^(-1)

length(freq)
log(I[1:(T/2)])


# FINAL CONFIRMATION EVER TO FIX IN MY FUCKING MIND

# We can use either:
freq <- 0:(T/2 - 1)/T
freq2 <- 2*pi*(0:(T/2 - 1))/T
  
  
  
  n.tilda <- n/2
k <- 0:((n.tilda) - 1)
freq2 <- 2*pi*k/n
  
ker = kernel("modified.daniell", c(19,19))

sum(ker$coef)

str(ker)
plot(k, main = "Weights")
plot(freq,ar1.spec(0.5,sigma2=1,num.freq=T/2),col="red",lwd = 3,type="l",
     main = "2M+1 = 77, Daniell Weights",ylab = "Power Spectrum", xlab = "Frequency")
spec.pgram(xt,kernel=k,add=TRUE,lwd=3)




#### LET'S TRY AGAIN

# Data
x <- as.vector(AirPassengers)
n <- length(AirPassengers)

# Frequency:
n.tilda <- n/2
k <- 1:n.tilda
freq <- 2*pi*k/n


# Periodogram function:
periodogram <- function(data) {
  n <- length(data)
  I <- (abs(fft(data))^2)/(n*2*pi)
  return(I)
}


I <- periodogram(x)



# Function to evaluate the scaled kernel
scaled.kernel <- function(x, b) {
  out <- (1/b) * (1/sqrt(2*pi)) * exp(-((x/b)^2)/2)
  return(out)
}

# Auxiliary stuff

k.aux <- seq(from = -n.tilda, to = n - 1); length(k.aux)
I.aux <- c(I[n.tilda:1], 0,  I[1:(n - 1)]); length(I.aux)
freq.aux <- 2*pi*k.aux/n

x <- 0.2
temp <- freq.aux - x
temp2 <- scaled.kernel(temp, b = 0.2)

temp2 %*% I.aux / sum(temp2)
temp2

scaled.kernel(freq.aux - x, b = 0.2)

length(temp)

# Smoothed kernel periodogram:
smoothed.periodogram <- function(x, I, b) {
  
  n <- length(I)
  n.tilda <- n/2
  k.aux <- seq(from = -n.tilda, to = n - 1); length(k.aux)
  I.aux <- c(I[n.tilda:1], 0,  I[1:(n - 1)]); length(I.aux)
  freq.aux <- 2*pi*k.aux/n
  
  temp <- scaled.kernel(freq.aux - x, b )
  
  numerator <- temp %*% I.aux
  denominator <- sum(temp)
  
  return(numerator/denominator)
}

get.smoothed.periodogram <- function(freq, I, b) {
  out <- c()
  for(i in 1:length(freq)) {
    out[i] <- smoothed.periodogram(freq[i], I, b)
  }
  return(out)
}

smoothed.periodogram(0.2, I, b = 0.2)
get.smoothed.periodogram

ciao <- c()
for(i in 1:length(freq)) {
  ciao[i] <- smoothed.periodogram(freq[i], I, b = 0.08)
}


length(ciao)

plot(freq, log(get.smoothed.periodogram(freq, I, b = 0.10)), 
     type = "l")



