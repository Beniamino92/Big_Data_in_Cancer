# URL http://homepage.univie.ac.at/robert.kunst/USPOP.TSM

uspop <- USPOP
names(uspop)[1] <- "population"

T <- 21
time <- c(1:T)

uspop <- cbind(uspop, time)

plot(time, uspop$population, type = "o", col = "red")

# Fitting Quadratic Regression
quadratic.reg <- lm(uspop$population ~ uspop$time + I((uspop$time)^2))
summary(quadratic.reg)

# Adding fitted model
plot(time, uspop$population, type = "o", col = "red")
lines(time, fitted(quadratic.reg), col = "blue", lty = 2)

# Prediction
quadratic.prediction <- function(coeff, time) {
  return(coeff[1] + coeff[2] * time + coeff[3]*(time^2))
}

# Evaluating prediction
quadratic.prediction(quadratic.reg$coefficients, time = 22)

# Plotting residuals
plot(time, quadratic.reg$residuals, type = "o", col = "red")

