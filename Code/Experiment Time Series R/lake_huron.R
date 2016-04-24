lake <- LakeHuron
time.lake <- seq(1:length(lake))

plot(time.lake, lake, col = "red", type = "o")

# Linear model
lake.linear.model <- lm(lake ~ time.lake)
summary(lake.linear.model)

# Plotting fitted model
lines(time.lake, fitted(lake.linear.model), col = "blue", lty = 2)

# Prevision
fcast <- forecast(lake.linear.model, newdata=data.frame(time.lake = 99))
plot(fcast)

# Plotting residuals
plot(time.lake, lake.linear.model$residuals, type = "o")