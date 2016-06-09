setwd("C:/Users/Beniamino/Desktop/Mini_Project_1/Code/Inference_on_periodicity_of_Circadian_Time_Series")
source("DNFL_model.R")
setwd("C:/Users/Beniamino/Desktop/Mini_Project_1/Code/On_Boostrapping_Kernel_Spectral_Estimates")
source("Final_Kernel_Smoothing_Spectrum.R")

grid <- 0:2500/100


######### Model 1 ######## (true period around  23)
model1 <- get.DNFL(N = 600, omega = 200, mu = 6, 
                   alpha = 10, SDE = TRUE)
mRNA1 <- model1$M; mRNA1 <- mRNA1[201:600]

# ODE
model1.ODE <- get.DNFL(N = 600, omega = 200, mu = 6, 
                       alpha = 10, SDE = FALSE)
mRNA1.ODE <- model1.ODE$M
mRNA1.ODE <- mRNA1.ODE[201:600]
model1.periods <- get.periods(mRNA1.ODE); mean(model1.periods)

# Double plots
model1[[1]] <- model1[[1]][201:400]
model1[[2]] <- model1[[2]][201:400]

double.plot.DNFL(model1)
plot(grid, dgamma(grid, 10, (10/6)), type = "l", 
     xlab = "", ylab = "", lwd = 2)

# Spectral Resampling
I.1 <- periodogram(mRNA1)$I
freq.1 <- periodogram(mRNA1)$freq
freq.rad.1 <- periodogram(mRNA1)$freq.rad
c.val.1 <- get.c(freq.rad.1, I.1)

SR.model1 <- Spectral.Resampling(freq.rad.1, I.1, c.val.1, R = 1000)
avg.SR.model1 <- rowMeans(SR.model1)
plot(freq.1, avg.SR.model1, type = "l", col = "red")
1/freq.1[which(avg.SR.model1 == max(avg.SR.model1))] # 23.52
Spectral.Resampling.CI(freq.1, SR.model1)





######### Model 2 ######## (true period around 11.69)
model2 <- get.DNFL(N = 600, omega = 200,  mu = 2,
                   alpha = 10, SDE = TRUE)
mRNA2 <- model2$M; mRNA2 <- mRNA2[201:600]

# ODE
model2.ODE <- get.DNFL(N = 400, omega = 200, mu = 2, 
                       alpha = 10, SDE = FALSE)
mRNA2.ODE <- model2.ODE$M
mRNA2.ODE <- mRNA2.ODE[201:400]
model2.periods <- get.periods(mRNA2.ODE); mean(model2.periods)

# Double plots
model2[[1]] <- model2[[1]][201:400]
model2[[2]] <- model2[[2]][201:400]
double.plot.DNFL(model2)
plot(grid, dgamma(grid, 10, (10/2)), type = "l", 
     xlab = "", ylab = "", lwd = 2)

# Spectral Resampling
I.2 <- periodogram(mRNA2)$I
freq.2 <- periodogram(mRNA2)$freq
freq.rad.2 <- periodogram(mRNA2)$freq.rad
c.val.2 <- get.c(freq.rad.2, I.2)

SR.model2 <- Spectral.Resampling(freq.rad.2, I.2, c.val.2, R = 1000)
avg.SR.model2 <- rowMeans(SR.model2)
plot(freq.2, avg.SR.model2, type = "l", col = "red")
1/freq.2[which(avg.SR.model2 == max(avg.SR.model2))] # 11.76
Spectral.Resampling.CI(freq.2, SR.model2)



##### Model3 = Model 1 + Model 2 #####

model3 <- c()
model3[[1]] <- model1$M + model2$M
model3[[2]] <- model1$P + model2$P
names(model3) <- c("M", "P")

mRNA3 <- model3$M

# Double plots
model3[[1]] <- model3[[1]][201:400]
model3[[2]] <- model3[[2]][201:400]
double.plot.DNFL(model3)

# Spectral Resampling
I.3 <- periodogram(mRNA3)$I
freq.3 <- periodogram(mRNA3)$freq
freq.rad.3 <- periodogram(mRNA3)$freq.rad
c.val.3 <- get.c(freq.rad.3, I.3)

SR.model3 <- Spectral.Resampling(freq.rad.3, I.3, c.val.3, R = 1000)
avg.SR.model3 <- rowMeans(SR.model3)
plot(freq.3, avg.SR.model3, type = "l", col = "red")
1/freq.3[which(avg.SR.model3 == max(avg.SR.model3))] # 11.76
Spectral.Resampling.CI(freq.3, SR.model3)




# New functions for plotting CI. plus periodogram 
Spectral.Resampling.CI.2 <- function(freq, SR.sample, periodogram) {
  
  # Getting dimensions and setting frequencies
  n.tilda <- length(freq)
  n <- n.tilda*2
  
  
  # Creating confidence intervals for spectrum estimate 
  conf.int <- matrix(NA, nrow = n.tilda, ncol = 2)
  for(i in 1:n.tilda) {
    conf.int[i, ] <- quantile(SR.sample[i, ], probs = c(.05, .95))
  }
  # Average
  mean.spectrum.SR <- rowMeans(SR.sample)
  
  # Auxiliary data.frame for plotting
  dat <- data.frame(freq = freq, spec = mean.spectrum.SR, period = periodogram[1:n.tilda],
                    ci.up = conf.int[, 2], ci.low = conf.int[, 1])
  # Plotting mean and C.I
  p <- ggplot() +
    geom_line(data = dat, 
              aes(x = freq, y = spec, linetype = "c"),
              size = 1.0, colour = "red") +
    geom_line(data = dat, 
              aes(x = freq, y = period), 
              size = 0.7, colour = "blue") +
    geom_line(data = dat, 
              aes(x = freq, y = ci.up), 
              size = 0.6, colour = "black", linetype = "dashed") + 
    geom_line(data = dat, 
              aes(x = freq, y = ci.low), 
              size = 0.6, colour = "black", linetype = "dashed") +
    xlab("Frequency") + ylab("Power Spectrum") +
    theme(legend.position = "none", 
          plot.title = element_text(size = 19)) +
    scale_x_continuous(minor_breaks = seq(min(freq), max(freq), 0.05))
  
  
  return(list(spec = mean.spectrum.SR, 
              ci.up = conf.int[, 2], ci.low = conf.int[, 1], 
              plot = p))
}

Spectral.Resampling.CI.2(freq.1, SR.model1, I.1)
Spectral.Resampling.CI(freq.1, SR.model1)
1/freq.1[which(avg.SR.model1 == max(avg.SR.model1))]; freq.1[which(avg.SR.model1 == max(avg.SR.model1))]
1/freq.2[which(avg.SR.model2 == max(avg.SR.model2))]; freq.2[which(avg.SR.model2 == max(avg.SR.model2))]

