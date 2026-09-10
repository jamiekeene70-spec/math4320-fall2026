# Euler's method SIR model

S0 <- 145
I0 <- 5
R0 <- 0

t0 <- 0
tFinal <- 60
tStep <- 1

t <- seq(t0, tFinal, by = tStep)

N <- S0 + I0 + R0


gamma <- 1/7

beta <- c(0.05, 0.10, 0.125, 0.14, 0.15, 0.20, 0.30)

# Create empty plot
colors <- c("blue", "green", "orange", "red", "purple", "brown", "black")
plot(NULL,
     xlim = c(t0, tFinal),
     ylim = c(0, N),
     xlab = "Time",
     ylab = "Population infected",
     main = "SIR Model")
     

for (b in beta) {
  
  Scurr <- S0
  Icurr <- I0
  Rcurr <- R0
  
  recordS <- numeric(length(t))
  recordI <- numeric(length(t))
  recordR <- numeric(length(t))
  
  recordS[1] <- S0
  recordI[1] <- I0
  recordR[1] <- R0
  
  for (i in 1:(length(t) - 1)) {
    
    # SIR equations
    dS <- -b * Scurr * Icurr / N
    dI <- b * Scurr * Icurr / N - gamma * Icurr
    dR <- gamma * Icurr
    
    # Euler update
    recordS[i + 1] <- Scurr + tStep * dS
    recordI[i + 1] <- Icurr + tStep * dI
    recordR[i + 1] <- Rcurr + tStep * dR
    
    Scurr <- recordS[i + 1]
    Icurr <- recordI[i + 1]
    Rcurr <- recordR[i + 1]
  }
  
  # Add each curve to the graph
  lines(t, recordI, col = colors[which(beta == b)], lwd = 2)
}

# Legend
legend("topright",
       legend = paste("β =", beta),
       col = colors,
       lty = 1,
       lwd = 2)

