#SIR on RK4
S0 <- 145
I0 <- 5
R0 <- 0

N <- S0 + I0 + R0

t0 <- 0
tFinal <- 60
tStep <- 1

t <- seq(t0, tFinal, by = tStep)

gamma <- 1/7
beta <- c(0.05, 0.10, 0.125, 0.14, 0.15, 0.20, 0.30)

colors <- c("blue", "green", "orange", "red", "purple", "brown", "black")
plot(NULL,
     xlim = c(t0, tFinal),
     ylim = c(0, N),
     xlab = "Time",
     ylab = "Population infected",
     main = "SIR Model")

SIR_func <- function(sir_vec, beta, gamma, N){
  S <- sir_vec[1]
  I <- sir_vec[2]
  R <- sir_vec[3]
  
  dS <- -beta * I * S / N
  dI <- beta * I * S / N - gamma*I
  dR <- gamma * I
  
  return(c(dS, dI, dR))
}

for(b in beta){
    currvec <- c(S0, I0, R0)  
  
    recordS <- numeric(length(t))
    recordI <- numeric(length(t))
    recordR <- numeric(length(t))
    
    recordS[1] <- S0
    recordI[1] <- I0
    recordR[1] <- R0
    
    for (i in 1:(length(t) - 1)) {
      
      k1 <- SIR_func(currvec, b, gamma, N)
      
      k2 <- SIR_func(currvec + tStep*k1/2,
                b, gamma, N)
      
      k3 <- SIR_func(currvec + tStep*k2/2,
                b, gamma, N)
      
      k4 <- SIR_func(currvec + tStep*k3,
                b, gamma, N)
      
      currvec <- currvec + tStep/6 * (k1 + 2*k2 + 2*k3 + k4)
      
      recordS[i + 1] <- currvec[1]
      recordI[i + 1] <- currvec[2]
      recordR[i + 1] <- currvec[3]
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

