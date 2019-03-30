library(dplyr)
library(ggplot2)
rm(list = ls())

############################## Functions ##############################
# Calculate IS mean, sd (Unstandardized)
estimate <- function(x, w) {
  w <- w / sum(w)
  ismean <- sum(w*x) / sum(w)
  issd <- sqrt( sum(w*(x-ismean)^2) / (1-sum(w^2)) )
  return(c(mcmean = mean(x), mcsd = sd(x), ismean = ismean, issd = issd))
}

# Generate Real Data
paths <- function(n = 100) {
# hyperparameters  
  tau = 0.000001 ; eta = 0.005

# initial values
  x = c(0.01, numeric(n-1)) ; y = c(20, numeric(n-1))
  vx = c(0.002, numeric(n-2)) ; vy = c(-0.06, numeric(n-2))
  
  for (t in 2:n) {
    vx[t] <- rnorm(1, vx[t-1], sqrt(tau)) ; vy[t] <- rnorm(1, vy[t-1], sqrt(tau))
    x[t] <- x[t-1] + vx[t-1] ; y[t] <- y[t-1] + vy[t-1]
  }
  table <- data.frame(t = 1:n, x, y, vx, vy) %>% mutate(z = atan(y/x) + rnorm(n, 0, eta))
  return(table)
}

realxy <- paths(n = 50)
ggplot(realxy) + theme_light() + theme(legend.position = "") + 
  scale_color_continuous() + labs(x = "x", y = "y") +
  geom_path(aes(x, y), color = "lightgray", size = 2) + 
  geom_point(aes(x, y), size = 3, shape = 1) + 
  geom_point(aes(x[1], y[1]), color = "lightgray", shape = 20, size = 10) +
  geom_text(aes(label = "ORIGIN", x[1], y[1]))

z <- realxy$z


############################## Tracking ##############################   

weight <- function(zt, xt, yt) { dnorm(zt, atan(yt/xt), 0.005) }

init <- c(0.000001, 0, 19.9, 0.001, -0.05)
TRACK1 <- function(t, init, niter=1000, RS = T, MOVE = T) {

# initial values  
  tau0 <- init[1] ; x0 <- init[2] ; y0 <- init[3] ; vx0 <- init[4] ; vy0 <- init[5] 
  
  # objects
  vx <- matrix(nrow = t, ncol = niter) ; vy <- vx ; x <- vx ; y <- vx ; weights <- vx
  vx[1,] <- vx0 ; vy[1,] <- vy0 ; x[1,] <- x0 ; y[1,] <- y0 ; 
  
  tau <- matrix(tau0, t, niter) 
  weights[1,] <- 1 / niter 
  
  for (ti in 2:t) {
    
    for (n in 1:niter) {
      vx[ti, n] <- rnorm(1, vx[ti-1, n], sqrt(mean(tau[ti,]))) ; x[ti, n] <- x[ti-1, n] + vx[ti, n]
      vy[ti, n] <- rnorm(1, vy[ti-1, n], sqrt(mean(tau[ti,]))) ; y[ti, n] <- y[ti-1, n] + vy[ti, n]
    
      weights[ti, n] <- weights[ti-1, n] * weight(z[ti], x[ti, n], y[ti, n])
      
      if (n%%niter==0) { 
        print(paste("t =", ti, "(", round(x[ti, n], 3), ", ", round(y[ti, n], 3), ")"), sep = "")
      }
    }  
    weights[ti,] <- weights[ti,] / sum(weights[ti,])  
    
    # 1) resample step
    if (RS) {  
      resample <- sample(niter, size = niter, prob = weights[ti,], replace = T)
      
      vx[1:ti,] <- vx[1:ti, resample] ; vy[1:ti,] <- vy[1:ti, resample]
      x[1:ti,] <- x[1:ti, resample] ; y[1:ti,] <- y[1:ti, resample]
      
      weights[ti,] <- 1 / niter # reset weights
    } 
    # 2) move
    if (MOVE) {
      tau[ti,] <- rgamma(niter, 2, 0.5*(sum(diff(vx[ti,])^2) + sum(diff(vy[ti,])^2))) # ~ Exp(0.5(x^2+y^2))
    }

  }
  output <- list(x = x, y = y, tau = rowMeans(tau), weight = weights[t,])
  return(output)
}
  
true <- c(0.000001, 0.01, 20, 0.002, -0.06) 
init <- c(0.0001, 0.015, 20.005, 0.0, -0.05)

track1 <- TRACK1(t = 50, init = init, niter = 1000, RS = T, MOVE = T)
output <- data.frame(x = rowMeans(track1$x), 
                     y = rowMeans(track1$y))

track2 <- TRACK1(t = 50, init = init, niter = 1000, RS = T, MOVE = F)
output2 <- data.frame(x = rowMeans(track2$x), 
                      y = rowMeans(track2$y))

track3 <- TRACK1(t = 50, init = init, niter = 1000, RS = F, MOVE = F)
output3 <- data.frame(x = apply(track3$x, 1, function(x) weighted.mean(x, track3$weight)), 
                      y = apply(track3$y, 1, function(x) weighted.mean(x, track3$weight)))

track4 <- TRACK1(t = 50, init = init, niter = 1000, RS = F, MOVE = T)
output4 <- data.frame(x = rowMeans(track4$x), 
                      y = rowMeans(track4$y))

ggplot(realxy) + theme_light() + theme(legend.position = "top") + 
  labs(x = "x", y = "y") +
  scale_color_manual("", values = c(RESAMPLEMOVE = "blue", RESAMPLE = "orange", 
                                    MOVE = "green", NOT = "red")) +
  geom_point(aes(x, y), color = "gray", size = 3) + 
  geom_text(aes(label = "ORIGIN", x[1], y[1]), color = "black") +
  geom_point(aes(output$x, output$y), color = "blue", shape = 8) +
  geom_path(aes(output$x, output$y, color = "RESAMPLEMOVE"), size = 1) + 
  geom_path(aes(output2$x, output2$y, color = "RESAMPLE"), linetype = "dashed", size = 1) +
  # geom_path(aes(output3$x, output3$y, color = "NOT"), linetype = "dashed", size = 1) +
  geom_path(aes(output4$x, output4$y, color = "MOVE"), linetype = "dashed", size = 1)

ggplot(realxy) + theme_light() + theme(legend.position = "") + 
    labs(x = "x", y = "y") +
    geom_point(aes(x, y), color = "gray", size = 3) + 
    geom_text(aes(label = "ORIGIN", x[1], y[1]), color = "black") +
    geom_point(aes(output$x, output$y), color = "blue", shape = 8) +
    geom_path(aes(output$x, output$y, color = desc(t)), size = 1) 
  
plot(track1$tau[-1], type = "l", xlab = "t", ylab = "tau", main = "Posterior mean of Tau")
acf(track1$tau, main = "Series of Tau")

print(isestimate1 <- estimate(track1$y[50,], rep(1/1000, 1000)))
print(isestimate2 <- estimate(track2$y[50,], rep(1/1000, 1000)))
print(isestimate3 <- estimate(track3$y[50,], rep(1/1000, 1000)))

ismean2 <- apply(track2$x, 1, function(x) estimate(x, track2$weights)[1])
issd2 <- apply(track2$x, 1, function(x) estimate(x, track2$weights)[2])

ggplot() + geom_bar(aes(ismean1))
plot(ismean1)
