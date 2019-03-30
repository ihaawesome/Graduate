library(numDeriv)
library(dplyr)
library(ggplot2)

# ---------------------------------------------------------------------------------

  mybisec <- function(f, a, b, E=10^-8) {
    
    maxiter = 1000
    error = 1
    niter = 0
    xt = (a+b)/2
    
    while (niter <= maxiter & error >= E) {

      if (grad(f,a)*grad(f,xt) < 0) { b <- xt }
      else if (grad(f,a)*grad(f,xt) > 0) { a <- xt }
      
      xt_1 <- xt
      xt <- (a+b)/2
      
      error <- abs(xt-xt_1)
      niter <- niter + 1
    }
    
    return(data.frame(x = xt, fx = f(xt), niter))
  }


  mynewton <- function(f, x0, E=10^-8) {

    maxiter = 1000
    error = 1
    niter = 0
    xt = x0
    
    while (niter <= maxiter & error >= E)  {
      
      xt_1 <- xt
      xt <- xt - grad(f, xt)/as.double(hessian(f, xt))
      error <- abs(xt-xt_1)
      niter <- niter + 1
    }
    
    result = data.frame(x = xt, fx = f(xt), niter)
    return(result)
  }

  
  mysecant <- function(f, x0, x1, E=10^-8) {
    
    maxiter = 1000
    error = 1
    niter = 0
    xt = x1
    xt_1 = x0
    
    while (niter <= maxiter & error >= E) {
      
      xt_new <- xt - grad(f, xt)*(xt-xt_1)/(grad(f, xt)-grad(f, xt_1))
      xt_1 <- xt
      xt <- xt_new
      
      error <- abs(xt-xt_1)
      niter <- niter + 1
    }
    
    return(data.frame(x = xt, fx = f(xt), niter))
  }


  bisectime <- function(f, a, b) { system.time( for (i in 1:1000) mybisec(f, a, b) ) / 1000 }   
  newtontime <- function(f, x0) { system.time( for (i in 1:1000) mynewton(f, x0) ) / 1000 }
  secanttime <- function(f, x0, x1) { system.time( for (i in 1:1000) mysecant(f, x0, x1) ) / 1000 }  

# ---------------------------------------------------------------------------------

# function 1
  myf1 <- function(x) { 
    f = log(x)/(1+x)
    return(f)
  }

  mybisec(myf1, 1, 5) ; bisectime(myf1, 1, 5)
  mynewton(myf1, 1) ; newtontime(myf1, 1)
  mysecant(myf1, 1, 5) ; secantime(myf1, 1, 5)   
  optimize(myf1, c(1, 5), maximum = T)
  
  x_star = mynewton(myf1, 1)$x
  
  x = seq(0.1, 10, by = 0.1) ; fx = myf1(x)
  plot(x, fx, main = paste("f(x) = log(x)/(1+x)"), type = "l", lwd = 2)
  points(x_star, myf1(x_star), pch = 16, col = "red")
  abline(v = x_star, h = myf1(x_star), lty = c(2,2), col = c("red","black"))

# function 2
  myf2 <- function(x) {
    f = exp(x)*cos(x^2+1)
    return(f)
  }
 
  mybisec(myf2, 3, 5) ; bisectime(myf2, 3, 5)
  mynewton(myf2, 3) ; newtontime(myf2, 3)
  mysecant(myf2, 3, 5) ; secanttime(myf2, 3, 5)
  optimize(myf2, c(3, 5), maximum = T)
  
  x2_star = mynewton(myf2, 5)$x
  x2_newton1 = mynewton(myf2, 3)$x ; x2_bisec = mybisec(myf2, 3, 5) ; x2_secant = mysecant(myf2, 3, 5)
  
  x2 = seq(1, 5, length = 100) ; fx2 = myf2(x2)  
  plot(x2, fx2, main = paste("f(x) = exp(x)*cos(x^2+1)"), type = "l", lwd = 2)
  points(x2_star, myf2(x2_star), pch = 16, col = "red")
  points(x2_newton1, myf2(x2_newton1), pch = 16, col = "blue")
  points(x2_bisec, myf2(x2_bisec), pch = 16, col = "blue")
  points(x2_secant, myf2(x2_secant), pch = 16, col = "blue")
  abline(v = c(x2_star, x2_newton1, x2_bisec, x2_secant), 
         h = c(myf2(x2_star), myf2(x2_newton1), myf2(x2_bisec) , myf2(x2_secant)), 
         lty = c(2,2,2,2), col = c("red", "black", "black", "black", "black", "black", "black", "black"))
  
  myf3 <- function(x) {
    f = 1/(1+x^2)+sqrt(x)
    return(f)
  }
  
  x3 = seq(0, 3, length = 100) ; fx3 = myf3(x3)  
  plot(x3, fx3, main = paste("f(x) = 1/(1+x^2)+sqrt(x)"), type = "l", lwd = 2)
  
  mybisec(myf3, 0.1, 1) ; bisectime(myf3, 0.1, 1)
  mynewton(myf3, 0.1) ; newtontime(myf3, 0.1)
  mysecant(myf3, 0.1, 1) ; secanttime(myf3, 0.1, 1)
  mysecant(myf3, 0.1, 0.7) ; secanttime(myf3, 0.1, 0.7)
  optimize(myf3, c(0, 1), maximum = T)
  
  x3_star = mybisec(myf3, 0.1, 1)$x
  x3_local = mysecant(myf3, 0.1, 1)$x
  plot(x3, fx3, main = paste("f(x) = 1/(1+x^2)+sqrt(x)"), type = "l", lwd = 2)
  points(x3_star, myf3(x3_star), pch = 16, col = "red")
  points(x3_local, myf3(x3_local), pch = 16, col = "blue")
  abline(v = c(x3_star, x3_local), h = c(myf3(x3_star), myf3(x3_local)), 
         lty = c(2,2,2,2), col = c("red", "red", "black", "black"))
  
  
  
# exercise 2.1 -------------------------------------------------------------------
 
  sample = c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 
             3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)
  
# a
  
  llcauchy <- function(theta) {
    density = dcauchy(sample, location = theta, scale = 1)
    loglik = sum(log(density))
    return(loglik)
  }
  
  tb <- data.frame(theta = seq(-50, 100, length = 300)) %>% mutate(ltheta = 0)
  for (i in 1:300) { tb$ltheta[i] <- llcauchy(tb$theta[i]) }
  
  plot(tb$theta, tb$ltheta, type = "l", lwd = 2, 
       xlab = "theta", ylab = "l(theta)", main = "log-likelihood of Cauchy(theta, 1)")
  
  start = c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
  newton21 = data.frame(theta = numeric(length(start)), loglik = numeric(length(start)), 
                        niter = numeric(length(start)))
  for (j in 1:3) {
    for (i in 1:length(start)) {
      newton21[i,j] <- ifelse(is(try(mynewton(llcauchy, x0 = start[i]), silent = T),"try-error"), 
                              t(rep(NA,3)), mynewton(llcauchy, x0 = start[i])[j]) 
    }
  }
  try(newtontime(llcauchy, x0 = start[1]), silent = T)
  newtontime(llcauchy, x0 = start[2])
  newtontime(llcauchy, x0 = start[3])
  newtontime(llcauchy, x0 = start[4])
  newtontime(llcauchy, x0 = start[5])
  newtontime(llcauchy, x0 = start[6])
  try(newtontime(llcauchy, x0 = start[7]), silent = T)
  newtontime(llcauchy, x0 = start[8])
  newtontime(llcauchy, x0 = start[9])
 
  mean(sample) 
  mynewton(llcauchy, x0 = mean(sample))
  newtontime(llcauchy, x0 = mean(sample)) 
  
  plot(tb$theta, tb$ltheta, type = "l", lwd = 2, 
       xlab = "theta", ylab = "l(theta)", main = "log-likelihood of Cauchy(theta, 1)")
  abline(v = newton21$theta, col = "red", lty = 2)
  abline(v = mynewton(llcauchy, x0 = mean(sample))$x, col = "blue", lty = 2)


# b
  mybisec(llcauchy, -1, 1) ; bisectime(llcauchy, -1, 1)
  
# d  
  mysecant(llcauchy, -2, 1) ; secanttime(llcauchy, -2, 1)
  mysecant(llcauchy, -3, 3) ; secanttime(llcauchy, -3, 3)
  mysecant(llcauchy, -1, 1) ; secanttime(llcauchy, -1, 1)

# e
  set.seed(100)
  sample2 = rnorm(20, 7, 1)
  
  llnorm <- function(theta) {

    loglik = c()
    for (i in 1:length(theta)) {
      density = dnorm(sample2, mean = theta[i], sd = 1)
      loglik[i] = sum(log(density))
    }
    return(loglik)
  }
  
  tb2 <- data.frame(theta = seq(-10, 20, length = 100)) %>% mutate(ltheta = llnorm(theta))
  plot(tb2$theta, tb2$ltheta, type = "l", lwd = 2, 
       xlab = "theta", ylab = "l(theta)", main = "log-likelihood of N(theta, 1)")

  optimize(llnorm, c(5,10), maximum = T) # 7.107867
  
  mybisec(llnorm, 5, 10) ; bisectime(llnorm, 5, 10)
  mynewton(llnorm, 5) ; newtontime(llnorm, 5)
  mysecant(llnorm, 5, 10) ; secanttime(llnorm, 5, 10)
  
  mean(sample2)
  mynewton(llnorm, mean(sample2)) ; newtontime(llnorm, mean(sample2))
  
  plot(tb2$theta, tb2$ltheta, type = "l", lwd = 2, 
       xlab = "theta", ylab = "l(theta)", main = "log-likelihood of N(theta, 1)")
  abline(v = mean(sample2), col = "red", lty = 2)
  
# -------------------------------------------------------------------------------------------------
  
# exercise 2.2

  data = c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96, 
           2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52, 2.50)
      
  ex2.2 <- function(x, theta) {
    fx = ifelse((x>=0 & x<=2*pi), (1-cos(x-theta))/(2*pi), 0)
  }

  llex2.2 <- function(theta) {
      density = ex2.2(data, theta)
      loglik = sum(log(density))
    return(loglik)
  }

# a  
  theta = seq(-pi, pi, length = 200)
  ltheta = c()
  for (i in 1:200) { ltheta[i] = llex2.2(theta[i]) }  

  plot(theta, ltheta, type = "l", lwd = 2, main = "log-likelihood of f(x)")  

# b
  mme = asin(mean(data)-pi)  
# c
  mynewton(llex2.2, mme) ; newtontime(llex2.2, mme) 
  mynewton(llex2.2, -2.7) ; newtontime(llex2.2, -2.7) 
  mynewton(llex2.2, 2.7) ; newtontime(llex2.2, 2.7)

# d
  newton22_x = c() ; newton22_fx = c()
  for (i in 1:length(theta)) {
    newton22_x[i] <- ifelse(is(try(mynewton(llex2.2, x0 = theta[i]), silent = T),"try-error"), 
                            NA, mynewton(llex2.2, x0 = theta[i])$x)
    newton22_fx[i] <- ifelse(is(try(mynewton(llex2.2, x0 = theta[i]), silent = T),"try-error"), 
                             NA, mynewton(llex2.2, x0 = theta[i])$fx)    
  }
  newton22 = data.frame(theta, ltheta, theta_star = newton22_x, ltheta_star = newton22_fx)
  newton22 <- filter(newton22, theta_star >= -pi & theta_star <= pi)  
  
  plot(theta, ltheta, type = "l", lwd = 2, 
       xlab = "theta", ylab = "l(theta)", main = "log-likelihood of f(x)")  
  points(newton22$theta_star, newton22$ltheta_star, col = "red", pch = 16)
  
  ggplot(newton22) + geom_line(aes(theta, ltheta), size = 0.8) + 
    geom_point(aes(theta, ltheta_star, color = factor(round(theta_star,4))), size = 1.5) + theme_bw() + 
    labs(color = "theta_star")
  
# e
  mynewton(llex2.2, x0 = newton22$theta[176]) ; newtontime(llex2.2, x0 = newton22$theta[176])
  mynewton(llex2.2, x0 = newton22$theta[177]) ; newtontime(llex2.2, x0 = newton22$theta[177])
  