library(dplyr)
library(ggplot2)
options("scipen" = 100)
  
# Rieman  
  myrieman <- function(f, a, b, start=1, E=10^-8)
  { 
    niter = 0 ; maxiter = 100 ; error = 1 
    
    # initial value  
    R = (b-a)*f(a)
    k <- start
    
    while (error >= E & niter <= maxiter)
    { R_0 <- R
    
    n <- 2^k ; h <- (b-a)/n
    i <- 0:(n-1)
    R <- h * sum( f(a+i*h) )
    
    error <- abs(R - R_0) / abs(R_0 + 10^-3)  
    niter <- niter + 1
    k <- k + 1
    
    print(paste("error = ", error, " niter = ", niter, sep = ""))
    }
    
    return(R)
  }
  
  
# Trapezoidal
  mytrapezoid <- function(f, a, b, start=1, E=10^-8)
  { 
    niter = 0 ; maxiter = 100 ; error = 1
    
    # initial value  
    t = (b-a)*f(a)
    k <- start 
    
    while (error >= E & niter <= maxiter)
    { t_0 <- t
    
    n <- 2^k ; h <- (b-a)/n
    i <- 1:(n-1)
    t <- h/2*f(a) + h*sum( f(a+i*h) ) + h/2*f(b)
    
    error <- abs(t - t_0) / abs(t_0 + 10^-3) 
    niter <- niter + 1
    k <- k + 1
    
    print(paste("error = ", error, " niter = ", niter, sep = ""))
    }
    
    return(t)
  }
  
  
# Simpson
  mysimpson <- function(f, a, b, start=1, E=10^-8)
  { 
    niter = 0 ; maxiter = 100 ; error = 1 
    
    # initial value  
    S = (b-a)*f(a)
    k <- start
    
    while (error >= E & niter <= maxiter)
    { S_0 <- S
    
    n <- 2^k ; h <- (b-a)/n
    i <- 1:(n/2)
    S <- (h/3) * sum( f(a+(2*i-2)*h) + 4*f(a+(2*i-1)*h) + f(a+(2*i)*h) )
    
    error <- abs(S - S_0) / abs(S_0 + 10^-3) 
    niter <- niter + 1
    k <- k + 1
    
    print(paste("error = ", error, " niter = ", niter, sep = ""))
    }
    
    return(S)
  }
  

# 5.3 Bayesian -----
# prior : mu ~ Cauchy(5, 2)
# likelihood : xbar ~ N(mu, 3^2/7)
  data1 <- c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)

# a) -----
  prop <- function(mu) { 
    prior <- dcauchy(mu, 5, 2)
    likelihood <- dnorm(mean(data1), mu, 3/sqrt(7))
    
    return(prior*likelihood) }

# true  
  k <- 1/integrate(prop, -Inf, Inf)$value ; k

# numerical int.    
  k1 <- 1/myrieman(prop, -100, 100) # system.time( for (i in 1:100) myrieman(prop, -100, 100) )
  k2 <- 1/mytrapezoid(prop, -100, 100) # system.time( for (i in 1:100) mytrapezoid(prop, -100, 100) )
  k3 <- 1/mysimpson(prop, -100, 100) # system.time( for (i in 1:100) mysimpson(prop, -100, 100) )
  round(c(k1, k2, k3), 5)

  
# b) -----
  post <- function(mu) { prop(mu) * 7.84654 }

  b <- integrate(post, 2, 8)$value
  b1 <- myrieman(post, 2, 8, E = 0.0001) # system.time( for (i in 1:100) myrieman(post, 2, 8, E = 0.0001) )
  b2 <- mytrapezoid(post, 2, 8, E = 0.0001) # system.time( for (i in 1:100) mytrapezoid(post, 2, 8, E = 0.0001) )
  b3 <- mysimpson(post, 2, 8, E = 0.0001) # system.time( for (i in 1:100) mysimpson(post, 2, 8, E = 0.0001) ) 
  round(c(b1, b2, b3), 5)
  
  b11 <- myrieman(post, 2, 8, E = 0.00001) # system.time( for (i in 1:100) myrieman(post, 2, 8, E = 0.00001) )
  b21 <- mytrapezoid(post, 2, 8, E = 0.00001) # system.time( for (i in 1:100) mytrapezoid(post, 2, 8, E = 0.00001) )
  
  round(c(b11, b21), 5)
  
  
# c) -----
  u <- function(mu) { exp(mu) / (1 + exp(mu)) }
  transf <- function(u) { post(log(u/(1-u))) * 1/(u*(1-u)) }
  
  c1 <- myrieman(transf, u(3), 1-1e-08) # system.time( for (i in 1:100) myrieman(transf, u(3), 1-1e-08) )
  c2 <- mytrapezoid(transf, u(3), 1-1e-08) # system.time( for (i in 1:100) mytrapezoid(transf, u(3), 1-1e-08) )
  c3 <- mysimpson(transf, u(3), 1-1e-08) # system.time( for (i in 1:100) mysimpson(transf, u(3), 1-1e-08) )
  round(c(c1, c2, c3), 5)
  integrate(transf, u(3), 1)
      
# d) -----
  u2 <- function(mu) { 1/mu }
  transf2 <- function(u) { post(1/u) * (-1/u^2) }
  
  d1 <- myrieman(transf2, u2(3), 0) # system.time( for (i in 1:100) myrieman(transf2, u2(3), 0) )
  d2 <- mytrapezoid(transf2, u2(3), 0+1e-08) # system.time( for (i in 1:100) mytrapezoid(transf2, u2(3), 0+1e-08) ) 
  d3 <- mysimpson(transf2, u2(3), 0+1e-08) # system.time( for (i in 1:100) mysimpson(transf2, u2(3), 0+1e-08) )
  round(c(d1, d2, d3), 5)

  
# 5.4 Romberg -----
  myromberg <- function(f, a, b, k) 
  { n <- 2^k ; h <- (b-a)/n
    i <- 1:(n-1)
    t <- h * ( f(a)/2 + ifelse(k==0, 0, sum(f(a+i*h)))  + f(b)/2 )
    return(t)
  }
  
  a = 10
  f <- function(x) { 1/x }
  
  m = 6
  result <- matrix(nrow = m+1, ncol = m+1)
  for (i in 0:m) { result[i+1, 1] <- myromberg(f, (a-1)/a, a-1, i) }
  for (i in 1:m) { for (j in 1:i) { result[i+1, j+1] <- (4^i*result[i+1, j] - result[i, j]) / (4^i - 1) } }
  
  c(result[m+1, m+1], log(a))
  
  Q <- matrix(0, nrow = m+1, ncol = m+1)
  for (i in 1:m) for (j in 1:j) Q[i+1, j+1] <- (result[i+1, j+1] - result[i, j+1]) / (result[i+2, j+1] - result[i+1, j+1])
   
  
  
# 6.1 Gamma -----
  tfunc <- function(y, r) { a = r-1/3 ; b = 1/sqrt(9*a) ; t = a*(1+b*y)^3 ; return(t) }
  qfunc <- function(y, r) { a = r-1/3 ; q = exp( a*log(tfunc(y, r)/a) - tfunc(y, r) + a ) ; return(q) }
  efunc <- function(y) { exp(-y^2/2) }
  
  n = 10000 ; r = 1
  set.seed(0) ; z <- rnorm(n) ; u <- runif(n) 
  
  x <- tfunc(z, r) ; u <- u[x > 0] ; z <- z[x > 0] ; x <- x[x > 0]
  
  accept <- (u <= qfunc(z, r)/efunc(z)) ; sum(accept)/n 
  z.keep <- z[accept]
  x.keep <- tfunc(z.keep, r)
  
  set.seed(0) ; qq <- data.frame(x = sort(x.keep), qgamma = sort(rgamma(length(x.keep), r, 1)))
  ggplot(qq) + theme_test() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    geom_smooth(aes(x, qgamma), method = "lm", size = 2, se = F, color = "pink") + geom_point(aes(x, qgamma), size = 1) + 
    labs(title = "Gamma Q-Q Plot", subtitle = paste("r = ", r, " ; accentance rate = ", sum(accept)/n*100, "%", sep = ""))

  gg <- data.frame(x, qx = qfunc(x, r), ex = efunc(x))
  ggplot(filter(gg, x <= 5)) + theme_test() + 
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    geom_area(aes(x, ex), fill = "gray") + geom_area(aes(x, qx), fill = "white") + 
    geom_line(aes(x, ex), color = "gray") + geom_line(aes(x, qx), linetype = "dashed") +
    labs(y = "density", title = "q(x) and e(x)")
  
  
# 6.2 Bayesian Sampling -----
  data2 <- c(8, 3, 4, 3, 1, 7, 2, 6, 2, 7)
  Lfunc <- function(lambda) 
  { lik = c()
    for (i in 1:length(lambda)) { lik[i] <- prod(dpois(data2, lambda[i])) }
    return(lik) 
  }
  
  n = 10000
  set.seed(0) ; l <- rlnorm(n, log(4), 0.5) ; u <- runif(n)
  
  accept <- (u <= Lfunc(l)/Lfunc(4.3)) ; sum(accept)/n
  
  gg <- data.frame(l, prior = dlnorm(l, log(4), 0.5)) %>% mutate(envelope = prior*Lfunc(4.3)*10^11)
  for (i in 1:nrow(gg)) { gg[i, "target"] = gg[i, "prior"]*Lfunc(gg[i, "l"])*10^11 }
  ggplot(filter(gg, l <= 20)) + 
    theme_test() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "") +
    geom_area(aes(l, envelope, fill = "envelope")) + geom_area(aes(l, target, fill = "target")) +
    geom_line(aes(l, envelope), size = 1) + geom_line(aes(l, target), size = 1, linetype = "dashed") +
    scale_fill_manual("", values = c(envelope = "gray", target = "white")) +
    labs(title = "posterior & envelope", x = "lambda", y = "density",
         subtitle = paste("acceptance rate = ", sum(accept)/n*100, "%", sep = ""))
  
  