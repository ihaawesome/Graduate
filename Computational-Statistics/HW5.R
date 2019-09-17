options("scipen" = 100)
library(dplyr)
library(ggplot2)

# (1) Peppered Moths (EM)
  myem5 <- function(E=10^-10)
  {
    nc = 85 ; ni = 196 ; nt = 341 ; n = nc + ni + nt
  # initial
    ncc = nc/3 ; nci = nc/3 ; nct = nc/3 ; nii = ni/2 ; nit = ni/2
    p = (2*ncc+nci+nct)/(2*n) ; q = (2*nii+nci+nit)/(2*n) ; r = 1 - p - q
    
    niter = 0 ; maxiter = 1000
    error = 1
    
    while (error >= E & niter <= maxiter)
    { ncc_0 <- ncc ; nci_0 <-nci ; nct_0 <-nct ; nii_0 <- nii ; nit_0 <- nit
      p_0 <- p ; q_0 <- q ; r_0 <- r
    # E-step
      ncc <- nc*p^2 / (p^2 + 2*p*q + 2*p*r)
      nci <- nc*2*p*q / (p^2 + 2*p*q + 2*p*r)
      nct <- nc*2*p*r / (p^2 + 2*p*q + 2*p*r)
      nii <- ni*q^2 / (q^2 + 2*q*r)
      nit <- ni*2*q*r / (q^2 + 2*q*r)
    # M-step
      p <- (2*ncc + nci + nct) / (2*n) ; q <- (2*nii + nci + nit) / (2*n) ; r <- 1 - p - q
      
      error <- sqrt( sum( (c(p, q, r) - c(p_0, q_0, r_0))^2 ) ) / sqrt( sum(c(p_0, q_0, r_0)^2) )
      niter <- niter + 1
      
      print(paste("error = ", error, "niter = ", niter, sep = ""))
    }
    
    loglik <- nc*log(p^2 + 2*p*q + 2*p*r) + ni*log(q^2 + 2*q*r) + 2*nt*log(r)
    output <- list(MLE = c(p, q, r), loglik = loglik, freq = round(c(ncc, nci, nct, nii, nit)))
    return(output)  
  }

  myem5() ; system.time( for (i in 1:100) myem5() )
  
  
# (2) Numerical Integration
  
# (2-1) Rieman  
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
      # if (k < 25) { k <- k + 1 } else { error <- 0 ; print("k is too large!") }
      
      print(paste("error = ", error, " niter = ", niter, sep = ""))
    }
    
    return(R)
  }
  
  
# (2-2) Trapezoidal
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
  
  
# (2-3) Simpson
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


# test
  myrieman(function(x) -x^3+1, -1, 1) ; system.time ( for (i in 1:3) myrieman(function(x) -x^3+1, -1, 1) )
  mytrapezoid(function(x) -x^3+1, -1, 1) ; system.time ( for (i in 1:1000) mytrapezoid(function(x) -x^3+1, -1, 1) )
  mysimpson(function(x) -x^3+1, -1, 1) ; system.time ( for (i in 1:1000) mysimpson(function(x) -x^3+1, -1, 1) )
  integrate(function(x) -x^3+1, lower = -1, upper = 1) ; system.time( for (i in 1:1000) integrate(function(x) -x^3+1, lower = -1, upper = 1) )
  
# normal density  
  myfunc1 <- function(x) dnorm(x, 1, 3)

  grim <- data.frame(x = seq(-8, 10, length = 1000)) %>% mutate(y = myfunc1(x))
  ggplot(grim) + theme_light() + geom_line(aes(x, y), size = 1) +
    geom_area(data = filter(grim, (x>=1 & x<=5)), aes(x, y), fill = "lightblue", alpha = 0.5) +
    geom_vline(aes(xintercept = 1), linetype = "dashed") + geom_vline(aes(xintercept = 5), linetype = "dashed")
    
  v11 <- myrieman(myfunc1, 1, 5) # ; system.time ( for (i in 1:3) myrieman(myfunc1, 1, 5) ) 
  v12 <- mytrapezoid(myfunc1, 1, 5) # ; system.time ( for (i in 1:1000) mytrapezoid(myfunc1, 1, 5) ) 
  v13 <- mysimpson(myfunc1, 1, 5) # ; system.time ( for (i in 1:1000) mysimpson(myfunc1, 1, 5) ) 
  v10 <- integrate(myfunc1, lower = 1, upper = 5)$value  # ; system.time ( for (i in 1:1000) integrate(myfunc1, 1, 5) )
  
# 험악한 함수
  myfunc2 <- function(x) ( log(1+x) + 5*x ) / ( 1 + x^2 )  
  x <- seq(5, 20, length = 1000)
  plot(x, myfunc2(x), type = "l", lwd = 2)
  
  grim <- data.frame(x = seq(0, 5, length = 1000)) %>% mutate(y = myfunc2(x))
  ggplot(grim) + theme_light() + geom_line(aes(x, y), size = 1) +
    geom_area(data = filter(grim, (x>=0.5 & x<=2.5)), aes(x, y), fill = "pink", alpha = 0.7) +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed") + geom_vline(aes(xintercept = 2.5), linetype = "dashed")
  
  v21 <- myrieman(myfunc2, 0.5, 2.5)
  v22 <- mytrapezoid(myfunc2, 0.5, 2.5)
  v23 <- mysimpson(myfunc2, 0.5, 2.5)
  v20 <- integrate(myfunc2, lower = 0.5, upper = 2.5)$value
  
# 더 험악  
  myfunc3 <- function(x) ( log(1+x) + 5*cos(x) ) / sqrt( 1 + x^2 )  
  x <- seq(5, 20, length = 1000)
  plot(x, myfunc3(x), type = "l", lwd = 2)

  grim <- data.frame(x = seq(5, 20, length = 1000)) %>% mutate(y = myfunc3(x))
  ggplot(grim) + theme_light() + geom_line(aes(x, y), size = 1) +
    geom_area(data = filter(grim, (x>=7 & x<=16)), aes(x, y), fill = "orange", alpha = 0.5) +
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 7), linetype = "dashed") + geom_vline(aes(xintercept = 16), linetype = "dashed")
  
  v31 <- myrieman(myfunc3, 7, 16, E=0.00000002) ; system.time( myrieman(myfunc3, 7, 16, E=0.00000002) )
  v32 <- mytrapezoid(myfunc3, 7, 16) ; system.time( for (i in 1:1000) mytrapezoid(myfunc3, 7, 16) )
  v33 <- mysimpson(myfunc3, 7, 16) ; system.time( for (i in 1:1000) mysimpson(myfunc3, 7, 16) )
  v30 <- integrate(myfunc3, 7, 16)$value ; system.time( for (i in 1:1000) integrate(myfunc3, 7, 16) )
  
  v41 <- myrieman(function(x) myfunc3(x)*myfunc1(x), -1+10^-8, 20, E=10^-6) 
  v42 <- mytrapezoid(function(x) myfunc3(x)*myfunc1(x), -1+10^-8, 20, E=10^-6)
  v43 <- mysimpson(function(x) myfunc3(x)*myfunc1(x), -1+10^-8, 20, E=10^-6) 
  v40 <- integrate(function(x) myfunc3(x)*myfunc1(x), -1, Inf)$value 
  
