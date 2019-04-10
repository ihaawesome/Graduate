setwd("D:/STAT/CS/HW8")
library(dplyr)
library(ggplot2)

# 1. Option Pricing -----
  S0 = 100 ; K = 102 ; sig = 0.3 ; N = 50 ; r = 0.05 ; n = 100000

# analytic sol. of theta  
  c3 = 1 + 1/N
  c2 = sig * ( c3*N/1095 * (1+1/(2*N)) )^0.5
  c1 = 1/c2 * ( log(S0/K) + (c3*N/730)*(r-sig^2/2) + c3*sig^2*N/1095*(1+1/(2*N)) )
  theta = S0 * pnorm(c1) * exp( -N*(r+c3*sig^2/6)*(1-(1/N))/730 ) - K * pnorm(c1-c2) * exp(-r*N/365)
      
# simul.
  mycomp1.1 <- function(S0, K, sig, N, r, n=100000) {
    St = matrix(nrow = n, ncol = N) ; A = numeric(n) ; A.geo = numeric(n)
    for (i in 1:n) {
      Zt <- rnorm(N)
      St[i,] <- S0 * cumprod(exp((r-sig^2/2)/365 + sig*Zt/sqrt(365)))
      A[i] <- exp(-r*N/365) * max(0, mean(St[i,])-K) 
      A.geo[i] <- exp(-r*N/365) * max(0, prod(St[i,])^(1/N)-K) 
    }
    return(c(mu.mc = mean(A), theta.mc = mean(A.geo)))
  }
  
  nsim = 100
  cal1 <- data.frame(mu.mc = numeric(nsim), theta.mc = numeric(nsim)) 
  for (i in 1:nsim) { cal1[i,] <- mycomp1.1(100, 102, 0.3, 50, 0.05) }
  cal1 <- cal1 %>% mutate(mu.cv1 = mu.mc - (theta.mc - theta))
  cal1 <- cal1 %>% mutate(mu.cv2 = mu.mc - cov(mu.mc, theta.mc)/var(theta.mc) * (theta.mc - theta))
  apply(cal1, 2, mean) ; apply(cal1, 2, sd)   
  cor(cal1$mu.mc, cal1$theta.mc) # 0.9999
  
# 6.3 ------
  rm(list = ls())
  
  qfunc <- function(x) { exp(-(abs(x))^3/3) }
  const <- integrate(qfunc, -Inf, Inf)$value
  q <- seq(-4, 4, length = 100)

  plot(q, 3*dnorm(q), type ="l", col = "blue") ; lines(q, qfunc(q), type = "l")  
  
# a) IS w/ std.weights
  n = 10000
  a <- rnorm(n) ; usiw <- qfunc(a)/(dnorm(a)) ; siw <- usiw/sum(usiw)
  mu.is <- sum(a^2*siw)
  
  ggplot() + scale_fill_brewer("", palette = "Set1") + theme_test() + 
    labs(x = "y", y = "density", title = "Importance Sampling") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
    geom_area(aes(a, dnorm(a), fill = "ISF"), alpha = 0.3, color = "black") + 
    geom_area(aes(a, qfunc(a)/const, fill = "target"), alpha = 0.3, color = "black")
  
  for (i in 1:1000) {
    a <- rnorm(n) ; usiw <- qfunc(a)/(dnorm(a)) ; siw <- usiw/sum(usiw)
    mu.is[i] <- sum(a^2*siw)
  }
  mean(mu.is) ; sd(mu.is)
  
  
# b) RS
  mycomp2.2 <- function(n = 10000) {
    y <- rnorm(n) ; u <- runif(n)
    keep <- (u <= qfunc(y)/(3*dnorm(y))) ; x <- y[keep]
    return(x)
  }
  b <- mycomp2.2() ; mu.rs <- mean(b^2)
  
  ggplot() + scale_fill_manual("", values = c("e(y)" = "pink", "q(y)" = "white")) + theme_test() + 
    labs(x = "y", y = "propotional density", title = "Rejection Sampling") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
    geom_area(aes(a, 3*dnorm(a), fill = "e(y)"), color = "black") + 
    geom_area(aes(a, qfunc(a), fill = "q(y)"), color = "black")
  
# c) Philippe and Robert
  mycomp2.3 <- function(x.rs) {
    c1 <- sum(diff(sort(x.rs))*sort(x.rs)[-length(x.rs)]^2*qfunc(sort(x.rs)[-length(x.rs)]))
    c2 <- sum(diff(sort(x.rs))*qfunc(sort(x.rs)[-length(x.rs)]))
    return(c1/c2)
  } 
  mu.pnr <- mycomp2.3(b)
  
# d)
  nsim = 1000 ; cal2.4 <- data.frame(a = numeric(nsim), b = numeric(nsim), c = numeric(nsim))
  
  for (i in 1:nsim) { 
    x <- rnorm(n) ; usiw <- qfunc(x)/(dnorm(x)) ; siw <- usiw/sum(usiw)
    x.rs <- mycomp2.2() ; cal2.4[i,] <- c(sum(x^2*siw), mean(x.rs^2), mycomp2.3(x.rs)) } 
  
  apply(cal2.4, 2, mean) ; apply(cal2.4, 2, sd)
  system.time( for (i in 1:nsim) { x <- rnorm(n) ; usiw <- qfunc(x)/(dnorm(x)) ; siw <- usiw/sum(usiw) ; sum(x^2*siw) } )
  system.time( for (i in 1:nsim) { x.rs <- mycomp2.2() ; mean(x.rs^2) } ) 
  system.time( for (i in 1:nsim) mycomp2.3(x.rs) )
  
  
# 6.6 -----
  rm(list = ls())
# a. size
  myh <- function(z) { ifelse(z >= 1.645, 1, 0) }
  
  mycomp3.1 <- function(n = 1000, nsim = 1000) {
    lambda.is <- 2 + 1.645*sqrt(2/25)
    
    cal <- data.frame(matrix(numeric(nsim)*5, ncol = 5))
    colnames(cal) <- c("MC", "AS", "IS.ustd", "IS.std", "wbar")
    
    for (i in 1:nsim) {
      x <- matrix(rpois(25*n, 2), ncol = 25)
      z <- apply(x, 1, function(x) (mean(x)-2)/sqrt(2/25))
      
      mu1 <- mean( myh(z) ) # naive MC
      mu2 <- sum( myh(z[1:(n/2)]) + myh(-z[1:(n/2)]) ) / n # AS
      
      x.is <- matrix(rpois(25*n, lambda.is), ncol = 25)
      xbar.is <- apply(x.is, 1, mean)
      z.is <- apply(x.is, 1, function(x) (mean(x)-2)/sqrt(2/25))
      
      usiw <- dnorm(xbar.is, 2, sqrt(2/25)) / dnorm(xbar.is, lambda.is, sqrt(lambda1/25))
      siw <- usiw / sum(usiw)
      
      mu3.1 <- mean( myh(z.is)*usiw ) # IS w/ unstd.weights
      mu3.2 <- sum( myh(z.is)*siw ) # IS w/ std.weights
      
      cal[i,] <- c(mu1, mu2, mu3.1, mu3.2, mean(usiw))
    }
    cal <- cal %>% transmute(MC, AS, IS.ustd, IS.std, 
                             ISCV = IS.ustd - cov(IS.ustd, wbar)/var(wbar)*(wbar-1), wbar)
    return(cal)
  }
  
  table3.1 <- mycomp3.1()
  info3.1 <- data.frame(mean = apply(table3.1, 2, mean), sd = apply(table3.1, 2, sd), 
                        lower = apply(table3.1, 2, function(x) sort(x)[25]), 
                        upper = apply(table3.1, 2, function(x) sort(x)[975]))[-6,]
  
  
  ggplot(table3.1) + theme_test() + scale_color_brewer("", palette = "Set1") +
    geom_density(aes(MC, color = "MC")) + geom_density(aes(AS, color = "AS")) + 
    geom_density(aes(IS.ustd, color = "IS*")) + geom_density(aes(IS.std, color = "IS")) + geom_density(aes(ISCV, color = "ISCV"))
 
# b. power
  mycomp3.2 <- function(n = 1000, nsim = 1000, lambda) {
    b <- 2 + 1.645*sqrt(2/25)
    lambda.is <- lambda + 1.645*sqrt(lambda/25)
    
    cal <- data.frame(matrix(numeric(nsim)*5, ncol = 5))
    colnames(cal) <- c("MC", "AS", "IS.ustd", "IS.std", "wbar")
    
    for (i in 1:nsim) {
      
      x <- matrix(rpois(n*25, lambda), ncol = 25) # under H1: lambda = 2.2 > 2
      z <- apply(x, 1, function(x) (b-mean(x))/(sd(x)/5))
  
      mu1 <- mean( 1 - pnorm(z) )
      mu2 <- sum( 1 - pnorm(z[1:(n/2)]) + pnorm(-z[1:(n/2)]) ) / n # AS
      
      x.is <- matrix(rpois(n*25, lambda.is), ncol = 25)
      xbar.is <- apply(x.is, 1, mean)
      z.is <- apply(x.is, 1, function(x) (b-mean(x))/(sd(x)/5))
      
      usiw <- dnorm(xbar.is, lambda, sqrt(lambda/25)) / dnorm(xbar.is, lambda.is, sqrt(lambda.is/25))
      siw <- usiw / sum(usiw)
      
      mu3.1 <- mean( (1-pnorm(z.is))*usiw ) # IS w/ unstd.weights
      mu3.2 <- sum( (1-pnorm(z.is))*siw ) # IS w/ std.weights
      
      cal[i,] <- c(mu1, mu2, mu3.1, mu3.2, mean(usiw))
      if (i %% 10 == 0) print(i)
    }
    cal <- cal %>% transmute(MC, AS, IS.ustd, IS.std, 
                             ISCV = IS.ustd - cov(IS.ustd, wbar)/var(wbar)*(wbar-1), wbar)
    return(cal)
  }
  
  mylambda <- list() ; l <- seq(2.2, 4, length = 5)
  for (i in 1:5) { mylambda[[i]] <- mycomp3.2(lambda = l[i]) } 

  table3.2 <- data.frame(mylambda[[1]], mylambda[[2]], mylambda[[3]], mylambda[[4]], mylambda[[5]]) %>% select(-6*(1:5))
  table3.2 <- table3.2[c(seq(1, 21, by = 5), seq(2, 22, by = 5), seq(3, 23, by = 5), seq(4, 24, by = 5), seq(5, 25, by = 5))]
  colnames(table3.2) <- c(paste("MC", 1:5, sep = ""), paste("AS", 1:5, sep = ""),
                          paste("IS.ustd", 1:5, sep = ""), paste("IS.std", 1:5, sep = ""), paste("ISCV", 1:5, sep = ""))

  table3.2.1 <- data.frame(lambda = rep(l, 5), 
                           mean = apply(table3.2, 2, mean), sd = apply(table3.2, 2, sd), 
                           lower = apply(table3.2, 2, function(x) sort(x)[25]), 
                           upper = apply(table3.2, 2, function(x) sort(x)[975]))

  ggplot(table3.2.1[11:15,]) + geom_line(aes(lambda, mean), color = "pink", size = 1) + 
    geom_point(aes(lambda, mean)) + geom_errorbar(aes(lambda, ymin = lower, ymax = upper), size = 0.5, width = 0.02) + 
    labs(y = "power", title = "ISCV") + theme_test() + theme(plot.title = element_text(hjust = 0.5))

  
# 6.7 -----
  rm(list = ls())
  S0 = 50 ; K = 52 ; sig = 0.5 ; N = 30 ; r = 0.05

# a. european
  n = 10000 ; nsim = 1000 ; C.mc <- c()
  for (i in 1:nsim) {  
    z <- rnorm(n)
    ST <- S0 * exp( (r-sig^2/2)*N/365 + sig*z*sqrt(N/365) )
    C <- exp(-r*N/365) * ifelse(ST>=K, ST-K, 0)  
    C.mc[i] <- mean(C)
  }
  mean(C.mc) ; sd(C.mc)
  ggplot() + geom_density(aes(C.mc), fill = "pink", alpha = 0.5) + theme_test()
  
# b. asian
  A.cv <- data.frame(mu.mc = numeric(nsim), theta.mc = numeric(nsim))
  for (i in 1:nsim) { A.cv[i,] <- mycomp1.1(50, 52, 0.5, 30, 0.05) }
  
  c3 = 1 + 1/N
  c2 = sig * ( c3*N/1095 * (1+1/(2*N)) )^0.5
  c1 = 1/c2 * ( log(S0/K) + (c3*N/730)*(r-sig^2/2) + c3*sig^2*N/1095*(1+1/(2*N)) )
  theta = S0 * pnorm(c1) * exp( -N*(r+c3*sig^2/6)*(1-(1/N))/730 ) - K * pnorm(c1-c2) * exp(-r*N/365)
  
  A.cv <- A.cv %>% mutate(mu.cv1 = mu.mc - (theta.mc - theta),
                          mu.cv2 = mu.mc - cov(mu.mc, theta.mc)/var(theta.mc) * (theta.mc - theta))
  apply(A.cv, 2, mean) ; apply(A.cv, 2, sd)   

# c & d.    
  mycomp4 <- function(S0, K, sig, N, r, n = 10000) {
    St = matrix(nrow = n, ncol = N) ; St.as = matrix(nrow = n, ncol = N)
    A.mean = numeric(n) ; A.geo = numeric(n) ; A.as = numeric((n/2))
      for (i in 1:(n/2)) {
        Zt <- rnorm(N)
        St.as[i,] <- S0 * cumprod(exp((r-sig^2/2)/365 + sig*Zt/sqrt(365)))
        St.as[n/2+i,] <- S0 * cumprod(exp((r-sig^2/2)/365 - sig*Zt/sqrt(365)))
      }
      for (i in 1:n) { A.as[i] <- exp(-r*N/365) * max(0, mean(St.as[i,])-K) }
    return(mean(A.as))
  }
  
  nsim = 1000
  table4 <- data.frame(mu.as = numeric(nsim), mu.mc = numeric(nsim), theta.mc = numeric(nsim))
  for (i in 1:nsim) {
    table4[i, 1] <- mycomp4(50, 52, 0.5, 30, 0.05, n = 10000)
    table4[i, 2:3] <- mycomp1.1(50, 52, 0.5, 30, 0.05, n = 10000)
    print(i)
  }  
  table4 <- table4 %>% mutate(mu.cv1 = mu.mc - (theta.mc - theta), 
                              mu.cv2 = mu.mc - cov(mu.mc, theta.mc)/var(mu.mc) * (theta.mc - theta))
  apply(table4, 2, mean) ; apply(table4, 2, sd)
  
  ggplot(table4) + theme_test() + scale_fill_brewer("method", palette = "Blues", direction = -1) + 
    scale_linetype_manual("method", values = c(MC = 1, AS = 2, CV = 3)) +
    labs(x = "estimate") + theme(legend.position = "top") +
    geom_density(aes(mu.mc, linetype  = "MC", fill = "MC"), alpha = 0.5) + 
    geom_density(aes(mu.cv1, linetype  = "CV", fill = "CV"), alpha = 0.5) +
    geom_density(aes(mu.as, linetype  = "AS", fill = "AS"), alpha = 0.5)
    