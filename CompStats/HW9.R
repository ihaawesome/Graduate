library(dplyr)
library(ggplot2)

##### 1. Example 7.2 Mixture Dist #####
  data1 <- read.table("C:/Users/HG/Desktop/18-2/CS/Textbook/Datasets/mixture.dat", header = T) ; y <- data1$y
#  f <- function(d) { d*dnorm(y, 7, 0.5) + (1-d)*dnorm(y, 10, 0.5) }
  L <- function(d) { f = d*dnorm(y, 7, 0.5) + (1-d)*dnorm(y, 10, 0.5) ; return(prod(f)) }
  
  n = 100000
# 1-1. independent chain    
  mychain1 <- function(a, b, n = 100000) {
    d = c(runif(1, 0, 1), numeric(n-1))
    for (i in 2:n) {
      d_ <- rbeta(1, a, b) ; R <- L(d_)/L(d[i-1])
      d[i] <- ifelse(rbinom(1, 1, min(1,R)) == 1, d_, d[i-1]) }
    return(d)
  }
  
# 1) Prior: Beta(1, 1) = Unif(0, 1)
  d1 <- mychain1(1, 1)
  ggplot() + geom_line(aes(1:n, d1), color = "darkblue") + theme_test() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "iter", y = "delta", title = "Beta(1,1)")
  ggplot() + theme_test() + theme(plot.title = element_text(hjust = 0.5)) +
    geom_histogram(aes(d1[-(1:500)]), fill = "darkblue", color = "white") + labs(x = "delta", title = "Beta(1,1)")
  
# 2) Prior: Beta(2, 10) 
  d2 <- mychain1(2, 10)
  ggplot() + geom_line(aes(1:n, d2), color = "darkred") + theme_test() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "iter", y = "delta", title = "Beta(2,10)")
  ggplot() + theme_test() + theme(plot.title = element_text(hjust = 0.5)) +
    geom_histogram(aes(d2[-(1:500)]), fill = "darkred", color = "white") +
    labs(x = "delta", title = "Beta(2,10)")

  indep <- data.frame(d1, d2)[-(1:500),] 
  
  info1 <- data.frame(mean = apply(indep, 2, mean), sd = apply(indep, 2, sd),
                      lowerCI = apply(indep, 2, function(x) quantile(x, 0.025)), 
                      upperCI = apply(indep, 2, function(x) quantile(x, 0.975)))
  info1 %>% round(5)
  summary(indep)
  
# 1-2. random walk chain  
  ilogit <- function(u) { exp(u)/(1+exp(u)) }
  
  mychain2 <- function(b, n = 100000) {
    u = c(runif(1, 0, 1), numeric(n-1))
    for (i in 2:n) {
      u_ <- u[i-1] + runif(1, -b, b) 
      R <- L(ilogit(u_))/L(ilogit(u[i-1])) * abs(ilogit(u_)/(1+exp(u_)))/abs(ilogit(u[i-1])/(1+exp(u[i-1])))
      u[i] <- ifelse(rbinom(1, 1, min(1,R)) == 1, u_, u[i-1]) }
    d <- ilogit(u) ; return(d)
  }
  d3 <- mychain2(1)
  d4 <- mychain2(0.01)
  
  ggplot() + geom_line(aes(1:n, d3), color = "darkblue") + theme_test() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "iter", y = "delta", title = "Unif(-1,1)") + lims(y = c(0.4, 0.9))
  ggplot() + geom_line(aes(1:n, d4), color = "darkred") + theme_test() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "iter", y = "delta", title = "Unif(-0.01,0.01)") + lims(y = c(0.4, 0.9))

  ggplot() + theme_test() + theme(plot.title = element_text(hjust = 0.5)) +
    geom_histogram(aes(d3[-(1:500)]), fill = "darkblue", color = "white") +
    labs(x = "delta", title = "Unif(-1,1)") + lims(x = c(0.5, 0.9))
  ggplot() + theme_test() + theme(plot.title = element_text(hjust = 0.5)) +
    geom_histogram(aes(d4[-(1:500)]), fill = "darkred", color = "white") +
    labs(x = "delta", title = "Unif(-0.1,0.1)")  + lims(x = c(0.5, 0.9)) 
  
  random <- data.frame(d3, d4)[-(1:500),]
  info2 <- data.frame(mean = apply(random, 2, mean), sd = apply(random, 2, sd),
                      lowerCI = apply(random, 2, function(x) quantile(x, 0.025)), 
                      upperCI = apply(random, 2, function(x) quantile(x, 0.975)))
  info2 %>% round(5)
  rbind(summary(d3[-(1:500)]), summary(d4[-(1:500)])) %>% round(5)
  
  
  
##### 2. Problems 7.5 clicical trial #####
  rm(list = ls())
  data2 <- read.table("C:/Users/HG/Desktop/18-2/CS/Textbook/Datasets/breastcancer.dat", header = T)
  x <- data2$recurtime ; C <- data2$censored ; H <- data2$treatment
  
# a. summary & plot
  summary(x[H==1&C==0])
  summary(x[H==0&C==0])
  summary(x[H==1&C==1])
  summary(x[H==0&C==1])

  a <- ggplot(data2) + theme_test() + scale_fill_brewer("treatment", palette = "Set1") + theme(plot.title = element_text(hjust = 0.5))
  a + geom_density(aes(recurtime, group = treatment, fill = factor(treatment)), alpha = 0.5) + facet_wrap(~censored)
  
# b-c. sampling      
  gtheta <- function(tau, a=3, b=1, c=60, d=120) { 
    rgamma(1, shape = a+sum((1-C)*(1-H))+sum((1-C)*H)+1, rate = c+sum((1-H)*x)+(d+sum(H*x))*tau) }
  gtau <- function(theta, a=3, b=1, c=60, d=120) { 
    rgamma(1, shape = b+sum((1-C)*H)+1, rate = (d+sum(H*x))*theta) }
  
  mychain3 <- function(n = 10000, a=3, b=1, c=60, d=120) {
    theta <- c(gtheta(1, a, b, c, d), numeric(n-1)) ; tau <- c(gtau(theta[1], a, b, c, d), numeric(n-1))
    for (i in 2:n) { theta[i] <- gtheta(tau[i-1], a, b, c, d) ; tau[i] <- gtau(theta[i], a, b, c, d) }
    return(data.frame(theta, tau))
  }
  
  n = 100000 ; burnin = 1:500  
  gibbs1 <- mychain3(n) ; theta1 <- gibbs1$theta ; tau1 <- gibbs1$tau
  
  plt <- ggplot() + theme_test() + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_brewer("", palette = "Set1")
  plt + geom_line(aes(1:n, theta1)) + labs(x = "iter", y = "theta", title = "theta")
  plt + geom_line(aes(1:n, tau1)) + labs(x = "iter", y = "tau", title = "tau")
  
  acf(theta1[-burnin], main="theta") ; acf(tau1[-burnin], main="tau")
  plt + geom_histogram(aes(x = theta1[-burnin], y = ..density..), color = "white") + labs(x = "theta", title = "theta")
  plt + geom_histogram(aes(x = tau1[-burnin], y = ..density..),  color = "white") + labs(x = "tau", title = "tau")
  
# d.summary statistics
  gibbs1 <- gibbs1[-(1:500),] ; theta1 <- theta1[-(1:500)] ; tau1 <- tau1[-(1:500)]
  summary(gibbs1)
  summ1 <- data.frame(mean = apply(gibbs1, 2, mean), sd = apply(gibbs1, 2, sd),
                      lowerCI = apply(gibbs1, 2, function(x) sort(x)[length(x)*0.025]), 
                      upperCI = apply(gibbs1, 2, function(x) sort(x)[length(x)*0.975]))
  summ1 %>% round(5)

# e. tau plot
  a=3 ; b=1 ; c=60 ; d=120
  
  prior <- function(theta, tau, a=3, b=1, c=60, d=120) { theta^a * tau^b * exp(-c*theta-d*tau*theta) }
  L <- function(theta, tau) { 
    theta^(sum((1-C)*(1-H))+sum((1-C)*H)) * tau^(sum((1-C)*H)) * exp(-sum((1-H)*x)*theta-sum(H*x)*tau*theta) }
  post <- function(theta, tau, a=3, b=1, c=60, d=120) { prior(theta, tau, a, b, c, d) * L(theta, tau) }
  
  tau <- seq(0, 6, length = 1000)
  prior_tau <- prior(mean(theta1), tau)  / integrate(prior, theta = mean(theta1), lower = 0, upper = Inf)$value
  post_tau <- post(mean(theta1), tau) / integrate(post, theta = mean(theta1), lower = 0, upper = Inf)$value
  
  plt + scale_color_brewer("", palette = "Set1", direction = -1) +
    geom_area(aes(tau, prior_tau, fill = "prior"), color = "black", alpha = 0.3) + 
    geom_area(aes(tau, post_tau, fill = "posterior"), color = "black", alpha = 0.3) +
    labs(y = "density")

# g.
# g-1. half
  gibbs2 <- mychain3(n, a=3/2, b=1/2, c=60/2, d=120/2) ; theta2 <- gibbs2$theta ; tau2 <- gibbs2$tau
  
  plt + geom_line(aes(1:n, theta2)) + labs(x = "iter", y = "theta", title = "theta")
  plt + geom_line(aes(1:n, tau2)) + labs(x = "iter", y = "tau", title = "tau")
  acf(theta2[-burnin]) ; acf(tau2[-burnin])
  plt + geom_histogram(aes(theta2[-burnin]), color = "white") + labs(x = "theta", title = "theta")
  plt + geom_histogram(aes(tau2[-burnin]), color = "white") + labs(x = "tau", title = "tau")
  
  gibbs2 <- gibbs2[-burnin,] ; theta2 <- gibbs2$theta ; tau2 <- gibbs2$tau
  summary(gibbs2)
  summ2 <- data.frame(mean = apply(gibbs2, 2, mean), sd = apply(gibbs2, 2, sd),
                      lowerCI = apply(gibbs2, 2, function(x) quantile(x, 0.025)), 
                      upperCI = apply(gibbs2, 2, function(x) quantile(x, 0.975)))
  summ2 %>% round(5)

# g-2. double  
  gibbs3 <- mychain3(n, a=3*2, b=1*2, c=60*2, d=120*2) ; theta3 <- gibbs3$theta ; tau3 <- gibbs3$tau
  plt + geom_line(aes(1:n, theta3)) + labs(x = "iter", y = "theta", title = "theta")
  plt + geom_line(aes(1:n, tau3)) + labs(x = "iter", y = "tau", title = "tau")
  acf(theta3[-burnin]) ; acf(tau3[-burnin])
  plt + geom_histogram(aes(theta3[-burnin]), color = "white") + labs(x = "theta", title = "theta")
  plt + geom_histogram(aes(tau3[-burnin]), color = "white") + labs(x = "tau", title = "tau")
  
  gibbs3 <- gibbs3[-burnin,] ; theta3 <- gibbs3$theta ; tau3 <- gibbs3$tau
  summary(gibbs3)
  summ3 <- data.frame(mean = apply(gibbs3, 2, mean), sd = apply(gibbs3, 2, sd),
                      lowerCI = apply(gibbs3, 2, function(x) quantile(x, 0.025)), 
                      upperCI = apply(gibbs3, 2, function(x) quantile(x, 0.975)))
  summ3 %>% round(5)
  
  tau <- seq(0, 5, length = 1000)
  prior_tau1 <- prior(mean(theta1), tau)  / integrate(prior, theta = mean(theta1), lower = 0, upper = Inf)$value
  post_tau1 <- post(mean(theta1), tau) / integrate(post, theta = mean(theta1), lower = 0, upper = Inf)$value
  prior_tau2 <- prior(mean(theta2), tau, a=3/2, b=1/2, c=60/2, d=120/2) / 
    integrate(prior, theta = mean(theta2), a=3/2, b=1/2, c=60/2, d=120/2, lower = 0, upper = Inf)$value
  post_tau2 <- post(mean(theta2), tau, a=3/2, b=1/2, c=60/2, d=120/2) / 
    integrate(post, theta = mean(theta2), a=3/2, b=1/2, c=60/2, d=120/2, lower = 0, upper = Inf)$value
  prior_tau3 <- prior(mean(theta3), tau, a=3*2, b=1*2, c=60*2, d=120*2) / 
    integrate(prior, theta = mean(theta3), a=3*2, b=1*2, c=60*2, d=120*2, lower = 0, upper = Inf)$value
  post_tau3 <- post(mean(theta3), tau, a=3*2, b=1*2, c=60*2, d=120*2) / 
    integrate(post, theta = mean(theta3), a=3*2, b=1*2, c=60*2, d=120*2, lower = 0, upper = Inf)$value
  
  compare <- data.frame(par = c(rep("origin", length(tau)), rep("half", length(tau)), rep("double", length(tau))), 
                        dist = c(rep("1.prior", 3*length(tau)), rep("2.posterior", 3*length(tau))),
                        density = c(prior_tau1, prior_tau2, prior_tau3, post_tau1, post_tau2, post_tau3))
  
  ggplot(compare) + theme_test() + 
    scale_linetype_discrete("") + scale_color_manual("", values = c(origin = "black", half = "blue", double = "red")) +
    geom_line(aes(rep(tau, 6), density, color = par, linetype = par)) + facet_wrap(~dist) + labs(x = "tau") + theme(legend.position = "top")
  
  
  
  
  