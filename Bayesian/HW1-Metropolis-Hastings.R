setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')
library(coda)
library(tidyverse)

# dataList
mu0 = 10 ; sigsq0 = 25 ; a = 0.5 ; b = 1
x = c(10, 13, 15, 11, 9, 18, 20, 17, 23, 21)
dataList <- list(x = x, mu0 = mu0, sigsq0 = sigsq0, a = a, b = b)

# posterior kernel
post.norm <- function(theta, dataList) {
  # retrieve data from dataList
  x = dataList$x
  mu0 = dataList$mu0
  sigsq0 = dataList$sigsq0
  a = dataList$a
  b = dataList$b
  
  # calculate likelihood
  mu <- theta[1] ; sigsq <- theta[2]
  f <- exp( -0.5*length(x)*log(sigsq) -0.5*sum((x-mu)^2)/sigsq -
             0.5*(mu-mu0)^2/sigsq0-(a+1)*log(sigsq) -b/sigsq )
  return(f)
}

# metropolis algorithm
metro.norm <- function(nsim, nburn, delta, dataList, initList) {
  # initial values of mu and log.sigsq
  mu = initList$mu
  log.sigsq = log(initList$sigsq)
  
  theta.curr <- c(mu, log.sigsq)
  p <- length(theta.curr)

  # --- start iterations
  par.samples <- matrix(0, nsim, p)
  for (iter in 1:(nsim+nburn)) {
    z <- rnorm(p)
    theta.prop <- z*delta + theta.curr
    mu.curr <- theta.curr[1] ; sigsq.curr <- exp(theta.curr[2])
    mu.prop <- theta.prop[1] ; sigsq.prop <- exp(theta.prop[2])
    alpha <- (post.norm(c(mu.prop, sigsq.prop), dataList) / post.norm(c(mu.curr, sigsq.curr), dataList)) *
             (sigsq.prop / sigsq.curr)
    if (runif(1) < alpha) { 
      theta.next <- theta.prop 
    } else { 
      theta.next <- theta.curr 
    }
    theta.curr <- theta.next
    
    if (iter > nburn) par.samples[iter-nburn,] <- c(theta.next[1], exp(theta.next[2]))
  }
  # --- end iterations
  
  return(par.samples)
}

# draw random initial values
init.random <- function(x) {
  resampledX <- sample(x, replace = T)
  mu.init <- mean(resampledX)
  sigsq.init <- var(resampledX)
  return(list(mu = mu.init, sigsq = sigsq.init))
}

# check convergence
# 경로그림과 사후밀도함수
check.plots <- function(mcmc.samples) {
  mu.samples <- mcmc.samples[,1,]
  sigsq.samples <- mcmc.samples[,2,]
  
  cols <- c('gray50', 'lightblue', 'pink')
  par(mfrow = c(2, 2))
  plot(mu.samples[,1], col = cols[1], type = "l", xlab = "iteration", ylab = "", main = quote(mu))
  lines(mu.samples[,2], col = cols[2])
  lines(mu.samples[,3], col = cols[3])
  plot(density(mu.samples[,1]), col = cols[1], xlab = "", ylab = "posterior density", main = quote(mu))
  lines(density(mu.samples[,2]), col = cols[2])
  lines(density(mu.samples[,3]), col = cols[3])
  plot(sigsq.samples[,1], col = cols[1], type = "l", xlab = "iteration", ylab = "", main = quote(sigma^2))
  lines(sigsq.samples[,2], col = cols[2])
  lines(sigsq.samples[,3], col = cols[3])
  plot(density(sigsq.samples[,1]), col = cols[1], xlab = "", ylab = "posterior density", main = quote(sigma^2))
  lines(density(sigsq.samples[,2]), col = cols[2])
  lines(density(sigsq.samples[,3]), col = cols[3])
}

# gelman 상수
check.gelman <- function(mcmc.samples) {
  samples.1 <- mcmc(mcmc.samples[,,1])
  samples.2 <- mcmc(mcmc.samples[,,2])
  samples.3 <- mcmc(mcmc.samples[,,3])
  codaSamples <- mcmc.list(list(samples.1, samples.2, samples.3))
  gelman <- gelman.diag(codaSamples)
  return(gelman)
}

# 채택확률
check.accept <- function(mcmc.samples) {
  Metro.draws <- mcmc(mcmc.samples[,,1])
  accept.rate <- 1 - rejectionRate(Metro.draws)
  return(accept.rate)
} 

# 자기상관함수
check.acf <- function(mcmc.samples) {
  mu.samples <- mcmc.samples[,1,]
  sigsq.samples <- mcmc.samples[,2,]
  
  par(mfrow = c(1, 2))
  acf(mu.samples[,1], main = quote(mu))
  acf(sigsq.samples[,1], main = quote(sigma^2))
}

# 95% HPD
infer.HPD <- function(mcmc.samples) {
  mcmc.samples.combined <- rbind(mcmc.samples[,,1], mcmc.samples[,,2], mcmc.samples[,,3])
  estimate <- apply(mcmc.samples.combined, 2, mean)
  HPD <- apply(mcmc.samples.combined, 2, function(x) quantile(x, c(0.025, 0.975)))
  out <- rbind(estimate, HPD)
  return(out)
}


##### Exercises #####

# 1) preliminary run
nchains = 3
nsim = 2000 ; nburn = 0
p = 2
mcmc.samples.one <- array(0, c(nsim, p, nchains))

delta.one <- 1
for (ich in 1:nchains) {
  initList <- init.random(x)
  mcmc.samples.one[,,ich] <- metro.norm(nsim, nburn, delta.one, dataList, initList)
}

# mu와 sigsq의 추정치와 표준편차
c(mean(mcmc.samples.one[,1,1]), mean(log(mcmc.samples.one[,2,1])))
c(var(mcmc.samples.one[,1,1]), var(log(mcmc.samples.one[,2,1])))

(delta.mu <- sd(mcmc.samples.one[,1,]))
(delta.sigsq <- sd(log(mcmc.samples.one[,2,])))


# 2)
(delta.est <- c(delta.mu, delta.sigsq) * sqrt(2.4))

# 3) metro.norm 함수에 delta 대신 delta.each를 사용한다 (vector)

# 4) quick MH & acceptance rate
mcmc.samples.est <- array(0, c(nsim, p, nchains))
for (ich in 1:nchains) {
  initList <- init.random(x)
  mcmc.samples.est[,,ich] <- metro.norm(nsim, nburn, delta.est, dataList, initList)
}
check.accept(mcmc.samples.est) 
check.plots(mcmc.samples.est)


# 5)
delta.change <- c(3.3, 1.1)
mcmc.samples.change <- array(0, c(nsim, p, nchains))
for (ich in 1:nchains) {
  initList <- init.random(x)
  mcmc.samples.change[,,ich] <- metro.norm(nsim, nburn, delta.change, dataList, initList)
}
check.accept(mcmc.samples.change) 
check.plots(mcmc.samples.change)


# 6) 
nchains = 3
nsim = 25000 ; nburn = 0

#####
mcmc.samples.change <- array(0, c(nsim, p, nchains))
for (ich in 1:nchains) {
  initList <- init.random(x)
  mcmc.samples.change[,,ich] <- metro.norm(nsim, nburn, delta.change, dataList, initList)
}

c(mean(mcmc.samples.change[-(1:5000),1,1]), mean(log(mcmc.samples.change[-(1:5000),2,1])))
c(var(mcmc.samples.change[-(1:5000),1,1]), var(log(mcmc.samples.change[-(1:5000),2,1])))

check.plots(mcmc.samples.change[-(1:5000),,])
check.gelman(mcmc.samples.change[-(1:5000),,])
check.accept(mcmc.samples.change[-(1:5000),,])
check.acf(mcmc.samples.change[-(1:5000),,])
(HPD.optim <- infer.HPD(mcmc.samples.change[-(1:5000),,]))

par(mfrow = c(1, 2))
plot(density(mcmc.samples.change[-(1:5000),1,1]), col = 'gray50', main = quote(mu), xlab = '')
abline(v = HPD.optim[1,1], col = 'steelblue')
abline(v = HPD.optim[-1,1], col = 'red', lty = 2)
plot(density((mcmc.samples.change[-(1:5000),2,1])), col = 'gray50', main = quote(sigma^2), xlab = '')
abline(v = (HPD.optim[1,2]), col = 'steelblue')
abline(v = (HPD.optim[-1,2]), col = 'red', lty = 2)

check.plots(mcmc.samples.change[(1:1000),,])
check.gelman(mcmc.samples.change[(1:1000),,])

#####
mcmc.samples.one <- array(0, c(nsim, p, nchains))
for (ich in 1:nchains) {
  initList <- init.random(x)
  mcmc.samples.one[,,ich] <- metro.norm(nsim, nburn, delta = delta.one, dataList, initList)
}
c(mean(mcmc.samples.one[-(1:5000),1,1]), mean(log(mcmc.samples.one[-(1:5000),2,1])))
c(var(mcmc.samples.one[-(1:5000),1,1]), var(log(mcmc.samples.one[-(1:5000),2,1])))

check.plots(mcmc.samples.one[-(1:5000),,])
check.gelman(mcmc.samples.one[-(1:5000),,])
check.accept(mcmc.samples.one[-(1:5000),,])
check.acf(mcmc.samples.one[-(1:5000),,])
(HPD.one <- infer.HPD(mcmc.samples.one[-(1:5000),,]))

par(mfrow = c(1, 2))
plot(density(mcmc.samples.one[-(1:5000),1,1]), col = 'gray50', main = quote(mu), xlab = '')
abline(v = HPD.one[1,1], col = 'steelblue')
abline(v = HPD.one[-1,1], col = 'red', lty = 2)
plot(density(mcmc.samples.one[-(1:5000),2,1]), col = 'gray50', main = quote(sigma^2), xlab = '')
abline(v = HPD.one[1,2], col = 'steelblue')
abline(v = HPD.one[-1,2], col = 'red', lty = 2)

check.plots(mcmc.samples.one[(1:1000),,])
check.gelman(mcmc.samples.one[(1:1000),,])

