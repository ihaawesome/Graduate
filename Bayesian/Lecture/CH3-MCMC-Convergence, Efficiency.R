library(coda)

##### CH3. MCMC의 수렴진단, 정확도, 효율 #####

##### 3.1 MCMC의 수렴진단 #####
##### Ex 3.1
n <- 20 ; xbar <- 4 ; sigma <- 1
mu0 <- 0 ; sigma0 <- 10
dataList <- list(n = n, xbar = xbar, sigma = sigma, mu0 = mu0, sigma0 = sigma0)

# Function to Compute Posterior Kernel
post.kernel <- function(mu, dataList) {
  n <- dataList$n
  xbar <- dataList$xbar
  sigma <- dataList$sigma
  mu0 <- dataList$mu0
  sigma0 <- dataList$sigma0
  post.kernel <- exp( -0.5*(((xbar-mu)*sqrt(n)/sigma)**2 + ((mu-mu0)/sigma0)**2 ))
  return(post.kernel)
}

# Function to Perform Random Walk Metropolis 
Metro <- function(nsim, mu.init, delta, dataList) {
  mu.samples <- mu.init
  mu.curr <- mu.init
  
  for (iter in 2:nsim) {
    mu.prop <- rnorm(1, mu.curr, delta)
    alpha <- min(1, post.kernel(mu.prop, dataList) /
                    post.kernel(mu.curr, dataList)) # 정규분포 q 동일해서 약분됨
    mu.next <- ifelse(runif(1) < alpha, mu.prop, mu.curr)
    
    mu.samples <- rbind(mu.samples, mu.next)
    mu.curr <- mu.next
  }
  return(mu.samples)
}

# Try Multi-Chains
delta <- 0.2
nsim <- 10000 ; nburn <- 500 ; nchains <- 3

mu.Samples <- matrix(nrow = nsim, ncol = nchains)
for (i in 1:nchains) {
  mu.init <- rnorm(1, mu0, 2) # 분산을 좀 크게 준다 (무정보)
  mu.Samples[,i] <- Metro(nsim, mu.init, delta, dataList)
}


# Fig 3.1 3개 체인의 경로그림
mycol <- c('gray50', 'steelblue', 'pink')
plot(mu.Samples[-(1:nburn),1], xlab = 'iteration', ylab = 'mu', type = 'l', main = '', col = mycol[1])
for (i in 2:3) lines(mu.Samples[-(1:nburn),i], col = mycol[i]) 

# Fig 3.2 사후밀도함수 추정
par(mfrow = c(1,2))
# 준비단계
plot(density(mu.Samples[(1:nburn),1]), xlab = 'mu', ylab = 'posterior', main = '', col = mycol[1])
for (i in 2:3) lines(density(mu.Samples[(1:nburn),i]), col = mycol[i]) 
# 준비단계 이후
plot(density(mu.Samples[-(1:nburn),1]), xlab = 'mu', ylab = 'posterior', main = '', col = mycol[1])
for (i in 2:3) lines(density(mu.Samples[-(1:nburn),i]), col = mycol[i]) 

# Fig 3.3 Gelman Plot
# 준비단계
mu.codaSamples <- mcmc(mu.Samples[(1:nburn),])
mu.codaSamples <- mcmc.list(list(mu.codaSamples[,1], mu.codaSamples[,2], mu.codaSamples[,3]))
gelman.plot(mu.codaSamples)
gelman.diag(mu.codaSamples)
# 준비단계 이후 
mu.codaSamples <- mcmc(mu.Samples[-(1:nburn),])
mu.codaSamples <- mcmc.list(list(mu.codaSamples[,1], mu.codaSamples[,2], mu.codaSamples[,3]))
gelman.plot(mu.codaSamples)
gelman.diag(mu.codaSamples)


##### 3.2 MCMC 추정치의 오차 추정 #####
##### Ex 3.1 (이어서)
par(mfrow = c(1,1))
ACF <- acf(mu.Samples[-(1:nburn),1], main = '')$acf

which(ACF < 0.05)
(ESS <- (nsim-nburn) / (1 + 2*sum(ACF[2:15]))) 
# 실제 표본수는 9500개지만 효용표본수는 1050개

# 사후추정치의 변동성
(SD <- sd(mu.Samples[-(1:nburn),1]))
(MCSE <- SD / sqrt(ESS))
