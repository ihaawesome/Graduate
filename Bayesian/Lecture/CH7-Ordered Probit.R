setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')
library(rjags)
library(mvtnorm)
library(truncnorm)

##### CH7. 순서적 프로빗 모형의 베이지안 추론 #####
##### 7.1 순서적 프로빗 모형 #####
##### 7.2 Gibbs 적용 #####
##### Ex 7.1
data <- read.table('Data/stress.txt', header = F)
stress <- data[,1]
SES <- data[,2]
mental <- data[,3]
ftable(stress, SES, mental)

y <- mental
n <- length(y)
X <- cbind(rep(1, n), stress, SES, stress*SES)
p <- ncol(X)
K <- 5

# prior
gamma0 <- rep(0, K) ; sigma0 <- rep(10, K)

# initial
beta <- rep(0, p) # beta zero
freq.y <- data.frame(table(y))$Freq 
prob <- freq.y / sum(freq.y)
gamma <- rep(0, K)
gamma[K] <- Inf
for(i in 2:(K-1)) {
  gamma[i] <- qnorm(sum(prob[1:i]), 0, 10)
  if (gamma[i] <= gamma[i-1]) gamma[i] <- gamma[i-1] + 0.1 
}
z <- rep(0, n)

# start MCMC
nsim <- 10000 ; nburn <- 5000 ; nthin <- 1
beta.Samples <- NULL
gamma.Samples <- NULL
for(iter in 1:(nsim*nthin+nburn)) {
  # generate z
  for (i in 1:n) {
    if (y[i] == 1) z[i] <- rtruncnorm(1, -Inf, gamma[1], X[i,]%*%beta, 1)
    for (k in 2:K) {
      if (y[i] == k) z[i] <- rtruncnorm(1, gamma[k-1], gamma[k], X[i,]%*%beta, 1)
    }
  }
  # generate beta
  var.beta <- n/(n+1) * solve(t(X)%*%X) 
  mean.beta <- var.beta %*% (t(X)%*%z)
  beta <- as.numeric(rmvnorm(1, mean.beta, var.beta))
  # generate gamma
  for (k in 2:(K-1)) {
    ak <- max(max(z[y==k]), gamma[k-1])
    bk <- min(min(z[y==(k+1)]), gamma[k+1])
    gamma[k] <- rtruncnorm(1, ak, bk, gamma0[k], sigma0[k])
  }      
  
  if (iter > nburn & iter%%nthin == 0) {
    beta.Samples <- rbind(beta.Samples, beta)
    gamma.Samples <- rbind(gamma.Samples, gamma) 
  }
}

# 수렴진단 
# thin = 1일 경우 gamma의 자기상관 매우 높다
par(mfrow = c(3,2))
for (i in 2:4) {
  plot(gamma.Samples[,i], type = 'l', xlab = 'iteration', ylab = '', main = paste0('gamma',i))
  acf(gamma.Samples[,i], main = paste0('gamma',i))
}
par(mfrow = c(4,2))
for (i in 1:p) {
  plot(beta.Samples[,i], type = 'l', xlab = 'iteration', ylab = '', main = paste0('beta',i-1))
  acf(beta.Samples[,i], main = paste0('beta',i-1))
}


##### 7.3 JAGS 활용 #####
modelString <- 'model
{
  for (i in 1:n) {
    y[i] ~ dcat(prob[i,1:K])
    mu[i] <- inprod(X[i,], beta[])
    prob[i,1] <- phi(gamma1-mu[i])
    prob[i,2] <- phi(gamma[2]-mu[i]) - phi(gamma1-mu[i])
    for (k in 3:(K-1)) {
      prob[i,k] <- phi(gamma[k]-mu[i]) - phi(gamma[k-1]-mu[i])
    }
    prob[i,K] <- 1 - phi(gamma[K-1] - mu[i])
  } 
  
  beta[1:p] ~ dmnorm(mu.beta, tau.beta)
  gamma1 <- 0
  mu.gamma <- 0
  tau.gamma <- 0.01
  for (k in 1:(K-1)) {
    gamma0[k] ~ dnorm(mu.gamma, tau.gamma)T(0,)
  }
  gamma[1:(K-1)] <- sort(gamma0)
}'
writeLines(modelString, 'ex7.1-ordered-probit.txt')
  
XtX <- t(X) %*% X
mu.beta <- as.numeric(solve(XtX) %*% t(X) %*% y) # rep(0, p)
tau.beta <- 1/n * XtX # diag(0.01, p)
dataList <- list(n = n, K = K, p = p, X = X, y = y, 
                 mu.beta = mu.beta, tau.beta = tau.beta)  

beta.init <- as.numeric(solve(XtX) %*% t(X) %*% y)
gamma.init <- (0:(K-2)) + 0.5 # gamma[1] = 0.5는 사용안함
initsList <- list(beta = beta.init, gamma0 = gamma.init)

ord.probit.Model <- jags.model('ex7.1-ordered-probit.txt', data = dataList, inits = initsList,
                               n.chains = 3, n.adapt = 1000)
update(ord.probit.Model, n.iter = 10000)
coda.Samples <- coda.samples(ord.probit.Model, variable.names = c('gamma[2:4]', 'beta'), thin = 1, n.iter = 10000)

# 채택확률 확인
1 - rejectionRate(coda.Samples)

# 수렴진단
# 경로그림과 자기상관
ESS <- effectiveSize(coda.Samples)
par(mfrow = c(4,2))
for (i in 1:4) {
  traceplot(coda.Samples[,i], main = variable.names(coda.Samples[[1]])[i])
  acf(coda.Samples[,i][[1]], main = paste0('ESS = ', round(ESS,2)[i]))
}
par(mfrow = c(3,2))
for (i in 5:7) {
  traceplot(coda.Samples[,i], main = variable.names(coda.Samples[[1]])[i])
  acf(coda.Samples[,i][[1]], main = paste0('ESS = ', round(ESS,2)[i]))
}

# 겔만 통계량
gelman.diag(coda.Samples)

# 사후밀도함수
mycol <- c('gray50', 'steelblue', 'orange', 'pink')
par(mfrow = c(1,2))
plot(density(as.matrix(coda.Samples[,1])), main = 'beta')
for (i in 2:4) lines(density(as.matrix(coda.Samples[,i])), main = 'beta', col = mycol[i])
plot(density(as.matrix(coda.Samples[,5])), main = 'gamma')
for (i in 5:7) lines(density(as.matrix(coda.Samples[,i])), main = 'beta', col = mycol[i-3])
