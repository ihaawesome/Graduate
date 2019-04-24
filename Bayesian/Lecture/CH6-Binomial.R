setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')
library(rjags)
library(mvtnorm)
library(truncnorm)

##### CH6. 이항자료의 베이지안 분석 #####

##### 6.1 프로빗 모형 
##### 6.1.2 잠재변수를 도입한 Gibbs Sampling
rmvnorm <- function(mu, Sig) {
  R <- t(chol(Sig))
  R %*% (rnorm(length(mu))) + mu
}
rtruncnorm <- function(mu, sig, l, u) {
  unif <- runif(1)
  truncn <- qnorm( unif*pnorm((u-mu)/sig) + (1-unif)*pnorm((l-mu)/sig))*sig + mu
  return(truncn)
}

##### Ex 6.1
nn <- c(4, 4, 5, 3, 5, 4, 4, 4) # insects
yy <- c(2, 1, 2, 3, 3, 3, 4, 4) # killed
xx <- c(1.69, 1.72, 1.75, 1.78, 1.81, 1.83, 1.86, 1.88) # dosage

K <- length(nn)
x <- rep(xx, nn) # each = nn
SF <- as.vector(rbind(yy, nn-yy))
y <- rep(rep(c(1,0), K), SF)

x <- x - mean(x)
n <- length(x)
X <- cbind(rep(1, n), x)
p <- ncol(X)

beta0 <- rep(0, p)
Sig0.inv <- diag(0, p) # non-informative (Inf)

# initial values
beta <- beta0

# sampling
nsim <- 10000 ; nburn <- 100
beta.Samples <- matrix(0, nsim, p)

Sig.beta <- solve(t(X)%*%X + Sig0.inv)
y.star <- rep(0, n)

# MCMC
for (iter in 1:(nsim+nburn)) {
  # generate y.star
  for (i in 1:n) {
    if (y[i] == 1) y.star[i] <- rtruncnorm(t(X[i,])%*%beta, 1, 0, Inf)
    else y.star[i] <- rtruncnorm(t(X[i,])%*%beta, 1, -Inf, 0)
  }
  # generate beta
  mu.beta <- Sig.beta %*% ( t(X)%*%y.star + Sig0.inv%*%beta0 )
  beta <- rmvnorm(mu.beta, Sig.beta)
  if (iter > nburn) beta.Samples[iter-nburn,] <- beta
}

# 수렴진단
# 경로그림과 자기상관
par(mfrow = c(2,2))
plot(beta.Samples[,1], type = 'l', main = 'beta0', xlab = 'iteration', ylab = '')
plot(beta.Samples[,2], type = 'l', main = 'beta1', xlab = 'iteration', ylab = '')
acf(beta.Samples[,1], main = 'beta0')
acf(beta.Samples[,2], main = 'beta1')

# 사후밀도함수와 HPD
beta1.HPD <- quantile(beta.Samples[,2], c(0.025, 0.975))
par(mfrow = c(1,1))
plot(density(beta.Samples[,2]), main = 'beta1', xlab = '', ylab = 'posterior')
abline(v = beta1.HPD, col = 'orange', lty = 2)


##### 6.1.2 JAGS를 이용한 사후표본 생성
##### Ex 6.1 (이어서)
modelString1 <- '
model
{
  for (i in 1:K) {
  y[i] ~ dbin(pi[i], n[i])
  probit(pi[i]) <- beta0 + beta1*(x[i]-mean(x[]))
  }
  beta0 ~ dnorm(mu0, invsig0)
  beta1 ~ dnorm(mu1, invsig1) 
}'
writeLines(modelString1, 'Lecture/ex6.1-probit-ver1.txt')
  
# vector version  
modelString2 <- '
model
{
  for (i in 1:length(y)) {
    y[i] ~ dbin(pi[i], n[i])
    probit(pi[i]) <- inprod(X[i,], beta[])
  }
  beta[1:length(beta0)] ~ dmnorm(beta0[], Sig0.inv[,])
}'  
writeLines(modelString2, 'Lecture/ex6.1-probit-ver2.txt')

### ver 1
n <- nn ; x <- xx ; y <- yy
K <- length(y)
dataList <- list(n = n, y = y, x = x, K = K, mu0 = 0, mu1 = 0, invsig0 = 0.0001, invsig1 = 0.0001)
initsList <- list(beta0 = 0, beta1 = 0)

jags.Model <- jags.model('Lecture/ex6.1-probit-ver1.txt', data = dataList, inits = initsList, n.chains = 3, n.adapt = 500)
update(jags.Model, n.iter = 1000)
coda.Samples <- coda.samples(jags.Model, variable.names = c('beta0', 'beta1'), n.iter = 10000)

# 채택확률 확인
1 - rejectionRate(coda.Samples) # Gibbs

# 수렴진단
convdiag.plot <- function(coda.Samples, var) {
  ESS <- effectiveSize(coda.Samples[,var])
  traceplot(coda.Samples[,var], main = var)
  acf(as.matrix(coda.Samples[,var]), main = paste0('ESS = ', round(ESS,2)))
}

par(mfrow = c(2,2))
convdiag.plot(coda.Samples, 'beta0')
convdiag.plot(coda.Samples, 'beta1')

### ver 2
X <- cbind(rep(1, K), x-mean(x))
p <- ncol(X)
beta0 <- rep(0, p)
Sig0.inv <- diag(0.0001, p)
dataList <- list(n = n, X = X, y = y, beta0 = beta0, Sig0.inv = Sig0.inv)
initsList <- list(beta = beta0)

jags.Model <- jags.model('ex6.1-probit-ver2.txt', data = dataList, inits = initsList, n.chains = 3, n.adapt = 500)
update(jags.Model, n.iter = 1000)
coda.Samples <- coda.samples(jags.Model, variable.names = c('beta'), n.iter = 10000)

# 채택확률 확인 
1 - rejectionRate(coda.Samples) # MH

# 수렴진단
# 경로그림과 자기상관
convdiag.plot(coda.Samples, 'beta[1]')
convdiag.plot(coda.Samples, 'beta[2]')
# 겔만 통계량
gelman.diag(coda.Samples)

# 사후밀도함수와 HPD
mcmc.Samples.combined <- NULL
for (i in 1:3) mcmc.Samples.combined <- rbind(mcmc.Samples.combined, mcmc(coda.Samples[[i]]))
beta.HPD <- apply(mcmc.Samples.combined, 2, quantile, prob = c(0.025, 0.975))
par(mfrow = c(1,2))
plot(density(as.matrix(coda.Samples[,1])), main = 'beta0')
abline(v = beta.HPD[,1], col = 'orange', lty = 2)
plot(density(as.matrix(coda.Samples[,2])), main = 'beta1')
abline(v = beta.HPD[,2], col = 'orange', lty = 2)


##### 6.2 로지스틱 모형

##### 6.2.2 랜덤워크 메트로폴리스
x <- xx ; y <- yy ; n <- nn
K <- length(y) 
x <- x - mean(x)
X <- cbind(rep(1, K), x)
p <- ncol(X)

glmFit <- glm(cbind(y, n-y) ~ x, family = binomial)
beta.mle <- as.numeric(coef(glmFit)) 

beta0 <- beta.mle
Sig0.inv <- diag(0, p) # non-informative
XtX <- t(X) %*% X
XtX.inv <- solve(XtX)

dataList <- list(n = n, X = X, y = y, beta0 = beta0, Sig0.inv = Sig0.inv)

# log-posterior function
log.post.kernel <- function(beta, dataList) {
  n <- dataList$n
  X <- dataList$X
  y <- dataList$y
  K <- length(y)
  beta0 <- dataList$beta0
  Sig0.inv <- dataList$Sig0.inv
  
  Xbeta <- X %*% beta
  p <- 1 / (1 + exp(-Xbeta))
  log.post <- sum(y*log(p) + (n-y)*log(1-p)) - 0.5*t(beta-beta0)%*%Sig0.inv%*%(beta-beta0)
  return(log.post)
}

nsim <- 10000 ; nburn <- 1000 ; nthin <- 1
beta.Samples <- matrix(0, nsim, p)

# Metropolis 
deltasq <- 9 ; naccept <- 0
beta.curr <- beta0
for (iter in 1:(nsim*nthin+nburn)) {
  beta.prop <- as.numeric(rmvnorm(beta.curr, deltasq*XtX.inv))
  log.alpha <- log.post.kernel(beta.prop, dataList) - log.post.kernel(beta.curr, dataList)
  if (log(runif(1)) < log.alpha) { beta.next <- beta.prop ; naccept <- naccept +1 } 
  else { beta.next <- beta.curr }
  
  if (iter > nburn & iter%%nthin == 0) { beta.Samples[(iter-nburn)/nthin,] <- beta.next }
  beta.curr <- beta.next
}

# 수렴 진단
par(mfrow = c(2,2))
plot(beta.Samples[,1], type = 'l', xlab = 'iteration', ylab = '', main = 'beta0', col = 'gray50')
plot(beta.Samples[,2], type = 'l', xlab = 'iteration', ylab = '', main = 'beta1', col = 'gray50')
acf(beta.Samples[,1], main = 'beta0')
acf(beta.Samples[,2], main = 'beta1')

# 사후 밀도함수
par(mfrow = c(1,2))
plot(density(beta.Samples[,1]), main = 'beta0')
plot(density(beta.Samples[,2]), main = 'beta1')


##### 6.2.3 독립 메트로폴리스-해스팅스
nsim <- 10000 ; nburn <- 1000 ; nthin <- 1
deltasq <- 2.5 ; naccept <- 0
beta.Samples <- matrix(0, nsim, p)

beta.curr <- beta0
for (iter in 1:(nsim*nthin+nburn)) {
  beta.prop <- as.numeric(rmvnorm(beta.mle, deltasq*XtX.inv))
  log.alpha <- log.post.kernel(beta.prop, dataList) - log.post.kernel(beta.curr, dataList) +
               ( -0.5/deltasq*(t(beta.curr-beta.mle) %*% XtX %*% (beta.curr-beta.mle)) + 
                  0.5/deltasq*(t(beta.prop-beta.mle) %*% XtX %*% (beta.prop-beta.mle)) )
  
  if (log(runif(1)) < log.alpha) { beta.next <- beta.prop ; naccept <- naccept +1 } 
  else { beta.next <- beta.curr }
  
  if (iter > nburn & iter%%nthin == 0) { beta.Samples[(iter-nburn)/nthin,] <- beta.next }
  beta.curr <- beta.next
}

# 수렴 진단
par(mfrow = c(2,2))
plot(beta.Samples[,1], type = 'l', xlab = 'iteration', ylab = '', main = 'beta0', col = 'gray50')
plot(beta.Samples[,2], type = 'l', xlab = 'iteration', ylab = '', main = 'beta1', col = 'gray50')
acf(beta.Samples[,1], main = 'beta0')
acf(beta.Samples[,2], main = 'beta1')

# 사후 밀도함수
par(mfrow = c(1,2))
plot(density(beta.Samples[,1]), main = 'beta0')
plot(density(beta.Samples[,2]), main = 'beta1')

# 랜덤워크 체인에서는 nthin = 1으로 했을 때 완전히 수렴하지 않은 것처럼 보인다
# 랜덤워크는 선택확률 보정 부분이 약분된다

# 독립 체인에서는 nthin = 1으로 하더라도 충분히 수렴 한 것으로 보인다
# 독립 체인은 beta.mle로부터 beta.prop을 뽑을 확률과 beta.mle로부터 beta.curr을 뽑을 확률이 다르므로
# alpha를 계산할 때 보정 부분을 곱해주어야 한다


##### 6.2.3 JAGS 활용
modelString <- 'model
{
  for (i in 1:length(y)) {
    y[i] ~ dbin(pi[i], n[i])
    logit(pi[i]) <- beta0 + beta1*(x[i]-mean(x[]))
  }
  beta0 ~ dnorm(mu0, invsig0)
  beta1 ~ dnorm(mu1, invsig1)
}'
writeLines(modelString, 'ex6.1-logit.txt')

dataList <- list(n = n, y = y, x = x, mu0 = 0, mu1 = 0, invsig0 = 0.0001, invsig1 = 0.0001)
initsList <- list(beta0 = 0, beta1 = 0)
  
jags.Model <- jags.model('ex6.1-logit.txt', data = dataList, inits = initsList, n.chains = 3, n.adapt = 500)
update(jags.Model, n.iter = 1000)
coda.Samples <- coda.samples(jags.Model, variable.names = c('beta0', 'beta1'), n.iter = 10000)

1 - rejectionRate(coda.Samples) # Gibbs
gelman.diag(coda.Samples)

par(mfrow = c(2,2))
traceplot(coda.Samples[,1], xlab = 'iteration', ylab = '', main = 'beta0')
traceplot(coda.Samples[,2], xlab = 'iteration', ylab = '', main = 'beta1')
acf(coda.Samples[,1][[1]], main = 'beta0')
acf(coda.Samples[,2][[2]], main = 'beta1')

par(mfrow = c(1,2))
plot(density(as.matrix(coda.Samples[,1])), main = 'beta0')
plot(density(as.matrix(coda.Samples[,2])), main = 'beta1')


