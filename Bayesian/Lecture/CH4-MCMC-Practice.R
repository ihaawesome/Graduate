setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')

##### CH4. MCMC 실습 #####

##### 4.1 Gibbs #####

##### 4.1.1 선형회귀모형 분석
##### Ex 4.1
dat <- read.csv('Data/immigrants.csv')
y <- dat$wage
n <- length(y)
X <- cbind(rep(1, n), dat$sp, dat$lit)
p <- ncol(X)

a <- 1 ; b <- 1
XtX <- t(X) %*% X
XtX.inv <- solve(XtX)
Xty <- t(X) %*% y

# Least Squares Estimates
beta.lse <- as.vector(XtX.inv %*% t(X) %*% y)
beta.hat <- beta.lse
sigsq.hat <- sum((y - X%*%beta.lse)**2 / (n-p)) # MSE

# prior 
beta0 <- beta.lse
sig0 <- diag(diag(XtX.inv)*sigsq.hat)*100 # 무정보를 위해 분산 크게
sig0.inv <- solve(sig0)
# 원래 prior 지정은 최대한 데이터와 독립적으로 하면 좋음
# 데이터에서 추정한 값 지정할 경우 중복사용


# Random Draw from Multivariate Normal 
rmvnorm <- function(n, mu, sig) {
  p <- length(mu)
  R <- chol(sig)
  z <- matrix(rnorm(n*p), n, p)
  tt <- z %*% R + matrix(mu, n, p, byrow = T)
  return(tt)
}
# or
rmvnorm <- mvtnorm::rmvnorm


# initial values
beta.init <- beta.lse
sigsq.init <- sigsq.hat

# Gibbs 
N <- 10000 ; nburn <- 1000
sigsq.Samples <- rep(0, N)
beta.Samples <- matrix(0, N, p)

beta <- beta.init
sigsq <- sigsq.init 

# start sampling
for (iter in 1:(N+nburn)) {
  # generate beta
  sig.beta <- solve(sig0.inv + XtX/sigsq)
  mu.beta <- sig.beta %*% (sig0.inv%*%beta0 + 1/sigsq*Xty)
  beta <- as.vector(rmvnorm(1, mu.beta, sig.beta))
  # generate sigsq
  MSE <- sum((y - X%*%beta)**2 / (n-p))
  sigsq <- 1/ rgamma(1, n/2+a, 1/2*MSE+b)
  
  if (iter > nburn) {
  beta.Samples[iter-nburn,] <- beta
  sigsq.Samples[iter-nburn] <- sigsq
  }
}

# 95% HPD intervals
(ci.beta <- apply(beta.Samples, 2, quantile, probs = c(0.025, 0.975)))
(ci.sigsq <- quantile(sigsq.Samples, probs = c(0.025, 0.975)))

mycol <- c('gray50', 'steelblue', 'pink')
# Fig 4.1 경로그림과 자기상관
par(mfrow = c(2,2))
for (i in 1:p) {
  plot(beta.Samples[,i], type = 'l', xlab = paste0('beta',i), ylab = '', col = mycol[1])
}
plot(sigsq.Samples, type = 'l', xlab = 'sigsq', ylab = '', col = mycol[1])

# Fig 4.2 사후밀도함수와 95% HPD 신뢰구간
par(mfrow = c(2,2))
for (i in 1:p) {
  plot(density(beta.Samples[,i]), type = 'l', col = mycol[1], 
       xlab = paste0('beta',i), ylab = 'posterior', main = '')
  abline(v = ci.beta[,i], col = mycol[2], lty = 2)
}
plot(density(sigsq.Samples), type = 'l', col = mycol[2],
     xlab = 'sigsq', ylab = 'posterior', main = '')
abline(v = ci.sigsq, col = mycol[2], lty = 2)


##### 4.1.2 제한된 다변량 정규분포로부터의 표본 추출
library(truncnorm)
##### Ex 4.2 
k <- 5
mu <- c(0, 1, 2, 3, 5)
Sig <- matrix(0.7, k, k) + diag(k)*0.3
A <- solve(Sig)

theta.init <- 0:4
theta.Samples <- matrix(0, N, k)

theta <- theta.init
for (iter in 1:(N+nburn)) {
  for (i in 1:k) {
    A.else <- A[i,-i]
    mu.else <- (theta-mu)[-i]
    condmean <- mu[i] - 1/A[i,i]
    condsd <- 1 / sqrt(A[i,i])
    a <- as.double(ifelse(i == 1, -Inf, theta[i-1]))
    b <- as.double(ifelse(i == k, Inf, theta[i+1]))
    theta[i] <- rtruncnorm(1, a, b, condmean, condsd)
  }
  if (iter > nburn) theta.Samples[iter-nburn,] <- theta
}

# Fig 4.3 사후표본의 산점도 (제한조건 확인)
par(mfrow = c(2,2))
for (i in 1:4) {
  plot(theta.Samples[,c(i,i+1)], col = mycol[2],
       xlab = paste0('theta',i), ylab = paste0('theta',i+1))
  lines(theta.Samples[,i], theta.Samples[,i], type = 'l')
}
# Fig 4.4 사후밀도함수
for (i in 1:4) plot(density(theta.Samples[,i]), xlab = paste0('theta',i), main = '')


##### 4.2 M-H #####




