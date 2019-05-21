setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')
library(rjags)
library(MASS)

##### CH9. 영과잉 카운트 자료의 분석
##### 9.1 영과잉 포아송 회귀모형 (ZIP)

##### Ex 9.1 ZIP #####
# Poisson 데이터와 ZIP 데이터
omega <- 0.4
n <- 200
x <- rnorm(n)
U <- rbinom(n, 1, omega) # Si
mu <- exp(rep(1, n) + x)
zip.y <- U*0 + (1-U)*rpois(n, mu)
mean.zip.y <- (1-omega)*mu
pois.y <- rpois(n, mean.zip.y)
table(pois.y)
table(zip.y)

# 포아송 자료와 영과잉 포아송 자료 비교
par(mfrow = c(2, 2))
hist(pois.y, freq = FALSE, breaks = c(-0.5:max(zip.y)), ylim = c(0, 0.6),
     main = '', xlab = 'poisson', col = 4)
hist(zip.y, freq = FALSE, breaks = c(-0.5:(max(zip.y)+1)), ylim = c(0, 0.6),
     main = '', xlab = 'ZIP', col = 4)
plot(exp(1+x), pois.y, xlim = c(0, 40), ylim = c(0, 30), 
     main = '', xlab = 'exp(1+x)', col = 4, pch = 16)
plot(exp(1+x), zip.y, xlim = c(0, 40), ylim = c(0, 30), 
     main = '', xlab = 'exp(1+x)', col = 4, pch = 16)


# 포아송 로그선형 모형 적합
glm.summary <- glm(zip.y ~ x, family = poisson(link = 'log'))
beta.pois <- as.vector(glm.summary$coefficients)
zip.y.hat.pois <- exp(beta.pois[1] + beta.pois[2]*x)
SSE.pois <- sum( (zip.y - zip.y.hat.pois)^2 /n )

# 영과잉 포아송 자료의 로그선형 모형 적합
par(mfrow = c(1, 2))
hist(zip.y, freq = FALSE, breaks = c(-0.5:(max(zip.y)+1)), ylim = c(0, 0.6),
     main = '', xlab = 'ZIP data', col = 4)
hist(zip.y.hat.pois, freq = FALSE, breaks = c(-0.5:(max(zip.y)+1)), ylim = c(0, 0.6),
     main = '', xlab = 'log-linear fit', col = 4)


##### 9.2 JAGS를 이용한 ZIP 모형의 베이지안 추론
##### Ex 9.1 JAGS (이어서) #####
# 0-1 방법
modelString <- '
model
{
  for (i in 1:length(y)) {
    Ones[i] ~ dbern(pi[i])
    pi[i] <- ( omega[i]*equals(y[i], 0) + (1-omega[i])*fd[i] ) / c
    fd[i] <- exp( y[i]*log(lambda[i]) - lambda[i] - loggam(y[i]+1) )
    log(lambda[i]) <- inprod(X[i,], beta[])
    logit(omega[i]) <- inprod(Z[i,], gamma[])
  }
  beta[1:p] ~ dmnorm(mu.beta[], Tau.beta[,])
  gamma[1:q] ~ dmnorm(mu.gamma[], Tau.gamma[,])
}
'
writeLines(modelString, 'model-ZIP.txt')

# 포아송 제로 방법 ( Poisson(0|mu) = exp(-mu) )
  # Zeros[i] ~ dpois(mu.pois[i])
  # mu.pois[i] <- -log.lik[i] + c # log.d = c
  # log.lik[i] <- log( omega[i]*equals(y[i], 0) + (1-omega[i])*fd[i] )
  # fd[i] <- exp( y[i]*log(lambda[i]) - lambda[i] - loggam(y[i]+1) )
  # log(lambda[i]) <- inprod(X[i,], beta[])
  # logit(omega[i]) <- inprod(Z[i,], gamma[])

# 합성분포 방법
  # y[i] ~ dpois(mu.pois[i])
  # mu.pois[i] <- (1-S[i])*lambda[i]
  # log(lambda[i]) <- inprod(X[i,], beta[])
  # S[i] ~ dbern(omega[i])


##### Ex 9.3 ZIP #####
data <- read.csv('Data/car_sample_300.csv')
y <- data$ninjured
table(y)

par(mfrow = c(1, 1))
hist(y, freq= FALSE, breaks = c(-0.5:(max(y)+1)), main = '',  xlab = 'ninjured', col = 4)

X <- model.matrix( ~ night + weekend + intersection + acc1 + acc2, data)
Z <- model.matrix( ~ acc1 + acc2 + rule, data)
p <- ncol(X)
q <- ncol(Z)

# (1)
modelString <- "
model
{
  for (i in 1:length(y)) {
    y[i] ~ dpois(mu.pois[i])
    mu.pois[i] <- (1-S[i])*lambda[i] + 1e-10*S[i]
    log(lambda[i]) <- inprod (X[i,], beta[])
 
    S[i] ~ dbern( omega[i] )
    logit(omega[i]) <- inprod (Z[i,], gamma[])
  }
 
  for (i in 1:p) {
    beta[i] ~ dnorm(mu.beta[i], Tau.beta[i])
  }
  for (i in 1:q) {
    gamma[i] ~ dnorm( mu.gamma[i], Tau.gamma[i])
  }
}
"
writeLines(modelString, "model-ZIP-mixture.txt")

# (2) prior parameters
mu.beta <- rep(0, p)
Tau.beta <- rep(0.01, p)
mu.gamma <- rep(0, q)
Tau.gamma <- rep(0.01, q)

glm.out <- glm(y ~ X -1, family = 'poisson')
beta.pois <- as.vector(glm.out$coefficients)

dataList <- list(p = p, q = q, y = y, X = X, Z = Z, 
                 mu.beta = mu.beta, Tau.beta = Tau.beta,
                 mu.gamma = mu.gamma, Tau.gamma = Tau.gamma)
initsList <- list(beta = beta.pois, gamma = mu.gamma)

jagsModel.zip <- jags.model('model-ZIP-mixture.txt', data = dataList, inits = initsList,
                            n.chains = 3, n.adapt = 1000)
update(jagsModel.zip, n.iter = 3000)
codaSamples <- coda.samples(jagsModel.zip, variable.names = c('beta', 'gamma'),
                            thin = 1, n.iter = 10000)
# (3) results
gelman.diag(codaSamples)
summary(codaSamples)
# 
(dic.zip <- dic.samples(jagsModel.zip, n.iter = 30000))


##### 9.3 영과잉 음이항 모형 (ZINB)
##### Ex 9.3 ZINB (이어서) #####
modelString <- '
model
{
  for (i in 1:length(y)) {
    
  }
}
'
writeLines(modelString, 'model-ZINB.txt')
