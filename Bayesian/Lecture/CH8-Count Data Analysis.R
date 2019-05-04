setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')
library(rjags)
library(ggmcmc)
library(ggplot2)

##### CH8. 카운트 자료의 분석 #####

##### 8.1 포아송 로그선형 모형 #####
##### Ex 8.1
mers <- read.table('Data/mers.txt', header = T)
head(mers)

day <- as.character(mers$day)
diag <- mers$diag_new
y <- diag
n <- length(y)
n.day <- (1:n)
plot(n.day, y, xlab = 'day', ylab = 'n_diag', type = 'l')
points(n.day, y)

x1 <- c(rep(1, 19), rep(0, n-19))
x2 <- c((1:19), rep(0, n-19))
x3 <- c(rep(0, 19), rep(1, n-19))
x4 <- c(rep(0, 19), (1:(n-19)))
X <- cbind(x1, x2, x3, x4)
data <- data.frame(y, X)
p <- ncol(X)

# JAGS
modelString <- 'model
{
  for (i in 1:length(y)) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- inprod(X[i,], beta[])
  }
  for (i in 1:p) { beta[i] ~ dnorm(mu.beta[i], Tau.beta[i]) }
}'
writeLines(modelString, 'model-8.1-pois.txt')
  
# prior parameters
mu.beta <- rep(0, p)
Tau.beta <- rep(0.01, p)

dataList <- list(p = p, y = y, X = X, mu.beta = mu.beta, Tau.beta = Tau.beta)
initsList <- list(beta = mu.beta)

# sampling
poisModel <- jags.model('model-8.1-pois.txt', data = dataList, inits = initsList, n.chains = 3, n.adapt = 1000)
update(poisModel, n.iter = 3000)
codaSamples <- coda.samples(poisModel, variable.names = c('beta'), n.iter = 10000, thin = 1)

# 수렴진단과 추정치 분포
gelman.diag(codaSamples)
summary(codaSamples)

# 추정치와 관측치 비교
mcmcSamples <- as.matrix(codaSamples)
beta.hat <- apply(mcmcSamples, 2, mean)
lambda.hat <- exp(X %*% beta.hat)
data.frame(n.day, y, lambda.hat) %>%
  ggplot() + theme_bw() +
  geom_point(aes(n.day, y)) +
  geom_line(aes(n.day, lambda.hat), color = 'red') +
  geom_point(aes(n.day, lambda.hat), color = 'red', shape = 4) +
  labs(x = 'days', y = 'diag')

# DIC
(dic.pois <- dic.samples(poisModel, n.iter = 10000))

# DIC를 이용한 모형 비교를 위해 day 변수의 제곱항을 추가해 다시 계산해보자
X <- cbind(x1, x2, x3, x4, x*x4)
p <- ncol(X)
mu.beta <- rep(0, p)
Tau.beta <- rep(0.01, p)
dataList <- list(p = p, y = y, X = X, mu.beta = mu.beta, Tau.beta = Tau.beta)
initsList <- list(beta = mu.beta)
poisModel2 <- jags.model('model-8.1-pois.txt', data = dataList, inits = initsList, n.chains = 3, n.adapt = 1000)
(dic.pois2 <- dic.samples(poisModel2, 10000))
# 제곱항이 없는 모형과 DIC 값이 거의 차이가 없으므로 두 모형의 적합도에 차이가 없다고 본다
# 더 간단한 모형이 더 적합하다고 볼 수 있다


##### Ex 8.2
KL <- read.csv('Data/KL1.csv')
head(KL)
table(KL$team1)
table(KL$team2)

except.team <- c('부산','성남','아산')
KL <- KL %>% filter(!(team1 %in% except.team) & !(team2 %in% except.team))
name.team1 <- data.frame(num.team1 = 1:12, team1 = setdiff(names(table(KL$team1)), except.team))
name.team2 <- data.frame(num.team2 = 1:12, team2 = setdiff(names(table(KL$team2)), except.team))
KL <- KL %>% merge(name.team1, by = 'team1') %>% merge(name.team2, by = c('team2'))
head(KL)

# JAGS
modelString <- 'model 
{
  for (i in 1:n) {
    y1[i] ~ dpois(lambda1[i])
    y2[i] ~ dpois(lambda2[i])
    log(lambda1[i]) <- mu + home + a[ht[i]] + d[at[i]]
    log(lambda2[i]) <- mu        + a[at[i]] + d[ht[i]]
  }
  for (k in 2:K) {
    a[k] ~ dnorm(0, 0.0001)
    d[k] ~ dnorm(0, 0.0001)
  }
  a[1] <- sum(a[2:K])
  d[1] <- sum(d[2:K])

  mu ~ dnorm(0, 0.0001)
  home ~ dnorm(0, 0.0001)
}'
writeLines(modelString, 'model-8.2-pois.txt')

y1 <- KL$score1 ; y2 <- KL$score2
ht <- KL$num.team1 ; at <- KL$num.team2
n <- nrow(KL)
K <- length(unique(KL$team1))
initsList <- list(mu = 0, home = 0)
dataList <- list(n = n, K = K, at = at, ht = ht, y1 = y1, y2 = y2)

poisModel <- jags.model('model-8.2-pois.txt', data = dataList, inits = initsList, n.chains = 3, n.adapt = 1000)


##### Ex 8.3



