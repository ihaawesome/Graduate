setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')
library(rjags)
library(dplyr)

##### 10장. MCMC를 이용한 베이지안 변수선택 #####

##### 10.1 Stochastic Search Variable Selection (SSVS) #####

##### 예 10.1 #####
# Simulation Data
n <- 200 ; k <- 8
x <- matrix(rnorm(n*k), n, k)

beta0.true <- 0
beta.true <- c()
for (m in 1:k) beta.true[m] <- 0.5**(m-1)
y <- beta0.true + x %*% beta.true + rnorm(n)
y <- as.vector(y)

# JAGS
modelString <- '
model {
  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], invsigsq)
    mu[i] <- beta0 + inprod(x[i,], beta[])
  }
  gamma0 ~ dbern(0.5)
  for (j in 1:k) { gamma[j] ~ dbern(0.5) }
  
  beta0 ~ dnorm(mu.beta0, tau.beta0)
  mu.beta0 <- 0
  tau.beta0 <- (1-gamma0)/0.0001 + gamma0/100
  for (j in 1:k) {
    beta[j] ~ dnorm(mu.beta[j], tau.beta[j])
    mu.beta[j] <- 0
    tau.beta[j] <- (1-gamma[j])/0.0001 + gamma[j]/100 # spike and slab
  }
  invsigsq ~ dgamma(0.01, 0.01)
}
'
write(modelString, file = 'model-SSVS.txt')


dataList <- list(n = n, k = k, y = y, x = x)
lm.out <- lm(y ~ x)
mu.beta0 <- coef(lm.out)[1]
mu.beta <- coef(lm.out)[-1]
initsList <- list(beta0 = mu.beta0, beta = mu.beta, gamma0 = 1, gamma = rep(1, k))

nChains <- 3 ; nIter <- 10000
jagsModel <- jags.model('model-SSVS.txt', data = dataList, inits = initsList,
                        n.chains = nChains, n.adapt = 3000)
update(jagsModel, n.iter = 10000)
codaSamples <- coda.samples(jagsModel, variable.names = c('gamma0','gamma','beta0','beta'),
                            n.chains = nChains, n.iter = nIter)

mcmcSamples <- as.matrix(codaSamples)
theta.hat <- colMeans(mcmcSamples)    # *** sample posterior mean
beta.hat <- theta.hat[1:k]
gamma.hat <- theta.hat[(k+2):(2*k+1)]

as.data.frame(rbind(beta.true, beta.hat, gamma.hat))


##### Correlations of variables
z <- rnorm(n)
for (i in 1:4) x[,i] <- x[,i] + 2*z 
# x1부터 x4까지 모두 같은 값을 더해줌으로써 상관관계를 만든다 

# 랜덤 시뮬레이션 데이터기 때문에 변수선택이 약간씩 달라질 수 있으므로
# 50번 반복하여 서로 다른 시뮬레이션 데이터에서 추정된 모수의 평균을 구한다

mysimul <- function(n = 200, k = 8) {
# simulation data 
  x <- matrix(rnorm(n*k), n, k)
  z <- rnorm(n)
  for (i in 1:4) x[,i] <- x[,i] + 2*z # correlation
  
  beta0.true <- 0
  beta.true <- c()
  for (m in 1:k) beta.true[m] <- 0.5**(m-1)
  y <- beta0.true + x %*% beta.true + rnorm(n)
  y <- as.vector(y)
  
# jags model
  dataList <- list(n = n, k = k, y = y, x = x)
  gammaInit <- rep(1, k)
  lm.out <- lm(y ~ x)
  mu.beta0 <- coef(lm.out)[1]
  mu.beta <- coef(lm.out)[-1]
  initsList <- list(beta0 = mu.beta0, beta = mu.beta, gamma0 = 1, gamma = gammaInit)
  
  jagsModel <- jags.model('model-SSVS.txt', data = dataList, inits = initsList,
                          n.chains = 3, n.adapt = 3000)
  update(jagsModel, n.iter = 10000)
  codaSamples <- coda.samples(jagsModel, variable.names = c('gamma0','gamma','beta0','beta'),
                              n.chains = 3, n.iter = 10000)
  return(codaSamples)
}

nSim <- 50
simulList <- list()
for(iter in 1:nSim) simulList[[iter]] <- mysimul()

simulMeans <- lapply(simulList, function(l) c(quantile(as.matrix(l), probs = c(0.025, 0.975)))
simulSummary <- NULL
for (i in 1:nSim) simulSummary <- rbind(simulSummary, simulMeans[[i]])



##### 10.5 실제 자료 분석 : Baseball Salary Data 
##### 예 10.3 #####

dat <- read.csv('Data/baseball.txt')
colnames(dat)[1] <- 'y'

y <- log(dat$y)
n <- nrow(dat)
x <- dat[,-1]
X <- cbind(rep(1,n), x)
k <- ncol(X)

summary(lm.out <- lm(y ~ ., data.frame(y, x)))
summary(lm.step <- MASS::stepAIC(lm.out, direction = 'both'))

mu.beta <- lm.out$coefficients[1:k]
var.beta <- diag(vcov(lm.out))[1:k]

## pseudo prior ###
pseudo.beta.mean <- mu.beta
pseudo.beta.var <- var.beta

### prior #####
prior.beta.mean <- rep(0, k)
prior.beta.var <- var.beta*100 

# JAGS
modelString <- '
model { 
  for (j in 1:k) { gbeta[j]<- gamma[j]*beta[j] }
  for (i in 1:n) {
    y[i] ~ dnorm( mu[i], invsigsq)
    mu[i] <- inprod(X[i, 1:k], gbeta[1:k])
  }

  for (j in 1:k) { gamma[j] ~ dbern(0.5) }

  for (j in 1:k) {
    beta[j]~ dnorm( mu.b[j], tau.b[j])
    mu.b[j] <- gamma[j]*prior.beta.mean[j] + (1-gamma[j])*pseudo.beta.mean[j]
    tau.b[j] <- gamma[j]/prior.beta.var[j] + (1-gamma[j])/pseudo.beta.var[j]
  }
  invsigsq ~ dgamma(0.01, 0.01)
}
'
write(modelString, 'model-GVS-baseball.txt')

dataList <- list(n = n, k = k, X = X, y = y, 
                 pseudo.beta.mean = pseudo.beta.mean, pseudo.beta.var = pseudo.beta.var,
                 prior.beta.mean = prior.beta.mean, prior.beta.var = prior.beta.var)
initsList <- list(beta = mu.beta, gamma = rep(1, k))

jagsModel <- jags.model('model-GVS-baseball.txt', data = dataList, inits = initsList,
                        n.chains = 3, n.adapt = 3000)
update(jagsModel, n.iter = 10000)
codaSamples <- coda.samples(jagsModel, variable.names = c('gamma','beta'), n.iter = 30000)

theta.Samples <- as.matrix(codaSamples)
beta.Samples <- theta.Samples[,1:k]
gamma.Samples <- theta.Samples[,(k+1):(k+k)]

# 변수들의 동시선택:
# 벡터변수 gamma의 가능한 조합에 대한 빈도수와 상대도수를 계산한 후 
# 이를 내림차순으로 정렬하여 가장 빈도수가 높은 조합을 gamma의 추정치로 선택한다.
library(data.table)
m <- gamma.Samples
mm <- as.data.table(m)[, .N, by = eval(paste0('gamma[', seq_len(ncol(m)), ']'))]
colnames(mm) <- c(paste0('g', 0:16), 'N')
mm.order <- order(mm$N, decreasing = T)
mm$N <- round(mm$N/(3*30000), 4)
gamma.hat <- as.numeric(mm[which.max(mm$N)])
gamma.hat <- gamma.hat[1:k]
gamma.hat

df <- as.data.frame(gamma.Samples)
colnames(df) <- paste0('gamma',0:16)
df.count <- df %>% group_by_all %>% count %>% arrange(desc(n)) %>% ungroup
gamma.hat <- as.numeric(df.count[1,1:17])
gamma.hat

# as.data.frame(df.count)[1:10,]
# 
# gamma.Samples.collapsed <- apply(gamma.Samples, 1, function(x) paste(x, collapse = ' '))
# gamm.hat.collapsed <- paste(gamma.hat, collapse = ' ')
# id.selected <- which(gamma.Samples.collapsed == gamma.hat.collapsed)


##### 사후평균으로 변수 선택하는 방법
beta.hat.mean <- colMeans(beta.Samples)
gamma.hat.mean <- colMeans(gamma.Samples)




