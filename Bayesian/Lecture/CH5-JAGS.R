setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')
library(rjags)

##### CH5. JAGS #####
##### Ex 5-1

# JAGS 모형 설정
modelString <- "model
{
  for (i in 1:n) {
    x[i] ~ dnorm(mu, invsigsq) # fixed node
  }
  
  mu ~ dnorm(mu0, invsigsq0) # random node
  invsigsq ~ dgamma(a, b)    # random node
  
  sigsq <- 1/invsigsq        # fixed node 
  mu0 <- 10                  # fixed node
  invsigsq0 <- 1/25          # fixed node
  a <- 0.5
  b <- 1
}"
writeLines(modelString, "model-ex5.1.txt")

# 리스트 형식으로 입력할 데이터를 저장
dataList <- list(n = 10, x = c(10, 13, 15, 11, 9, 18, 20, 17, 23, 21)) 

# 리스트 형식으로 랜덤 노드의 초기치를 설정
# prior 값 사용
initsList <- list(mu = 10, invsigsq = 1/25)  

# 랜덤 초기치 사용 
initsList <- function(x) {
  resampledX <- sample(x, replace = T)
  muInit <- sum(resampledX) / length(resampledX)
  invsigsqInit <- (1/sd(resampledX))^2*0.9999 + 0.01
  return(list(mu = muInit, invsigsq = invsigsqInit))  
}

# JAGS 실행
jagsModel <- jags.model(file = "model-ex5.1.txt", data = dataList, inits = initsList, 
                        n.chains = 3, n.adapt = 500)       # initial은 꼭 주진 않아도 됨
                      # 다중 체인의 수, optimal proposal density 등을 찾기 위한 초기튜닝 횟수

# 예비단계
update(jagsModel, n.iter = 500) # burn.in하고 저장은 하지 않는다
# 복잠한 모형에서는 burn-in 몇 천번 정도는 해야 됨

# MCMC - update 이후에 이어서 실행
codaSamples <- coda.samples(jagsModel, variable.names = c("mu", "sigsq"), n.iter = 5000)
                                     # 저장하고 싶은 variable

#####
# 수렴진단
# 경로그림과 자기상관계수
par(mfrow = c(2,2))
traceplot(codaSamples[,"mu"], main = "mu")
acf(codaSamples[,"mu"][[1]], main = "mu")
traceplot(codaSamples[,"sigsq"], main = "sigsq")
acf(codaSamples[,"sigsq"][[1]], main = "sigsq")

# Gelman 통계량
gelman <- gelman.diag(codaSamples)
gelman1 <- as.matrix(gelman$psrf)
if (max(gelman1) > 1.1) cat ("Warning : Gelman Shrink Factor > 1.1\n")
gelman2 <- gelman$mpsrf
if(gelman2 > 1.1) cat("Warning : Gelman Multivariate Shrink Factor > 1.1\n")

# ESS
# check MCMC Efficiency
nChains <- 3
mcmcSamples.combined <- NULL
for (ich in 1:nChains) mcmcSamples.combined <- rbind(mcmcSamples.combined, mcmc(codaSamples[[ich]]))
ESS <- effectiveSize(mcmcSamples.combined)
cat("Effective Sample Size = ", ESS) 

muSamples <- as.matrix(codaSamples[,"mu"])
sigsqSamples <-as.matrix(codaSamples[,"sigsq"])

# 사후추론
par(mfrow = c(1,2))
plot(density(muSamples), main = "", xlab = bquote(mu), ylab = "posterior density")
plot(density(sigsqSamples), main = "", xlab = bquote(sigma^2), ylab = "posterior density")

# 채택확률
acceptRate <- 1 - rejectionRate(codaSamples)
acceptRate


##### Ex 5-2
dat <- read.csv("Data/immigrants.csv")
y <- dat$wage
n <- length(y)
X <- cbind(rep(1,n), dat$sp, dat$lit)
p <- ncol(X)

a <- 1
b <- 1

XtX <- t(X) %*% X
XtX.inv <- solve(XtX)
Xty <- t(X) %*% y
beta.hat <- as.vector(XtX.inv %*% Xty)
sigsq.hat <- sum((y - X %*% beta.hat)^2)/(n-p)

beta0 <- beta.hat
Sig0 <- diag(diag(XtX.inv))*sigsq.hat*100
Sig0.inv <- solve(Sig0)

#####
# JAGS 모형
modelString <- "
model 
{
  for (i in 1:length(y)) {
    y[i] ~ dnorm(inprod(X[i,], beta[]), invsigsq)
  }

  beta[1:length(beta0)] ~ dmnorm(beta0[], Sig0.inv[,])
  invsigsq ~ dgamma(a, b)

  sigsq <- 1/invsigsq
}
"
writeLines(modelString, "model-ex5.2-reg.txt")

dataList <- list(X = X, y = y, a = a, b = b, beta0 = beta0, Sig0.inv = Sig0.inv)
initsList <- list(beta = beta.hat, invsigsq = 1/sigsq.hat)

nChains <- 3

jagsModel <- jags.model("model-ex5.2-reg.txt", data = dataList, inits = initsList,
                        n.chains = nChains, n.adapt = 500)

update(jagsModel, n.iter = 1000) 
codaSamples <- coda.samples(jagsModel, variable.names = c("beta", "sigsq"),
                            n.iter = 30000)
para.names <- variable.names(codaSamples[[1]])

#####
# 경로그림과 자기상관
par(mfrow = c(4,2))
for (i in 1:4) {
  traceplot(codaSamples[,i], main = "", ylab = para.names[i])
  acf(codaSamples[,i][[1]], main = para.names[i])
}
# savePlot("figure-5.7", type = c("bmp"), device = dev.cur())

MCMCSamples <- as.matrix(codaSamples)
(HPD <- round(apply(MCMCSamples, 2, quantile, probs = c(0.025, 0.975)), 4))

# 사후밀도함수와 95% HPD
par(mfrow = c(2,2))
for (i in 1:4) {
  plot(density(MCMCSamples[,i]), main = "", xlab = para.names[i], col = "gray50")
  abline(v = HPD[,i], col = "steelblue") 
}
# savePlot("figure-5.8", type = c("bmp"), device = dev.cur())
