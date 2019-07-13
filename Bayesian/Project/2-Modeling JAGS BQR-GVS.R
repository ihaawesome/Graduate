library(tidyverse)
library(readxl)
library(glue)
library(quantreg)
library(rjags)
theme_set(theme_light() + theme(panel.grid = element_blank()))

# Load Data
MYDATA <- read.csv('DATA2018-FINAL-FOOD.csv', stringsAsFactors = F)

# 설명변수 Label
myvar <- c('주말매출비율','11~14시매출비율','14~17시매출비율','17~21시매출비율','21~24시매출비율',
           '여성매출비율','10대매출비율','20대미출비율','30대매출비율','점포수',
           '총상주인구','총유동인구','월평균소득','총소비')
names(myvar) <- colnames(MYDATA)[-(1:9)]

# JAGS model
modelString <- "
model { 
  for (j in 1:K) { gbeta[j] <- gamma[j]*beta[j] }

  for (i in 1:length(y)) {
    fd[i] <- p*(1-p)*exp(- check[i] )
    check[i] <- e[i]*(p-(1-step(e[i])))
    e[i] <- y[i] - inprod(x[i,1:K], gbeta[1:K]) 
  
    z[i] ~ dbern(pi[i]) 
    pi[i] <- fd[i]/10000
  }
  
  for (j in 1:K) { gamma[j] ~ dbern(0.5) }
  
  for(j in 1:K) { 
    beta[j] ~ dmnorm( mu[j], tau[j] )
    mu[j] <- (1- gamma[j]) * pseudo.mean.beta[j]
    tau[j] <- gamma[j]/100 + (1-gamma[j])/pseudo.var.beta[j]
  }
}
"
writeLines(modelString, "model-BQR-GVS.txt")


# 데이터 정의
# ------------------------------------------------------------------------------------------------
# 1000개만 sampling
set.seed(100)
sid <- sample(nrow(MYDATA), 1000)
MYDATA <- MYDATA[sid,]

x <- as.matrix(MYDATA[,-c(1:9)])
for (i in 1:ncol(x)) x[,i] <- (x[,i]-mean(x[,i]))/sd(x[,i])
x <- cbind(INTERCEPT = rep(1, nrow(x)), x)
x <- data.matrix(x)
head(x)
K <- ncol(x)
y <- MYDATA$SALES_MONTH_AMT / MYDATA$N_STORE / 1e06

# dataList, initsList (quantreg 추정치 사용) 
make.list <- function(p) {
  rq.out <- rq.fit.br(x, y, p, ci = T)
  pseudo.mean.beta <- rq.out$coefficients[1:K,1]
  pseudo.sd.beta <- (rq.out$coefficients[1:K,3]-rq.out$coefficients[1:K,2])/(2*1.96)
  pseudo.var.beta <- pseudo.sd.beta^2
  
  z <- rep(1, length(y))
  dataList = list(p = p, K = K, y = y, x = x, z = z,
                  pseudo.mean.beta = pseudo.mean.beta,
                  pseudo.var.beta = pseudo.var.beta)
  gammaInit <- rep(0, K)
  initsList <- list(beta = pseudo.mean.beta, gamma = gammaInit)
  
  return(list(dataList = dataList, initsList = initsList))
}


# JAGS 실행
# ------------------------------------------------------------------------------------------------
nAdapt <- 10000 ; nUpdate <- 10000 ; nIter <- 20000 ; nChains <- 3
p.set <- c(0.1, 0.25, 0.5, 0.75, 0.9)

# 변수 동시선택 계산
select.all <- function(gamma.Samples) {
  gamma.df <- as.data.frame(gamma.Samples)
  freq <- gamma.df %>% group_by_all() %>% summarise(FREQ = n()) %>% arrange(desc(FREQ))
  colnames(freq) <- c(colnames(x), 'FREQ')
  freq$FREQ <- freq$FREQ / nrow(gamma.df)
  top <- unlist(freq[1,])
  return(top)
}

# p에 따라 각각 샘플링 반복   
iterate.GVS <- function(p) {
  myList <- make.list(p)
  dataList <- myList$dataList
  initsList <- myList$initsList
  
  jagsModel <- jags.model(file = "model-BQR-GVS.txt", data = dataList, inits = initsList,
                          n.chains = nChains, n.adapt = nAdapt)
  update(jagsModel, n.iter = nUpdate)
  codaSamples <- coda.samples(jagsModel, variable.names = c("beta","gamma"), 
                              thin = 10, n.iter = nIter)
  return(codaSamples)
}

# 저장 공간
beta.hat <- matrix(0, length(p.set), K)
beta.cl <- matrix(0, length(p.set), K)
beta.cu <- matrix(0, length(p.set), K)

gamma.hat.mean <- matrix(0, length(p.set), K)
gamma.hat.all <- matrix(0, length(p.set), K+1)

codaSamples.set <- list()
mcmcSamples.set <- list()
beta.Samples.set <- list()
gamma.Samples.set <- list()

# ------------------------------------------------------------------------------ Iteratation ***
for (i in 1:length(p.set)) {
  p <- p.set[i]
  
  codaSamples.set[[i]] <- iterate.GVS(p)
  mcmcSamples.set[[i]] <- as.matrix(codaSamples.set[[i]])
  beta.Samples.set[[i]] <- mcmcSamples.set[[i]][,1:K]
  gamma.Samples.set[[i]] <- mcmcSamples.set[[i]][,-(1:K)]
  
  beta.hat[i,] <- apply(beta.Samples.set[[i]], 2, quantile, 0.5)
  beta.cl[i,] <- apply(beta.Samples.set[[i]], 2, quantile, 0.025)
  beta.cu[i,] <- apply(beta.Samples.set[[i]], 2, quantile, 0.975)
  
  gamma.hat.all[i,] <- select.all(gamma.Samples.set[[i]])
  gamma.hat.mean[i,] <- colMeans(gamma.Samples.set[[i]])
}


# 결과 요약 
# ----------------------------------------------------------------------------------------------
beta.Result <- as.data.frame(beta.hat) %>% gather(var, hat) %>% 
  mutate(p = rep(p.set, K), var = rep(variable.names(x), each = length(p.set))) %>%
  mutate(var = factor(var, levels = variable.names(x)))
beta.Result <- beta.Result %>% cbind(
  as.data.frame(beta.cl) %>% gather(var, lower) %>% select(-var),
  as.data.frame(beta.cu) %>% gather(var, upper) %>% select(-var)
) %>% filter(var != 'INTERCEPT')


# 분위에 따른 회귀계수, 신뢰구간
resultPlot1 <- beta.Result %>% 
  mutate(var = as.character(var)) %>%
  mutate(var = factor(myvar[var], levels = myvar)) %>%
  ggplot() + geom_line(aes(p, hat), size = 1, color = 'dodger blue') + 
  geom_line(aes(p, lower), size = 1, color = 'lightblue') + 
  geom_line(aes(p, upper), size = 1, color = 'lightblue') + 
  geom_hline(aes(yintercept = 0), linetype = 2) + 
  facet_wrap(~var, ncol = 5, scales = 'free') +
  labs(y = 'beta')

# 수렴 진단
for (i in 1:length(p.set)) gelman.diag(codaSamples.set[[i]][,1:K])

# ACF plots
par(mfrow = c(3, 5))
for (i in 1:length(p.set)) {
  for (j in 2:K) acf(beta.Samples.set[[i]][,j], main = myvar[j-1])
} 

# trace plot
par(mfrow = c(4, 4))
for (i in 1:length(p.set)) {
  for (j in 1:K) { 
    plot(mcmcSamples.set[[i]][1:2000,j], main = colnames(x)[j], 
         type = 'l', col = 'gray30', ylab = '')
    lines(mcmcSamples.set[[i]][2001:4000,j], col = 'steelblue')
    lines(mcmcSamples.set[[i]][4001:6000,j], col = 'salmon')
  }
}

# density
par(mfrow = c(3, 5))
for (i in 1:length(p.set)) {
  for (j in 2:K) { 
    plot(density(beta.Samples.set[[i]][,j]), main = myvar[j-1], xlab = '') 
    abline(v = beta.cl[i,j], lty = 2, col = 'steelblue')
    abline(v = beta.cu[i,j], lty = 2, col = 'steelblue')
  }
}

# 최종 변수선택 & 사후 추정치
# -----------------------------------------------------------------------------------------------
# top-1 gamma에 해당하는 표본들만 추출함
N <- nrow(gamma.Samples.set[[1]])
beta.Final.set <- list()
mys <- c()
for (i in 1:length(p.set)) { 
  for (ni in 1:N) { mys[ni] <- all(gamma.Samples.set[[i]][ni,] == gamma.hat.all[i,1:K]) }
  beta.Final.set[[i]] <- beta.Samples.set[[i]][mys,]
}

# 사후 평균
postMean <- sapply(beta.Final.set, colMeans)
rownames(postMean) <- c('INTERCEPT', myvar)
postMean*t(gamma.hat.all[,1:K])

# 사후 표준편차
postSD <- sapply(beta.Final.set, function(p) apply(p, 2, sd)) 
rownames(postSD) <- c('INTERCEPT', myvar)
postSD*t(gamma.hat.all[,1:K])

# 사후 95% 신뢰구간
temp <- lapply(beta.Final.set, function(p) apply(p, 2, quantile, c(0.025, 0.975)))
postCI <- t(temp[[1]])
for (a in 2:5) postCI <- cbind(postCI, t(temp[[a]]))
postCI

# summary plot
postMean <- as.data.frame(postMean)
postMean <- postMean %>% mutate(var = rownames(postMean))
colnames(postMean) <- p.set
postMean <- postMean %>% gather(p, beta, -var)

postCI <- as.data.frame(postCI, make.names = F)
colnames(postCI) <- str_c('p', rep(p.set, each = 2), colnames(postCI), sep = '_')
postCI <- postCI %>% mutate(var = c('INTERCEPT', myvar))
postCI <- postCI %>% gather(pp, value, -var)
postCI <- postCI %>% separate(pp, into = c('out','p','quantile'), sep = '_') %>% select(-out)
postCI1 <- postCI %>% spread(key = 'quantile', 'value') %>% 
  mutate(var = factor(var, levels = c('INTERCEPT', myvar))) %>% arrange(p, var)
postCI1 <- postCI1 %>% left_join(postMean)
postCI1 <- postCI1 %>% mutate(var = factor(var, c('INTERCEPT', myvar)))
postCI1 <- postCI1 %>% mutate(gamma = as.vector(t(gamma.hat.all[,1:K])))

# 분위에 따른 계수와 신뢰구간 비교
resultPlot2 <- postCI1 %>%
  filter(var != 'INTERCEPT') %>%
  ggplot() + geom_col(aes(p, beta, fill = factor(gamma))) + 
  geom_errorbar(aes(x = p, ymin = `2.5%`, ymax = `97.5%`), size = 1) + 
  facet_wrap(~var, scales = 'free', ncol = 5) + 
  scale_fill_manual(values = c('0' = 'gray90', '1' = 'dodger blue')) +
  theme(legend.position = '')



