##### CH2. MCMC #####

##### 2.1 Sample-based Inference #####
##### Ex 2.1 
a <- 1 ; b <- 1
n <- 40 ; x <- 15

# True Posterior: theta ~ Beta(a+x, b+n-x)
aa <- a+x ; bb <- b+n-x
theta.grid <- seq(0, 1, by = 0.01)
post.true <- dbeta(theta.grid, aa, bb)
(true.mean <- aa / (aa+bb))
(true.sd <- sqrt((aa*bb) / ((aa+bb+1)*(aa+bb)**2)))

# Samples from the Posterior
theta.N500 <- rbeta(500, aa, bb)
theta.N5000 <- rbeta(5000, aa, bb)
theta.N50000 <- rbeta(50000, aa, bb)

# Plots
par(mfrow = c(2,2))
plot(theta.grid, post.true, type = 'l', col = 4, xlab = 'theta', ylab = 'true posterior')
text(0.7, 5, paste0('mean = ', round(true.mean,3), '\nsd = ', round(true.sd,3)))
hist(theta.N500, col = 'skyblue', xlab = 'N=500', ylab = '', xlim = c(0,1), main = '', nclass = 50)
text(0.7, 30, paste0('mean = ', round(mean(theta.N500),3), '\nsd = ', round(sd(theta.N500),3)))
hist(theta.N5000, col = 'skyblue', xlab = 'N=5000', ylab = '', xlim = c(0,1), main = '', nclass = 50)
text(0.7, 250, paste0('mean = ', round(mean(theta.N5000),3), '\nsd = ', round(sd(theta.N5000),3)))
hist(theta.N50000, col = 'skyblue', xlab = 'N=50000', ylab = '', xlim = c(0,1), main = '', nclass = 50)
text(0.7, 2500, paste0('mean = ', round(mean(theta.N50000),3), '\nsd = ', round(sd(theta.N50000),3)))

# 사후분포를 정확히 알 수 없는 경우가 많다
# 사후분포로부터 랜덤표본 생성(몬테칼로)이 가능하다면 충분히 많이 뽑아서 추론할 수 있다


##### 2.2 Metropolis-Hastings Algorithm #####
##### Ex 2.2 
theta <- 1:5
prob <- c(1, 3, 8, 5, 3) / 20

theta.curr <- 2 # 초기값

# alpha (채택확률) 계산 함수
prob.ratio <- function(theta1, theta2, theta, prob) {
  ind1 <- which(theta == theta1)
  ind2 <- which(theta == theta2)
  return(prob[ind1]/prob[ind2])
}

N <- 50000
theta.Samples <- c()
theta.Samples[1] <- theta.curr

# Simulation
for (iter in 2:N) {
  # 0.5의 확률로 왼쪽 오른쪽 중에 하나로 이동
  # 범위를 벗어나면 다시 현재 값
  theta.prop <- ifelse(runif(1) < 0.5, theta.curr +1, theta.curr -1) 
  if (theta.prop < 1 || theta.prop > 5) theta.prop <- theta.curr
  theta.prop <- round(theta.prop, 0)
  # 채택확률 계산
  alpha.star <- prob.ratio(theta.prop, theta.curr, theta, prob)
  alpha <- min(1, alpha.star)
  # Update
  theta.next <- ifelse(runif(1) < alpha, theta.prop, theta.curr)
  theta.Samples[iter] <- theta.next
  theta.curr <- theta.next
}

# Fig 2.4 초기스텝 100개의 경로그림
Ntrace <- 100
par(mfrow = c(1,1))
plot(theta.Samples[1:Ntrace], type = 'l', xlab = 'iteration', ylab = 'theta')
points(1:Ntrace, theta.Samples[1:Ntrace], pch = 19, col = 'steelblue')
# 초기 표본은 버리자 (수렴 전)

# Fig 2.5 실제 분포와 표본분포
par(mfrow = c(1,2))
barplot(prob, names.arg = theta, xlab = 'theta', ylab = 'prob', col = 'skyblue', sub = 'true probability')
prob.samp <- table(theta.Samples[-(1:Ntrace)]) / sum(theta.Samples[-(1:Ntrace)])
barplot(prob.samp, xlab = 'theta', ylab = 'prob', col = 'skyblue', sub = 'relative frequency of samples')


theta.Samples.2 <- c()
theta.Samples.2[1] <- theta.curr
for (iter in 2:N) {
  theta.prop <- sample(theta, 1) # 후보값을 theta space 아무데서나 뽑히도록 해보자 
  if (theta.prop < 1 || theta.prop > 5) theta.prop <- theta.curr
  theta.prop <- round(theta.prop, 0)
  # 채택확률 계산
  alpha.star <- prob.ratio(theta.prop, theta.curr, theta, prob)
  alpha <- min(1, alpha.star)
  # Update
  theta.next <- ifelse(runif(1) < alpha, theta.prop, theta.curr)
  theta.Samples.2[iter] <- theta.next
  theta.curr <- theta.next
}
# 실제분포와 비교 (거의 동일하다)
par(mfrow = c(1,2))
barplot(prob, names.arg = theta, xlab = 'theta', ylab = 'prob', col = 'skyblue', sub = 'true probability')
prob.samp <- table(theta.Samples.2[-(1:Ntrace)]) / sum(theta.Samples.2[-(1:Ntrace)])
barplot(prob.samp, xlab = 'theta', ylab = 'prob', col = 'skyblue', sub = 'relative frequency of samples')


# Independent Chain이 아닌 경우
prob.prop <- c(0.1, 0.1, 0.2, 0.3, 0.3) 
theta.Samples.3 <- c()
theta.Samples.3[1] <- theta.curr
for (iter in 2:N) {
  theta.prop <- sample(theta, 1, prob = prob.prop) # 후보값이 뽑히는 확률을 다르게 ㅎ
  if (theta.prop < 1 || theta.prop > 5) theta.prop <- theta.curr
  theta.prop <- round(theta.prop, 0)
  # 채택확률 계산
  alpha.star <- prob.ratio(theta.prop, theta.curr, theta, prob)
  alpha <- min(1, alpha.star)
  # Update
  theta.next <- ifelse(runif(1) < alpha, theta.prop, theta.curr)
  theta.Samples.3[iter] <- theta.next
  theta.curr <- theta.next
}
# 실제분포와 비교 (잘못된 proposal density 사용)
par(mfrow = c(1,2))
barplot(prob, names.arg = theta, xlab = 'theta', ylab = 'prob', col = 'skyblue', sub = 'true probability')
prob.samp <- table(theta.Samples.3[-(1:Ntrace)]) / sum(theta.Samples.3[-(1:Ntrace)])
barplot(prob.samp, xlab = 'theta', ylab = 'prob', col = 'skyblue', sub = 'relative frequency of samples')

# 후보 선택확률의 차이를 보정하는 방법(q function)
theta.Samples.4 <- c()
theta.Samples.4[1] <- theta.curr
for (iter in 2:N) {
  theta.prop <- sample(theta, 1, prob = prob.prop) 
  if (theta.prop < 1 || theta.prop > 5) theta.prop <- theta.curr
  theta.prop <- round(theta.prop, 0)
  # 채택확률 계산 
  alpha.star <- prob.ratio(theta.prop, theta.curr, theta, prob) /    ### 보정
                prob.ratio(theta.prop, theta.curr, theta, prob.prop) 
  alpha <- min(1, alpha.star)
  # Update
  theta.next <- ifelse(runif(1) < alpha, theta.prop, theta.curr)
  theta.Samples.4[iter] <- theta.next
  theta.curr <- theta.next
}
# 실제분포와 비교 (잘못된 proposal density 사용)
par(mfrow = c(1,2))
barplot(prob, names.arg = theta, xlab = 'theta', ylab = 'prob', col = 'skyblue', sub = 'true probability')
prob.samp <- table(theta.Samples.4[-(1:Ntrace)]) / sum(theta.Samples.4[-(1:Ntrace)])
barplot(prob.samp, xlab = 'theta', ylab = 'prob', col = 'skyblue', sub = 'relative frequency of samples')


##### 2.3 Gibbs Sampling Algorithm #####
##### Ex 2.3
nsim <- 3000 ; nburn <- 500 

# true parameters & data
mu0 <- 10 ; sigsq0 <- 25 ; a <- 0.5 ; b <- 1
x <- c(10, 13, 15, 11, 9, 18, 20, 17, 23, 21)
n <- length(x)

# initial value of mu, sigsq 
xbar <- mean(x) ; ssq <- var(x)
mu <- xbar ; sigsq <- ssq 

# Gibbs
theta.Samples <- matrix(nrow = nsim, ncol = 2)
for (iter in 1:nsim) {
  # generate mu
  condpost.sigsq <- 1 / (1/sigsq0 + n/sigsq)
  condpost.mu <- condpost.sigsq * (xbar/(sigsq/n) + mu0/sigsq0)
  mu <- rnorm(1, condpost.mu, sqrt(condpost.sigsq))
  
  # generate sigsq
  condpost.a <- a + n/2
  condpost.b <- b + 1/2*((n-1)*ssq + n*(xbar-mu)**2)
  sigsq <- 1 / rgamma(1, condpost.a, condpost.b)
  
  theta.Samples[iter,] <- c(mu, sigsq)
}

# Fig 2.7 준비단계 사후표본들의 경로그림
par(mfrow = c(1,3))
plot(theta.Samples[1:5,], type = 'n', xlab = 'mu', ylab = 'sigsq')
lines(theta.Samples[1:5,], lty = 2)
for(i in 1:5) text(theta.Samples[i,1], theta.Samples[i,2], i)
plot(theta.Samples[1:15,], type = 'n', xlab = 'mu', ylab = 'sigsq')
lines(theta.Samples[1:15,], lty = 2)
for(i in 1:15) text(theta.Samples[i,1], theta.Samples[i,2], i)
plot(theta.Samples[1:100,], type = 'n', xlab = 'mu', ylab = 'sigsq')
lines(theta.Samples[1:100,], lty = 2)
for(i in 1:100) text(theta.Samples[i,1], theta.Samples[i,2], i)

# Fig 2.8 초기 500개 제외 후 경로그림
par(mfrow = c(1,1))
plot(theta.Samples[-(1:nburn),], xlab = 'mu', ylab = 'sigsq', pch = 19)

# Fig 2.9 주변 사후밀도함수 추정
par(mfrow = c(1,2))
plot(density(theta.Samples[-(1:nburn),1]), xlab = 'mu', main = '')
plot(density(theta.Samples[-(1:nburn),2]), xlab = 'sigsq', main = '')

