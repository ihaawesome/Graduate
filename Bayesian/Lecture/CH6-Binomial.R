setwd('C:/Users/HK/Desktop/GitHub/Graduate/Bayesian')
library(rjags)

##### Ex 6.2 딱정벌레 #####
nn <- c(4, 4, 5, 3, 5, 4, 4, 4)
yy <- c(2, 1, 2, 3, 3, 3, 4, 4)
xx <- c(1.69, 1.72, 1.75, 1.78, 1.81, 1.83, 1.86, 1.88)


##### 6.2 로지스틱 모형 #####
##### 6.2.1 로지스틱 모형 #####

# Logistic M-H
x <- xx ; y <- yy ; n <- nn ; K <- length(y)
x <- x - mean(x)
X <- cbind(rep(1, K), x)
p <- ncol(X)

# log-posterior function
log.post <- function(beta, K, n, x, y, beta0, Sig0.inv) {
  Xbeta <- X %*% beta
  p <- 1 / (1+exp(-Xbeta))
  logpost <- sum(y*log(p) + (n-y)*log(1-p)) - 0.5*t(beta-beta0) %*% Sig0.inv %*% (beta-beta0)
}

#####
# Random Walk Metropolis-Hastings
glm.out <- glm(cbind(y,n-y) ~ x, family = binomial(link = "logit"))
beta.MLE <- as.vector(glm.out$coefficient)
beta <- beta.MLE
beta0 <- beta.MLE
Sig0.inv <- diag(0, p) # non-informative prior
XtX.inv <- solve(t(X)%*%X)

nsim <- 10000 ; nthin <- 1 ; nwarm <- 1000
beta.save <- matrix(0,nsim,p)
#####
# Independent Metropolis-Hastings
deltasq <- 9; nAccept <- 0
for (isim in 1:(nsim*nthin+nwarm)) {
  beta.star <- as.vector(rmvnorm(beta, deltasq*XtX.inv))
  log.alpha <- ( logpost(beta.star, K, n, X, y, beta0, Sig0.inv) - logpost(beta, K, n, X, y, beta0, Sig0.inv) ) +
               ( -0.5/deltasq*(t(beta-beta.MLE) %*% XtX %*% (beta-beta.MLE)) + 
                  0.5/deltasq*(t(beta.star-beta.MLE) %*% XtX %*% (beta.star-beta.MLE)) )
  u <- runif(1)
  if (log(u) < log.alpha) { beta <- beta.star ; nAccept <- nAccept +1 }
  if (isim > nwarm & isim%%nthin == 0) beta.save[(isim-nwarm)/nthin,] <- beta
}














