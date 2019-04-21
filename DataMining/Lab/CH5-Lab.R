library(ISLR)
library(MASS)
library(class)

########## Chapter 5: Cross-Validation and the Bootstrap  ########## 

##### 5.3.1 Validation Set Approach #####

set.seed(1)
train <- sample(392, 196)

attach(Auto)
lm.fit <- lm(mpg ~ horsepower, data = Auto, subset = train)
mean((mpg - predict(lm.fit, Auto))[-train]^2) # MSE

lm.fit2 <- lm(mpg ~ poly(horsepower, 2), data = Auto, subset = train)
lm.fit3 <- lm(mpg ~ poly(horsepower, 3), data = Auto, subset = train)
mean((mpg - predict(lm.fit2, Auto))[-train]^2)
mean((mpg - predict(lm.fit3, Auto))[-train]^2)

# change training set
set.seed(2)
train <- sample(392, 196)
lm.fit <- lm(mpg ~ horsepower, data = Auto, subset = train)
lm.fit2 <- lm(mpg ~ poly(horsepower, 2), data = Auto, subset = train)
lm.fit3 <- lm(mpg ~ poly(horsepower, 3), data = Auto, subset = train)
mean((mpg - predict(lm.fit, Auto))[-train]^2)
mean((mpg - predict(lm.fit2, Auto))[-train]^2)
mean((mpg - predict(lm.fit3, Auto))[-train]^2)


##### 5.3.2 LOOCV #####
coef(glm.fit <- glm(mpg ~ horsepower, data = Auto))
coef(lm.fit <- lm(mpg ~ horsepower, data = Auto))

library(boot)
glm.fit <-glm(mpg ~ horsepower, data = Auto)
cv.err <- cv.glm(Auto, glm.fit) # default K = n (LOOCV)
names(cv.err)
cv.err$delta 
# [1] standard CV estimate [2] bias-corrected

cv.error <- rep(0, 5)
for (i in 1:5) {
  glm.fit <- glm(mpg ~ poly(horsepower, i), data = Auto)
  cv.error[i] <- cv.glm(Auto, glm.fit)$delta[1]
}
cv.error 
# no clear improvement after quadratic fits
plot(cv.error, type = 'l')


##### 5.3.3 K-Fold CV #####
set.seed(17)
cv.error.10 <- rep(0, 10)
for (i in 1:10) {
  glm.fit <- glm(mpg ~ poly(horsepower, i), data = Auto)
  cv.error.10[i] <- cv.glm(Auto, glm.fit, K = 10)$delta[1]
}
plot(cv.error.10, type = 'l')


##### 5.3.4 Bootstrap #####

alpha.fn <- function(data, index) {
  X <- data$X[index]
  Y <- data$Y[index]
  alpha <- (var(Y)-cov(X,Y)) / (var(X)+var(Y)-2*cov(X,Y))
  return(alpha)
}
alpha.fn(Portfolio, 1:100)

set.seed(1)
alpha.fn(Portfolio, sample(100, 100, replace = T)) # BS 1
boot(Portfolio, alpha.fn, R = 1000) # 1000 Bootstrap Samples 

boot.fn <- function(data, index) { 
  coef(lm(mpg ~ horsepower, data = data, subset = index)) 
}
boot.fn(Auto, 1:392)

set.seed(1)
boot.fn(Auto, sample(392, 392, replace = T)) # BS1
boot.fn(Auto, sample(392, 392, replace = T)) # BS2 
boot(Auto, boot.fn, R = 1000) # 1000 BS's
summary(lm(mpg ~ horsepower, data = Auto))$coefficients # real

boot.fn2 <- function(data, index) {
  coef(lm(mpg ~ horsepower + I(horsepower^2), data = data, subset = index))
}  
boot(Auto, boot.fn2, R = 1000)
summary(lm(mpg ~ horsepower + I(horsepower^2), data = Auto))$coefficients

