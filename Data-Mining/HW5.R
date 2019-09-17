library(ISLR)
library(boot)
library(MASS)
library(dplyr)
library(ggplot2)

##### 5.2 #####
# (d)-(f)
prob <- function(n) 1 - (1-1/n)^n
prob(5)
prob(100)
prob(10000)

# (g)
n <- 1:100000
ggplot() + theme_bw() + labs(y = 'probability') +
  geom_line(aes(n, prob(n)), color = 'gray', size = 1) + 
  geom_point(aes(n, prob(n)), color = 'gray50', size = 1)
summary(prob(n))  

# (h)
store <- rep(0, 10000)
for (i in 1:10000) store[i] <- sum(sample(1:100, replace = TRUE) == 4) > 0
mean(store)



##### 5.5 #####
# (a) Logistic Regression
summary(glm5fit <- glm(default ~ income + balance, Default, family = binomial))

# (b) Validation Approach
# i. split
set.seed(10)
train.i <- sample(1:10000, 5000)
train <- Default[train.i,]
val <- Default[-train.i,]
# ii. fit
summary(glm5fit <- glm(default ~ income + balance, train, family = binomial))
# iii. CV
glm5prob <- predict(glm5fit, val, type = 'response')
glm5pred <- ifelse(glm5prob > 0.5, 'Yes', 'No')
table(glm5pred, val$default)
mean(glm5pred != val$default)

# (c) Repeat 3 times
repeatCV <- function() {
  train.i <- sample(1:10000, 5000)
  train <- Default[train.i,]
  val <- Default[-train.i,]
  glm5fit <- glm(default ~ income + balance, train, family = binomial)
  glm5prob <- predict(glm5fit, val, type = 'response')
  glm5pred <- ifelse(glm5prob > 0.5, 'Yes', 'No')
  return(mean(glm5pred != val$default))
}
repeatCV()

# (d) 
summary(glm5fit.d <- glm(default ~ student + income + balance, train, family = binomial))
glm5prob.d <- predict(glm5fit.d, val, type = 'response')
glm5pred.d <- ifelse(glm5prob.d > 0.5, 'Yes', 'No')
mean(glm5pred.d != val$default)
table(glm5pred.d, val$default)



##### 5.7 #####
# (a)
summary(glm7fit <- glm(Direction ~ Lag1 + Lag2, Weekly, family = binomial))
# (b)
summary(glm7fit.b <- glm(Direction ~ Lag1 + Lag2, Weekly, family = binomial, subset = -1))
predict(glm7fit.b, Weekly[1,], type = 'response')
Weekly[1,'Direction']
# (c)
store <- rep(0, nrow(Weekly))
for (i in 1:nrow(Weekly)) {
  fit <- glm(Direction ~ Lag1 + Lag2, Weekly, family = binomial, subset = -i)
  pred <- ifelse(predict(fit, Weekly[i,], type = 'response') > 0.5, 'Up', 'Down') 
  store[i] <- as.numeric(pred != Weekly[i,'Direction'])
}
# (d)
mean(store)
myloocv <- cv.glm(Weekly, glm7fit)
myloocv$delta



##### 5.9 #####
# (a)-(b)
attach(Boston)
mean(medv)
sd(medv) / sqrt(length(medv))

# (c)
boot9c <- function(data,index) mean(data[index])
(result9c <- boot(medv, boot9c, R = 1000))

# (d) 95% CI
c(result9c$t0 - 2*sd(result9c$t), result9c$t0 + 2*sd(result9c$t))
t.test(medv)

# (e)
median(medv)
# (f)
boot9f <- function(data,index) median(data[index])
(result9f <- boot(medv, boot9f, R = 1000))

# (g)
quantile(medv, 0.1)
# (h)
boot9h <- function(data, index) quantile(data[index], 0.1)
(result9h <- boot(medv, boot9h, R = 1000))
