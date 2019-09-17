
##### 9.6 Lab: Support Vector Machines #####
##### 9.6.1 Support Vector Classifier (linear)
library(e1071)

set.seed(1)
x <- matrix(rnorm(20*2), ncol = 2)
y <- c(rep(-1, 10), rep(1, 10))
x[y==1,] <- x[y==1,] + 1
plot(x, col = (3-y))

dat <- data.frame(x, y = as.factor(y))
svm.fit <- svm(y ~ ., data = dat, kernel = 'linear', cost = 10, scale = FALSE)
# scale: not to scale each feature (generally use scale = TRUE)

plot(svm.fit, dat)
svm.fit$index # support vectors (x shape in the plot)
summary(svm.fit)

# smaller cost -> margin gets wider
svm.fit <- svm(y ~ ., data = dat, kernel = 'linear', cost = 0.1, scale = FALSE)
plot(svm.fit, dat)
svm.fit$index # more support vectors

# cross-validation
set.seed(1)
tune.out <- tune(svm, y ~ ., data = dat, kernel = 'linear',
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1.5, 10, 100)))
summary(tune.out)
best.mod <- tune.out$best.model
summary(best.mod)

xtest <- matrix(rnorm(20*2), ncol = 2)
ytest <- sample(c(-1, 1), 20, replace = TRUE)
xtest[ytest==1,] <- xtest[ytest==1,] + 1
dat.test <- data.frame(xtest, y = as.factor(ytest))

ypred <- predict(best.mod, dat.test)
table(predict = ypred, true = ytest) ; mean(ypred != ytest)


svm.fit <- svm(y ~ ., data = dat, kernel = 'linear', cost = 0.01, scale = FALSE)
ypred <- predict(svm.fit, dat.test)
table(predict = ypred, true = ytest) ; mean(ypred != ytest)

# 2-classes linearly separable
x[y==1,] <- x[y==1,] + 0.5
plot(x, col = (y+5)/2, pch = 19)
dat <- data.frame(x, y = as.factor(y))
svm.fit <- svm(y ~ ., data = dat, kernel = 'linear', cost = 1e5)
summary(svm.fit)
table(true = dat$y, fitted = svm.fit$fitted)

svm.fit <- svm(y ~ ., data = dat, kernel = 'linear', cost = 1)
summary(svm.fit)
table(true = dat$y, fitted = svm.fit$fitted)
plot(svm.fit, dat)


##### 9.6.2 Support Vector Machine (radial, polynomial)
set.seed(1)
x <- matrix(rnorm(200*2), ncol = 2)
x[1:100,] <- x[1:100,] + 2
x[101:150,] <- x[101:150,] - 2 
y <- c(rep(1, 150), rep(2, 50))
dat <- data.frame(x, y = as.factor(y))

plot(x, col = y, pch = 19)

train <- sample(200, 100)
svm.fit <- svm(y ~ ., dat[train,], kernel = 'radial', gamma = 1, cost = 1)
plot(svm.fit, dat[train,])
summary(svm.fit)

svm.fit <- svm(y ~ ., dat[train,], kernel = 'radial', gamma = 1, cost = 1e5)
plot(svm.fit, dat[train,])

set.seed(1)
tune.out <- tune(svm, y ~ ., data = dat[train,], kernel = 'radial',
                 range = list(cost = c(0.1, 1, 10, 100, 100), gamma = c(0.5, 1:4)))
summary(tune.out)
summary(tune.out$best.model) # cost = 1, gamma = 0.5

ypred <- predict(tune.out$best.model, dat[-train,])
table(true = dat[-train,'y'], predicted = ypred)
mean(ypred != dat[-train,'y'])


##### 9.6.3 ROC Curves 
library(ROCR)
roc.plot <- function(pred, true, ...) {
  predob <- prediction(pred, true)
  perf <- performance(predob, 'tpr', 'fpr')
  plot(perf, ...)
}

svm.opt <- svm(y ~., data = dat[train,], kernel = 'radial',
               gamma = 0.5, cost = 1, decision.values = TRUE)
fitted <- attributes(predict(svm.opt, dat[train,], decision.values = TRUE))$decision.values

par(mfrow = c(1, 2))
roc.plot(fitted, dat[train,'y'], main = 'Training data')

# increase gamma for a flexible fit (smooth)
svm.flex <- svm(y ~., data = dat[train,], kernel = 'radial',
               gamma = 50, cost = 1, decision.values = TRUE)
fitted <- attributes(predict(svm.flex, dat[train,], decision.values = TRUE))$decision.values
roc.plot(fitted, dat[train,'y'], main = 'Training data', add = TRUE, col = 'red')

ypred <- predict(svm.opt, dat[-train,], decision.values = TRUE)
fitted <- attributes(ypred)$decision.values
roc.plot(fitted, dat[-train,'y'], main = 'Test data')
ypred <- predict(svm.flex, dat[-train,], decision.values = TRUE)
fitted <- attributes(ypred)$decision.values
roc.plot(fitted, dat[-train,'y'], main = 'Test data', add = TRUE, col = 'red')


##### 9.6.4 SVM with Multi-Classes
# add 1 class
set.seed(1)
x <- rbind(x, matrix(rnorm(50*2), ncol = 2))
y <- c(y, rep(0, 50))
x[y==0, 2] <- x[y==0, 2] + 2
dat <- data.frame(x = x, y = as.factor(y))

par(mfrow = c(1, 1))
plot(x, col = (y+1), pch = 19)

svm.fit <- svm(y ~ ., dat, kernel = 'radial', cost = 10, gamma = 1)
plot(svm.fit, dat)
table(true = dat$y, fitted = svm.fit$fitted)


##### 9.6.5 Application to Gene Expression Data
data(Khan, package = 'ISLR')
names(Khan)
dim(Khan$xtrain) ; dim(Khan$xtest)
length(Khan$ytrain) ; length(Khan$ytest)

table(Khan$ytrain)
table(Khan$ytest)

dat <- data.frame(x = Khan$xtrain, y = as.factor(Khan$ytrain))
out <- svm(y ~ ., dat, kernel = 'linear', cost = 10)
summary(out)
table(true = dat$y, fitted = out$fitted) # training error = 0

# when n << p, it is easy to find hyperplanes that fully separate the classes.

dat.test <- data.frame(x = Khan$xtest, y = as.factor(Khan$ytest)) 
ypred.test <-predict(out, dat.test)
table(true = dat.test$y, predicted = ypred.test) # 2 misclassification
