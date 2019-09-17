library(ISLR)
library(MASS)
library(leaps)
library(glmnet)
library(pls)
library(ggplot2)
library(reshape2)

##### Exercise 6.9 #####
fix(College)
X <- data.matrix(College[,-2])
y <- College$Apps
MSE.College <- function(pred, test) mean((College$Apps[test]-pred)**2)

# (a) split the data set
set.seed(1)
train <- sample(1:nrow(College), nrow(College)/2)
test <- (-train)

# (b) least squares
(lmFit <- lm(Apps ~ ., College[train,]))
lmPred <- predict(lmFit, newdata = College[test,])
(mse.lm <- MSE.College(lmPred, test))
round(coef(lmFit), 3) 

# (c) ridge
set.seed(1)
(cv.ridge <- cv.glmnet(X[train,], y[train], alpha = 0))
(best.lambda.ridge <- cv.ridge$lambda.min)
ridgeFit <- glmnet(X[train,], y[train], alpha = 0)
ridgePred <- as.numeric(predict(ridgeFit, type = "response", s = best.lambda.ridge, newx = X[test,]))
(mse.ridge <- MSE.College(ridgePred, test))
ridgeCoef <- predict(ridgeFit, type = "coefficients", s = best.lambda.ridge, newx = X[test,])[1:ncol(X),]
round(ridgeCoef, 3)

# (d) lasso
set.seed(1)
(cv.lasso <- cv.glmnet(X[train,], y[train], alpha = 1))
(best.lambda.lasso <- cv.lasso$lambda.min)
lassoFit <- glmnet(X[train,], y[train], alpha = 1)
lassoPred <- as.numeric(predict(lassoFit, type = "response", s = best.lambda.lasso, newx = X[test,]))
(mse.lasso <- MSE.College(lassoPred, test))
lassoCoef <- predict(lassoFit, type = "coefficients", s = best.lambda.lasso, newx = X[test,])[1:ncol(X),]
round(lassoCoef, 3) 

# (e) PCR
set.seed(1)
summary(pcrFit <- pcr(Apps ~ ., data = College[train,], scale = TRUE, validation = "CV"))
validationplot(pcrFit, val.type = 'MSEP', main = 'Validation Plot of PCR') # M = 16

pcrFit <- pcr(Apps ~ ., data = College[train,], ncomp = 16)
pcrPred <- predict(pcrFit, X[test,], ncomp = 16)
(mse.pcr <- MSE.College(pcrPred, test))

# (f) PLS
set.seed(1)
summary(plsFit <- plsr(Apps ~ ., data = College[train,], scale = TRUE, validation = "CV"))
validationplot(plsFit, val.type = 'MSEP', main = 'Validation Plot of PLS') # M = 10

plsFit <- plsr(Apps ~ ., data = College[train,], ncomp = 10)
plsPred <- predict(plsFit, X[test,], ncomp = 10)
(mse.pls <- MSE.College(plsPred, test))

# (g)
mse.summary <- data.frame(method = c('OLS', 'Ridge', 'Lasso', 'PCR', 'PLS'),
                          MSE = c(mse.lm, mse.ridge, mse.lasso, mse.pcr, mse.pls))
mse.summary$method <- factor(mse.summary$method, mse.summary$method)
ggplot(mse.summary) + theme_bw() + scale_fill_brewer(palette = 'Pastel1') +
  geom_col(aes(method, MSE, fill = method), color = 'black')
ggplot() +theme_bw() + geom_point(aes(lassoPred, y[test]), color = 'gray50', size = 2) + 
  geom_abline(aes(intercept = 0, slope = 1), col = 'lightblue', size = 1) +
  labs(x = 'Predicted', y = 'True y')


##### Exercise 6.10 #####
fix(Boston)
X <- data.matrix(Boston[,-14])
y <- Boston$medv
MSE.Boston <- function(pred, test) mean((Boston$medv[test]-pred)^2)

# (a)
set.seed(10)
train <- sample(1:nrow(Boston), nrow(Boston)/2)
test <- (-train)

# OLS
lmFit <- lm(medv ~ ., Boston[train,])
lmPred <- predict(lmFit, Boston[test,])
(mse.lm <- MSE.Boston(lmPred, test))
round(coef(lmFit), 3) 

# Best Subset
predict.regsubsets <- function(object, newdata, id, ...) {
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <- coef(object, id = id)
  xvars <- names(coefi)
  mat[,xvars] %*% coefi
}
# stepwise
regFit <- regsubsets(medv ~ ., Boston[train,], nvmax = 13)
regSummary <- summary(regFit)
par(mfrow = c(1,4))
plot(regSummary$rss, type = 'l', pch = 19, xlab = 'p', ylab = 'RSS', main = 'RSS')
which.min(regSummary$rss)
points(13, regSummary$rss[13], col = 'blue', pch = 16)
plot(regSummary$adjr2, type = 'l', pch = 19, xlab = 'p', ylab = 'Adj.Rsq', main = 'Adj.Rsq')
which.max(regSummary$adjr2)
points(11, regSummary$adjr2[11], col = 'blue', pch = 16)
plot(regSummary$cp, type = 'l', pch = 19, xlab = 'p', ylab = 'Cp', main = 'Cp')
which.min(regSummary$cp)
points(11, regSummary$cp[11], col = 'blue', pch = 16)
plot(regSummary$bic, type = 'l', pch = 19, xlab = 'p', ylab = 'BIC', main = 'BIC')
which.min(regSummary$bic)
points(7, regSummary$bic[7], col = 'blue', pch = 16)

plot(regFit, scale = 'Cp')
coef(regFit, 11)

dtrain <- Boston[train,]
k <- 10
set.seed(1)
folds <- sample(1:k, nrow(dtrain), replace = TRUE)
cv.errors <- matrix(NA, k, 13, dimnames = list(NULL, paste(1:13)))
for(j in 1:k) {
  best.fit <- regsubsets(medv ~ ., data = dtrain[folds!=j,], nvmax = 13)
  for(i in 1:13) {
    pred <- predict(best.fit, dtrain[folds==j,], id = i) 
    cv.errors[j,i] <- mean((dtrain$medv[folds==j]-pred)^2)
  }
}
(mean.cv.errors <- apply(cv.errors, 2, mean))

par(mfrow = c(1,1))
plot(mean.cv.errors, type = 'b', pch = 19, xlab = 'p', main = 'CV Error')
which.min(mean.cv.errors)
points(11, mean.cv.errors[11], col = 'blue', pch = 19)
abline(v = 11, lty = 2)

regPred <- predict(regFit, Boston[test,], 11)
(mse.best <- MSE.Boston(regPred, test))
round(coef(regFit, 11), 3) 

# Ridge
set.seed(10)
(cv.ridge <- cv.glmnet(X[train,], y[train], alpha = 0))
(best.lambda.ridge <- cv.ridge$lambda.min)
ridgeFit <- glmnet(X[train,], y[train], alpha = 0)
ridgePred <- as.numeric(predict(ridgeFit, type = "response", s = best.lambda.ridge, newx = X[test,]))
(mse.ridge <- MSE.Boston(ridgePred, test))

ridgeCoef <- (predict(ridgeFit, type = "coefficients", s = best.lambda.ridge, newx = X[test,]))[1:14,]
round(ridgeCoef, 3)

# Lasso
set.seed(10)
(cv.lasso <- cv.glmnet(X[train,], y[train], alpha = 1))
(best.lambda.lasso <- cv.lasso$lambda.min)
lassoFit <- glmnet(X[train,], y[train], alpha = 1)
lassoPred <- as.numeric(predict(lassoFit, type = "response", s = best.lambda.lasso, newx = X[test,]))
(mse.lasso <- MSE.Boston(lassoPred, test))

lassoCoef <- (predict(lassoFit, type = "coefficients", s = best.lambda.lasso, newx = X[test,]))[1:14,]
round(lassoCoef, 3)

# PCR
set.seed(10)
pcrFit <- pcr(medv ~ ., data = Boston[train,], scale = TRUE, validation = "CV")
summary(pcrFit)
validationplot(pcrFit, val.type = 'MSEP', main = 'Validation Plot of PCR') # M = 13

pcrFit <- pcr(medv ~ ., data = Boston[train,], ncomp = 13)
pcrPred <- predict(pcrFit, X[test,], ncomp = 13)
(mse.pcr <- MSE.Boston(pcrPred, test))

# PLS
set.seed(10)
plsFit <- plsr(medv ~ ., data = Boston[train,], scale = TRUE, validation = "CV")
summary(plsFit)
validationplot(plsFit, val.type = 'MSEP', main = 'Validation Plot of PLS') # M = 9

plsFit <- plsr(medv ~ ., data = Boston[train,], ncomp = 9)
plsPred <- predict(plsFit, X[test,], ncomp = 9)
(mse.pls <- MSE.Boston(plsPred, test))

# (b)
mse.summary <- data.frame(method = c('OLS', 'BestSub', 'Ridge', 'Lasso', 'PCR', 'PLS'),
                          MSE = c(mse.lm, mse.best, mse.ridge, mse.lasso, mse.pcr, mse.pls))
mse.summary$method <- factor(mse.summary$method, mse.summary$method)
ggplot(mse.summary) + theme_bw() + scale_fill_brewer(palette = 'Pastel1') +
  geom_col(aes(method, MSE, fill = method), color = 'black')
ggplot() +theme_bw() + geom_point(aes(lmPred, y[test]), color = 'gray50', size = 2) + 
  geom_abline(aes(intercept = 0, slope = 1), col = 'lightblue', size = 1) +
  labs(x = 'Predicted', y = 'True y')


