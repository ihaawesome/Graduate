library(e1071)
library(ROCR)
library(ISLR)
library(tidyverse)

myggplot <- function(...) {
  ggplot(...) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_color_brewer(palette = 'Set1') 
}
##### Exercise 9.4 #####

# Generate a simulated 2-class data set with 100 observations and 2 features.
set.seed(4)
mydata <- data.frame(matrix(rnorm(100*2), ncol = 2)) %>% 
  mutate(y = as.factor(ifelse(X1^2 + X2^2 < 1.5, 1, 2)))

myggplot(mydata, legend.name = 'class') + geom_point(aes(X1, X2, color = y, shape = y)) 
table(mydata$y)

# show that a SVM with a polynomial kernel or a radial kernel will outperform a classifier.
set.seed(4)
train <- sample(100, 50)

svm.linear <- svm(y ~ ., mydata[train,], kernel = 'linear', cost = 10, gamma = 1)
summary(svm.linear)
plot(svm.linear, mydata[train,])
table(true = mydata[train,'y'], pred = svm.linear$fitted)
mean(mydata[train,'y'] != svm.linear$fitted) # train error = 20%

svm.radial <- svm(y ~ ., mydata[train,], kernel = 'radial', cost = 10, gamma = 1)
summary(svm.radial)
plot(svm.radial, mydata[train,])
table(true = mydata[train,'y'], pred = svm.radial$fitted)
mean(mydata[train,'y'] != svm.radial$fitted) # train error = 0%

svm.poly <- svm(y ~ ., mydata[train,], kernel = 'polynomial', cost = 10, gamma = 1)
summary(svm.poly)
plot(svm.poly, mydata[train,])
table(true = mydata[train,'y'], pred = svm.poly$fitted)
mean(mydata[train,'y'] != svm.poly$fitted) # train error = 34%

# performance: radial > linear >= polynomial
pred.train <- data.frame(
  linear = svm.linear$fitted,
  radial = svm.radial$fitted,
  poly = svm.poly$fitted
)
apply(pred.train, 2, function(pred) mean(mydata[train,'y'] != pred))

pred.test <- data.frame(
  linear = predict(svm.linear, mydata[-train,]),
  radial = predict(svm.radial, mydata[-train,]),
  poly = predict(svm.poly, mydata[-train,])
)
apply(pred.test, 2, function(pred) mean(mydata[-train,'y'] != pred))


##### Exercise 9.5 #####
# (a) Generate a data set with n=500, p=2, 
# such that the observations belong to 2 classes with a quadratic decision boundary between them.
set.seed(500)
mydata <- data.frame(X1 = runif(500) - 0.5, X2 = runif(500) - 0.5) %>%
  mutate(y = as.factor(1*(X1^2 - X2^2 > 0)))

# (b) Plot
myggplot(mydata) + geom_point(aes(X1, X2, col = y, shape = y))

# (c) Logistic Regression
summary(myglm <- glm(y ~ ., mydata, family = binomial))
# (d) Plot (colored by predicted)
pred.glm <- as.factor(as.numeric(predict(myglm, type = 'response') > 0.5))
myggplot(mydata) + geom_point(aes(X1, X2, col = pred.glm, shape = y))

# (e)
summary(myglm.poly <- glm(y ~ poly(X1, 2) + poly(X2, 2), mydata, family = binomial))
# (f)
pred.glm <- as.factor(as.numeric(predict(myglm.poly, type = 'response') > 0.5))
myggplot(mydata) + geom_point(aes(X1, X2, col = pred.glm, shape = y))

# (g)
summary(mysvm.linear <- svm(y ~ ., mydata, kernel = 'linear'))
myggplot(mydata, legend.name = 'linear') + 
  geom_point(aes(X1, X2, col = predict(mysvm.linear), shape = y))
 
# (h) 
summary(mysvm.poly <- svm(y ~ ., mydata, kernel = 'polynomial'))
myggplot(mydata, legend.name = 'polynomial') + geom_point(aes(X1, X2, col = predict(mysvm.poly), shape = y)) 
summary(mysvm.radial <- svm(y ~ ., mydata, kernel = 'radial'))
myggplot(mydata, legend.name = 'radial') + geom_point(aes(X1, X2, col = predict(mysvm.radial), shape = y))

# (i)
pred.table <- data.frame(
  glm.linear = as.factor(as.numeric(predict(myglm, type = 'response') > 0.5)),
  glm.poly = as.factor(as.numeric(predict(myglm.poly, type = 'response') > 0.5)),
  svm.linear = predict(mysvm.linear),
  svm.poly = predict(mysvm.poly),
  svm.radial = predict(mysvm.radial)
)
apply(pred.table, 2, function(pred) mean(mydata$y != pred))


##### Exercise 9.8
data(OJ)
# (a)
set.seed(8)
train <- sample(nrow(OJ), 800)

# (b) Fit a support vector classifier using cost=0.01
summary(svm.linear <- svm(Purchase ~ ., OJ[train,], kernel = 'linear', cost = 0.01))

# (c) training and test error rates
misclassification <- function(svm.model) {
  out1 <- table(true = OJ[train,'Purchase'], predicted.train = predict(svm.model, OJ[train,]))
  out2 <- mean(OJ[train,'Purchase'] != predict(svm.model, OJ[train,]))
  out3 <- table(true = OJ[-train,'Purchase'], predicted.test = predict(svm.model, OJ[-train,]))
  out4 <- mean(OJ[-train,'Purchase'] != predict(svm.model, OJ[-train,]))
  return(list(train.mat = out1, train.error = out2, test.mat = out3, test.error = out4))
}
misclassification(svm.linear)

# (d) tune to select an optimal cost
set.seed(8)
tune.cost <- tune(svm, Purchase ~ ., data = OJ[train,], kernel = 'linear',
                  ranges = list(cost = c(0.01, seq(0.1, 10, by = 0.5))))
tune.cost
summary(tune.cost)
plot(tune.cost$performances[,1:2], type = 'l', pch = 19, lwd = 2)
abline(h = min(tune.cost$performance[,2]), lty = 2, col = 'blue')
# (e) predict using new value for cost
best.linear <- tune.cost$best.model
misclassification(best.linear)

# (f) repeat using radial kernel
summary(svm.radial <- svm(Purchase ~ ., OJ[train,], kernel = 'radial', cost = 0.01))
misclassification(svm.radial)

set.seed(8)
tune.cost.rad <- tune(svm, Purchase ~ ., data = OJ[train,], kernel = 'radial',
                  ranges = list(cost = c(0.01, seq(0.1, 10, by = 0.5))))
tune.cost.rad
summary(tune.cost.rad)
plot(tune.cost.rad$performances[,1:2], type = 'l', pch = 19, lwd = 2)
best.radial <- tune.cost.rad$best.model
misclassification(best.radial)

# (g) repat using polynomial kernel with degree=2
summary(svm.poly <- svm(Purchase ~ ., OJ[train,], kernel = 'polynomial', degree = 2, cost = 0.01))
misclassification(svm.poly)

set.seed(8)
tune.cost.poly <- tune(svm, Purchase ~ ., data = OJ[train,], kernel = 'polynomial', degree = 2,
                      ranges = list(cost = c(0.01, seq(0.1, 10, by = 0.5))))
tune.cost.poly
summary(tune.cost.poly)
plot(tune.cost.poly$performances[,1:2], type = 'l', pch = 19, lwd = 2)
best.poly <- tune.cost.poly$best.model
misclassification(best.poly)

# (h)
model.list <- list(svm.linear, svm.radial, svm.poly, best.linear, best.radial, best.poly)
data.frame(row.names = c('linear','radial','poly','tuned.liniear','tuned.radial','tuned.poly'),
  train = unlist(lapply(model.list, function(obj) misclassification(obj)[[2]])),
  test = unlist(lapply(model.list, function(obj) misclassification(obj)[[4]]))
)
tune.cost$best.performance
tune.cost.rad$best.performance 
tune.cost.poly$best.performance
