library(tree)
library(ISLR)
library(MASS)
library(randomForest)
library(gbm)

##### Chapter 8: Tree-based Methods #####
##### 8.3  Lab: Decision Trees #####
##### 8.3.1 Fitting 'Classification' Trees
attach(Carseats)
High <- ifelse(Sales <= 8, 'No', 'Yes') 
Carseats <- data.frame(Carseats, High) 

summary(tree.carseats <- tree(High ~ . -Sales, Carseats))
plot(tree.carseats)
text(tree.carseats, pretty = 0)

set.seed(2)
train <- sample(nrow(Carseats), 200)
Carseats.test <- Carseats[-train,]
High.test <- High[-train]

tree.carseats <- tree(High ~ . -Sales, Carseats, subset = train)
tree.pred <- predict(tree.carseats, Carseats.test, type = 'class')
table(tree.pred, High.test)

# cross-validation 
# to determine the optimal level of tree complexity

set.seed(3)
(cv.carseats <- cv.tree(tree.carseats, FUN = prune.misclass))
names(cv.carseats)
# size: number of terminal nodes of each tree considered
# dev : total deviance of each tree in the cost-complexity pruning sequence
# k   : cost-complexity parameter used

par(mfrow = c(1, 2))
plot(cv.carseats$size, cv.carseats$dev, type = 'b')
plot(cv.carseats$k, cv.carseats$dev, type = 'b')

# 9-node tree를 얻기 위한 가지치기
prune.carseats <- prune.misclass(tree.carseats, best = 9)
par(mfrow = c(1, 1))
plot(prune.carseats)
text(prune.carseats, pretty = 0)

tree.pred <- predict(prune.carseats, Carseats.test, type = 'class')
table(tree.pred, High.test)


##### 8.3.2 Fitting 'Regression' Trees
set.seed(1)
train <- sample(nrow(Boston), nrow(Boston)/2)
summary(tree.boston <- tree(medv ~ ., Boston, subset = train))  
# deviance = Sum of Squared Errors
plot(tree.boston)
text(tree.boston, pretty = 0)

cv.boston <- cv.tree(tree.boston) # default FUN = prune.tree
plot(cv.boston$size, cv.boston$dev, type = 'b') # 가장 복잡한 모형이 선택됨

prune.boston <- prune.tree(tree.boston, best = 5) # 5-node tree
plot(prune.boston)
text(prune.boston, pretty = 0)

yhat <- predict(tree.boston, Boston[-train,])
boston.test <- Boston[-train, 'medv']
plot(yhat, boston.test)
abline(0, 1)
mean((yhat - boston.test)^2) # MSE


##### 8.3.3 Bagging and Random Forest
# Bagging
set.seed(1)
(bag.boston <- randomForest(medv ~ ., data = Boston, subset = train, mtry = 13, importance = TRUE))
# bagging에서는 mtry = p (전체 설명변수를 모두 split variable로 사용)
yhat.bag <- predict(bag.boston, Boston[-train,])
plot(yhat.bag, boston.test)
abline(0, 1)
mean((yhat.bag - boston.test)^2) 

# Random Forest
set.seed(1)
(rf.boston <- randomForest(medv ~ ., data = Boston, subset = train, mtry = 6, importance = TRUE))
yhat.rf <- predict(rf.boston, Boston[-train,])
plot(yhat.rf, boston.test)
abline(0, 1)
mean((yhat.rf - boston.test)^2) 

importance(rf.boston)
varImpPlot(rf.boston)


##### 8.3.4 Boosting
set.seed(1)
summary(boost.boston <- gbm(
  medv ~ ., data = Boston[train,], distribution = 'gaussian', n.trees = 5000, interaction.depth = 4
))

# partial dependence plots for some important variables
par(mfrow = c(1, 2))
plot(boost.boston, i = 'rm')    # rm이 높을수록 y가 높아짐
plot(boost.boston, i = 'lstat') # lstat이 높을수록 y가 작아짐

yhat.boost <- predict(boost.boston, newdata = Boston[-train,], n.trees = 5000)
par(mfrow = c(1, 1))
plot(yhat.boost, boston.test)
abline(0, 1)
mean((yhat.boost - boston.test)^2)


