library(ISLR)
library(tree)
library(randomForest)
library(gbm)
library(glmnet)
library(tidyverse)

##### CH8 Exercises

##### Quick Functions
regFitPlot <- function(pred, true) {
  mse <- mean((true - pred)^2)
  plot(pred, true, pch = 19, col = 'gray', xlab = 'predicted',
       main = glue('test MSE = {round(mse, 3)}'))
  abline(0, 1, col = 'steelblue', lwd = 2)
  return(mse)
}

##### Exercise 8.8 #####
data(Carseats)
dim(Carseats)

# (a) split the dataset
set.seed(2)
train <- sample(nrow(Carseats), nrow(Carseats)/2)
carseats.test <- Carseats[-train,]
sales.test <- Carseats$Sales[-train]

# (b) fit a regression tree: plot, test MSE
summary(tree.carseats <- tree(Sales ~ ., Carseats, subset = train))  
plot(tree.carseats, col = 'gray')
text(tree.carseats, digits = 3, pretty = 0, font = 10)

pred.tree <- predict(tree.carseats, carseats.test)
(mse.tree <- regFitPlot(pred.tree, sales.test))

# (c) cv.tree
set.seed(20)
(cv.tree.carseats <- cv.tree(tree.carseats, FUN = prune.tree))

cv.table <- data.frame(cv.tree.carseats[1], cv.tree.carseats[2], cv.tree.carseats[3])
cv.table[which.min(cv.table$dev),]
plot(cv.tree.carseats)
ggplot(cv.table, aes(size, dev)) + theme_bw() +
  geom_point() + geom_line(col = 'gray') +
  geom_hline(aes(yintercept = min(dev)), linetype = 2, col = 'steelblue')
# 13-node tree
# pruning can improve the test MSE.
pred.tree <- predict(prune.tree(tree.carseats, best = 13), carseats.test)
(mse.tree <- regFitPlot(pred.tree, sales.test))

# (d) bagging
p <- ncol(Carseats) -1
set.seed(2)
(bag.carseats <- randomForest(Sales ~ ., Carseats, subset = train, mtry = p, importance = TRUE))
importance(bag.carseats)
varImpPlot(bag.carseats, pch = 19, main = 'Variable Importance')
# price의 Importance가 가장 높다 
pred.bag <- predict(bag.carseats, carseats.test)
(mse.bag <- regFitPlot(pred.bag, sales.test))

# (e) random forest (mtry = 3 = p/3)
set.seed(2)
(rf.carseats <- randomForest(Sales ~ ., Carseats, subset = train, importance = TRUE))
varImpPlot(rf.carseats, pch = 19, main = 'Variable Importance')

pred.rf <- predict(rf.carseats, carseats.test)
(mse.rf <- regFitPlot(pred.rf, sales.test))

c(mse.tree, mse.bag, mse.rf)


##### Exercise 8.9 #####
data(OJ)
dim(OJ)

# (a) split the data set
set.seed(9)
train <- sample(nrow(OJ), 800)
oj.test <- OJ[-train,]
purchase.test <- OJ$Purchase[-train]

# (b) fit a single tree
summary(tree.oj <- tree(Purchase ~ ., OJ, subset = train))

# (c) text output
tree.oj$frame

# (d) plot 
plot(tree.oj, col = 'gray')
text(tree.oj, digits = 3, pretty = 0, font = 10)

# (e) prediction
pred.tree <- predict(tree.oj, oj.test, type = 'class')
table(purchase.test, pred.tree)
mean(purchase.test != pred.tree)

# (f) cv.tree
set.seed(9)
(cv.tree.oj <- cv.tree(tree.oj, FUN = prune.misclass))

# (g) plot
cv.table <- data.frame(cv.tree.oj[1:3])
ggplot(cv.table, aes(size, dev)) + theme_bw() +
  geom_line(col = 'gray', size = 1) + geom_point(size = 2) + 
  geom_hline(aes(yintercept = min(dev)), linetype = 2, color = 'steelblue') +
  labs(y = 'deviance')

# (h) min error rate (deviance)
cv.table[which.min(cv.table$dev),] 

# (i)-(j)
summary(prune.oj <- prune.misclass(tree.oj, best = 6))

# (k) test error rate
pred.prune <- predict(prune.oj, oj.test, type = 'class')
table(purchase.test, pred.prune)
mean(purchase.test != pred.prune)


##### Exercise 8.10 Boosting
data(Hitters)
dim(Hitters)

# (a) 
sapply(Hitters, function(x) sum(is.na(x)))
hitters <- Hitters %>% na.omit %>% mutate(log.Salary = log(Salary)) %>% select(-Salary)

# (b)
train <- 1:200
hitters.train <- hitters[train,]
hitters.test <- hitters[-train,]
y.train <- hitters.train$log.Salary
y.test <- hitters.test$log.Salary

# (c)-(d)
lambda <- seq(0.001, 0.1, length = 100)
mse.lambda <- NULL
for (i in lambda) {
  bst.hitters <- gbm(
    log.Salary ~ ., data = hitters.train, 
    distribution = 'gaussian', n.trees = 1000, shrinkage = i
  )
  pred.train <- predict(bst.hitters, hitters.train, n.trees = 1000)
  pred.test <- predict(bst.hitters, hitters.test, n.trees = 1000)
  mse.train <- mean((y.train - pred.train)^2)
  mse.test <- mean((y.test - pred.test)^2)
  mse.lambda <- rbind(mse.lambda, c(mse.train, mse.test))
}
mse.lambda <- data.frame(lambda, mse.lambda)
colnames(mse.lambda) <- c('lambda', 'train', 'test')
mse.lambda %>% filter(test == min(test))
best.lmabda <- mse.lambda$lambda[which.min(mse.lambda$test)]

ggplot(mse.lambda) + theme_bw() +
  geom_line(aes(lambda, train, col = 'train'), size = 1) +
  geom_line(aes(lambda, test, col = 'test'), size = 1) +
  scale_color_manual('', values = c(train = 'gray50', test = 'red')) +
  labs(y = 'MSE') + theme(legend.position = 'top') +
  geom_vline(aes(xintercept = best.lmabda), linetype = 2)

# (e)
set.seed(10)
bst.hitters <- gbm(
  log.Salary ~ ., data = hitters.train, 
  distribution = 'gaussian', n.trees = 1000, shrinkage = best.lambda
)
pred.bst <- predict(bst.hitters, hitters.test, n.trees = 1000)

# linear regression (CH3)
summary(lm.hitters <- lm(log.Salary ~ ., hitters.train))
pred.lm <- predict(lm.hitters, hitters.test) %>% as.numeric

# lasso regression (CH6)
X.train <- model.matrix(lm.hitters)
X.test <- model.matrix(log.Salary ~ ., hitters.test)
lasso.hitters <- glmnet(X.train, y.train, alpha = 1)
set.seed(10)
cv.lasso <- cv.glmnet(X.train, y.train, alpha = 1)
pred.lasso <- predict(lasso.hitters, X.test, s = cv.lasso$lambda.min, exact = TRUE) %>% as.numeric

result <- data.frame(pred.bst, pred.lm, pred.lasso)
sapply(result, function(x) mean((y.test - x)^2))

# (f)
bst.tb <- summary(bst.hitters)
ggplot(bst.tb) + theme_bw() + coord_flip() +
  geom_col(aes(fct_reorder(var, (rel.inf)), rel.inf, fill = desc(rel.inf))) + 
  labs(x = 'variable') + theme(legend.position = '')

plot(bst.hitters, i = 'CAtBat', lwd = 2) # partial dependence plot

# (g) bagging
p <- ncol(hitters) -1
set.seed(10)
(bag.hitters <- randomForest(log.Salary ~ ., hitters.train, mtry = p, importance = TRUE))
pred.bag <- predict(bag.hitters, hitters.test)
mean((y.test - pred.bag)^2)

