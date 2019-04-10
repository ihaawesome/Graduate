library(ISLR)
library(class)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)

########## Chapter 4 Exercises ##########

##### Exercise 4.4 #####
meanProp <- function(X) {
  n <- nrow(X) ; p <- ncol(X) ; prop <- matrix(nrow = n, ncol = p)
  for (i in 1:n) for(j in 1:p) 
    prop[i,j] <- mean(between(X[,p], X[i,p]-0.05, X[i,p]+0.05))
  prop <- apply(prop, 1, prod)
  return(mean(prop))
}

n <- 1000

# (a) p = 1
set.seed(1)
X1 <- matrix(runif(n), ncol = 1)
meanProp(X1)

# (b) p = 2
set.seed(1)
X2 <- matrix(runif(n*2, 0, 1), ncol = 2)
meanProp(X2)

# (c) p = 100
set.seed(1)
X100 <- matrix(runif(n*100, 0, 1), ncol = 100)
meanProp(X100)

# (e)
(0.1) # p = 1
(0.1)^(1/2) # p = 2
(0.1)^(1/100) # p = 100



##### Exercise 4.10 #####
# (a)
dim(Weekly)
apply(Weekly[,-9], 2, quantile, probs = c(0, 0.5, 1))
apply(Weekly[,-9], 2, mean)
apply(Weekly[,-9], 2, sd)
table(Weekly[,9])
cor(Weekly[,-9])
pairs(Weekly[,-9], upper.panel = NULL, col = 'gray')

gather(Weekly, var, value, -Direction) %>%
  mutate(var = factor(var, levels = names(Weekly)[-9])) %>%
  ggplot() + theme_bw() +
  geom_boxplot(aes(x = Direction, y = value, group = Direction)) + 
  facet_wrap(~var, ncol = 4, scales = 'free')

# (b) GLM Full
glmFit <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume,
              data = Weekly, family = binomial)
summary(glmFit) # Lag2 

# (c) 
glmProb <- predict(glmFit, type = 'response')
glmPred <- ifelse(glmProb > 0.5, 'Up', 'Down')
table(glmPred, Weekly$Direction)
mean(glmPred == Weekly$Direction)

# (d) Split & only Lag2
train <- Weekly %>% filter(Year %in% (1990:2008))
test <- Weekly %>% filter(Year %in% (2009:2010))

glmFit.d <- glm(Direction ~ Lag2, data = train, family = binomial)
summary(glmFit.d)

glmProb.d <- predict(glmFit.d, newdata = test, type = 'response')
glmPred.d <- ifelse(glmProb.d > 0.5, 'Up', 'Down')
table(glmPred.d, test$Direction)
mean(glmPred.d == test$Direction)

# (e) LDA
ldaFit <- lda(Direction ~ Lag2, data = train)
ldaPred <- predict(ldaFit, newdata = test)$class
table(ldaPred, test$Direction)
mean(ldaPred == test$Direction)

# (f) QDA
qdaFit <- qda(Direction ~ Lag2, data = train)
qdaPred <- predict(qdaFit, newdata = test)$class
table(qdaPred, test$Direction)
mean(qdaPred == test$Direction)

# (g) KNN
trainX <- as.matrix(train$Lag2)
testX <- as.matrix(test$Lag2)
trainY <- as.matrix(train$Direction)

set.seed(1)
knnPred <- knn(trainX, testX, trainY, k = 1)
table(knnPred, test$Direction)
mean(knnPred == test$Direction)

# (i)
glmFit.i <- step(glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume,
                     data = train, family = binomial))
glmProb.i <- predict(glmFit.i, test, type = 'response')
glmPred.i <- ifelse(glmProb.i > 0.5, "Up", "Down")
mean(glmPred.i == test$Direction)

glmFit.i <- update(glmFit.i, ~ Lag1*Lag2)
glmProb.i <- predict(glmFit.i, test, type = 'response')
glmPred.i <- ifelse(glmProb.i > 0.5, "Up", "Down")
mean(glmPred.i == test$Direction)

glmFit.i <- update(glmFit.i, ~ poly(Lag2, 2))
glmProb.i <- predict(glmFit.i, test, type = 'response')
glmPred.i <- ifelse(glmProb.i > 0.5, "Up", "Down")
mean(glmPred.i == test$Direction)

ldaFit.i <- update(ldaFit, ~ I(Lag1^2) + Lag2)
ldaPred.i <- predict(ldaFit.i, newdata = test)$class
mean(ldaPred.i == test$Direction)

# KNN
p <- 100
for (i in 2:p) {
  set.seed(1)
  tmp <- knn(trainX, testX, trainY, k = i)
  knnPred <- data.frame(knnPred, tmp)
}
colnames(knnPred) <- paste0('p', 1:p)

knnAcc <- NULL
for(i in 1:p) knnAcc <- c(knnAcc, mean(knnPred[,i] == test$Direction))
which.max(knnAcc) ; knnAcc[which.max(knnAcc)]
ggplot()+ theme_bw() + 
geom_line(aes(1:p, knnAcc), col = 'steelblue', size = 1) +
  geom_vline(aes(xintercept = which(knnAcc == max(knnAcc))), lty = 2) +
  labs(x = 'K', y = 'Accuracy')



##### Exercise 4.11 #####
# (a)
myAuto <- Auto %>% 
  mutate(mpg01 = as.numeric(mpg >= median(mpg))) %>% select(-name, -mpg)

# (b)
# scatterplot
gather(myAuto, var, value, -mpg01) %>% 
  mutate(var = factor(var, levels = names(myAuto)[-8])) %>%
  ggplot() + theme_bw() + theme(legend.position = '') +
  geom_point(aes(x = factor(mpg01), y = value), color = 'gray50') + 
  facet_wrap(~var, ncol = 4, scales = 'free') +
  labs(x = 'mpg01')

# origin은 categorical이라서 scatterplot X

# boxplot
gather(myAuto, var, value, -mpg01) %>% 
  mutate(var = factor(var, levels = names(myAuto)[-8])) %>%
  filter(var != 'origin') %>%
  ggplot() + theme_bw() + theme(legend.position = '') +
  geom_boxplot(aes(x = factor(mpg01), y = value, group = mpg01)) + 
  facet_wrap(~var, ncol = 3, scales = 'free') +
  labs(x = 'mpg01')

# contingency
table(myAuto$origin, myAuto$mpg01)

# (c)
myAuto$origin <- as.factor(myAuto$origin)

set.seed(1)
tr <- sample(1:nrow(myAuto), nrow(myAuto)*0.8)
te <- setdiff(1:nrow(myAuto), tr)
train <- myAuto[tr,]
test <- myAuto[te,]

# (d) LDA
(ldaFit <- lda(mpg01 ~ . -origin -year, data = train))
ldaPred <- predict(ldaFit, test)$class
table(ldaPred, test$mpg01)
mean(ldaPred != test$mpg01) 

# (e) QDA
(qdaFit <- qda(mpg01 ~ . -origin -year, data = train))
qdaPred <- predict(qdaFit, test)$class
table(qdaPred, test$mpg01)
mean(qdaPred != test$mpg01)  

# (f) GLM
summary(glmFit <- glm(mpg01 ~ ., data = train, family = binomial))
glmProb <- predict(glmFit, test)
glmPred <- (ifelse(glmProb > 0.5, 1, 0))
table(glmPred, test$mpg01)
mean(glmPred != test$mpg01) 

# (g) KNN
trainX <- data.matrix(train[,-(7:8)])
testX <- data.matrix(test[,-(7:8)])
trainY <- data.matrix(train[,8])

set.seed(10)
knnPred <- knn(trainX, testX, trainY, k = 1)
table(knnPred, test$mpg01)
mean(knnPred != test$mpg01) 

p <- 100
for (i in 2:p) {
  set.seed(10)
  tmp <- knn(trainX, testX, trainY, k = i)
  knnPred <- data.frame(knnPred, tmp)
}
colnames(knnPred) <- paste0('p', 1:p)

knnError <- NULL
for (i in 1:p) knnError <- c(knnError, mean(knnPred[,i] != test$mpg01))
which.min(knnError)
knnError[knnError == min(knnError)] 

table(knnPred[,50], test$mpg01)
mean(knnPred[,50] != test$mpg01) 

ggplot() + theme_bw() + 
  geom_line(aes(1:p, knnError), col = 'steelblue', size = 1) +
  geom_vline(aes(xintercept = which.min(knnError)), lty = 2) +
  labs(x = 'K', y = 'Error')

