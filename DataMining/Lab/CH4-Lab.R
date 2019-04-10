library(ISLR)

########## Stock Market Data ##########
names(Smarket) 
dim(Smarket)
summary(Smarket)
cor(Smarket[,-9])

attach(Smarket)
plot(Volume) # increasing over time


##### Logistic Regression #####
glm.fits <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume, 
                data = Smarket, family = binomial)
summary(glm.fits)
coef(glm.fits)

glm.probs <- predict(glm.fits, type = 'response')
glm.probs[1:10]

contrasts(Direction) # dummy variable

glm.pred <- rep('Down', 1250)
glm.pred[glm.probs > 0.5] <- 'Up'

table(glm.pred, Direction)
mean(glm.pred == Direction)

train <- (Year < 2005)
Smarket2005 <- Smarket[!train,]
Direction2005 <- Direction[!train]
dim(Smarket2005)

# Using subset data
glm.fits <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume, 
                data = Smarket, family = binomial, subset = train)
glm.probs <- predict(glm.fits, Smarket2005, type = 'response')

glm.pred <- rep('Down', 252)
glm.pred[glm.probs > 0.5] <- 'Up'

table(glm.pred, Direction2005)
mean(glm.pred == Direction2005)
mean(glm.pred != Direction2005)

# Using just Lag1 and Lag2
glm.fits <- glm(Direction ~ Lag1 + Lag2, 
                data = Smarket, family = binomial, subset = train)
glm.probs <- predict(glm.fits, Smarket2005, type = "response")
glm.pred <- rep("Down", 252)
glm.pred[glm.probs > 0.5] <- 'Up'
table(glm.pred, Direction2005)

predict(glm.fits, newdata = data.frame(Lag1 = c(1.2, 1.5), Lag2 = c(1.1, -0.8)), type = "response")


##### LDA #####
library(MASS)

lda.fit <- lda(Direction ~ Lag1 + Lag2, data = Smarket, subset = train)
lda.fit
# prior probabilities
# group means
# coefficients (decision rule)

lda.pred <- predict(lda.fit, Smarket2005)
names(lda.pred)
lda.class <- lda.pred$class

table(lda.class, Direction2005)
mean(lda.class == Direction2005)

sum(lda.pred$posterior[,1] >= 0.5)
sum(lda.pred$posterior[,1] < 0.5)

lda.pred$posterior[1:20,1]
lda.class[1:20]
# 사후확률이 50%가 넘으면 Down으로 분류된다

# change threshord
sum(lda.pred$posterior[,1] > 0.9)


##### QDA #####
qda.fit <- qda(Direction ~ Lag1 + Lag2, data = Smarket, subset = train)
qda.fit
# no coefficients (quadratic)

qda.pred <- predict(qda.fit, Smarket2005)
qda.class <- qda.pred$class

table(qda.class, Direction2005)
mean(qda.class == Direction2005)


##### K-NN #####
library(class)
trainX <- cbind(Lag1, Lag2)[train,]
testX <- cbind(Lag1, Lag2)[!train,]
trainDirection <- Direction[train]

set.seed(1)
knn.pred <- knn(trainX, testX, trainDirection, k = 1) # 1-NN
table(knn.pred, Direction2005)
mean(knn.pred == Direction2005)


knn.pred <- knn(trainX, testX, trainDirection, k = 3) # 3-NN
table(knn.pred, Direction2005)
mean(knn.pred == Direction2005)



########## Caravan Insurance Data ##########
dim(Caravan) # 85 predictors...
attach(Caravan)
summary(Purchase) # y 

##### KNN #####
stdX <- scale(Caravan[,-86]) # standardize
apply(Caravan[,1:2], 2, var)
apply(stdX[,1:2], 2, var)

# Split the dataset
test <- 1:1000
trainX <- stdX[-test,]
testX <- stdX[test,]
trainY <- Purchase[-test]
testY <- Purchase[test]

set.seed(1)
knn.pred <- knn(trainX, testX, trainY, k = 1) # 1-NN
mean(testY != knn.pred) # misclassification
mean(testY != "No") # misclassification

table(knn.pred, testY)

knn.pred <- knn(trainX, testX, trainY, k = 3)
table(knn.pred, testY)

knn.pred <- knn(trainX, testX, trainY, k = 5)
table(knn.pred, testY)


##### Logistic Regression #####
glm.fits <- glm(Purchase ~ ., data = Caravan, family = binomial, subset = -test)
glm.probs <- predict(glm.fits, Caravan[test,], type = "response")
glm.pred <- rep("No", 1000)
glm.pred[glm.probs > 0.5] <- 'Yes'
table(glm.pred, testY) # 꽝

glm.pred <- rep("No", 1000)
glm.pred[glm.probs > 0.25] <- 'Yes' # threshold 수정
table(glm.pred, testY) 



