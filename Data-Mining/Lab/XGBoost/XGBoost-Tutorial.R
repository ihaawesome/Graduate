setwd('C:/Users/HK/Desktop/GitHub/Graduate/DataMining/XGBoost')
##### XGBoost Tutorial #####

# how to use Xgboost to build a model and make predictions
# gradient boosting framework: linear & tree learning

# Input Type: matrix, dgCMatrix, xgb.DMatrix (recommended)

# 1.2 Installation
library(xgboost)

# 1.3 Learning
# 1.3.2 Dataset loading
# use Mushroom data
data(agaricus.train, package = 'xgboost')
data(agaricus.test, package = 'xgboost')
train <- agaricus.train
test <- agaricus.test

str(train) # data (X ; dgCMatrix class), label (y)
dim(train$data) # 80%
dim(test$data)  # 20%

# 1.3.3 Basic training
# max_depth = depth of the trees
# nthread = the number of cpu threads to use
# nrounds
# = each round enhances the model by further reducing the difference 
#   between ground truth and prediction

# dgCMatrix class
bstSparse <- xgboost(
  data = train$data, label = train$label, 
  max_depth = 2, eta = 1, nthread = 2, nrounds = 2, 
  objective = 'binary:logistic'
)
# xgb.DMatrix class
dtrain <- xgb.DMatrix(data = train$data, label = train$label)
bst <- xgboost(
  data = dtrain, max_depth = 2, eta = 1, nthread = 2, nrounds = 2, 
  objective = 'binary:logistic', verbose = 2 # 0 (silence)
)

# 1.5 Perform the prediction
pred <- predict(bst, newdata = test$data)
print(length(pred))
print(head(pred)) # predicted probabilities

# 1.6 Transform the regression in a binary classification
# The only thing that XGBoost does is a regression.
# set the rule that if a specific observation is classified as 1.
prediction <- as.numeric(pred > 0.5)
print(head(prediction)) # predicted label

# 1.7 Measuring model performance
err <- mean(prediction != test$label) # misclassification rate
print(paste('test-error =', err))

# 1.8 Advanced features
dtrain <- xgb.DMatrix(data = train$data, label = train$label)
dtest <- xgb.DMatrix(data = test$data, label = test$label)

# 1.8.2 Measure learning progress with 'xgb.train'
# follow the progress of the learning after each round to evalutate an overfitting.
# cross-validation
watchlist <- list(train = dtrain, test = dtest)
bst <- xgb.train(
  data = dtrain, max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
  watchlist = watchlist, objective = 'binary:logistic'
)

# have some evaluation metrics
bst <- xgb.train(
  data = dtrain, max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
  watchlist = watchlist, eval_metric = 'error', eval_metric = 'logloss',
  objective = 'binary:logistic'
)

# 1.8.3 Linear boosting
# All the learnings we have performed were based on boosting trees.
# Second algorithm: linear boosting
bst <- xgb.train(
  data = dtrain, booster = 'gblinear', 
  max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
  watchlist = watchlist, eval_metric = 'error', eval_metric = 'logloss',
  objective = 'binary:logistic'
)
# to catch a linear link, liniear boosting is the best.
# to catch a non=linear link, decision trees can be much better.

# 1.8.4 Manipulating xgb.DMatrix
xgb.DMatrix.save(dtrain, 'dtrain.buffer') # save
dtrain2 <- xgb.DMatrix('dtrain.buffer')   # load

bst <- xgb.train(
  data = dtrain2, max_depth = 2, eta = 1, nthread = 2, nrounds = 2, 
  watchlist = watchlist, objective = "binary:logistic"
)

# 1.8.4.2 Information Extraction
label <- getinfo(dtest, 'label')
pred <- predict(bst, dtest)
err <- as.numeric(sum(as.integer(pred > 0.5) != label)) / length(label)
print(paste('test-error =', err))

# 1.8.5 View feature importance/influence
importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix)

# 1.8.5.1 View the trees
xgb.dump(bst, with_stats = TRUE)
xgb.plot.tree(model = bst)

# 1.8.5.2 Save and load models
xgb.save(bst, 'xgboost.model') # save to a local MODEL file
bst2 <- xgb.load('xgboost.model')
pred2 <- predict(bst2, test$data)
print(sum(abs(pred2 - pred))) # same

rawVec <- xgb.save.raw(bst) # save model to R's raw vector
bst3 <- xgb.load(rawVec)
pred3 <- predict(bst3, test$data)
print(sum(abs(pred3 - pred)))
