setwd('C:/Users/HK/Desktop/GitHub/Graduate/DataMining/XGBoost')
##### XGBoost: Parameter Tuning #####
library(data.table)
library(mlr)

# set variable names
setcol <- c('age', 'workclass', 'fnlwgt', 'education', 'education-num', 
            'marital-status', 'occupation', 'relationship', 'race', 'sex', 
            'capital-gain', 'capital-loss', 'hours-per-week', 'native-country', 'target')
# load data
train <- read.table('adult.DATA', header = F, sep = ',', col.names = setcol,
                    na.strings = c('?'), stringsAsFactors = F)
test <- read.table('adult.TEST', header = F, sep = ',', col.names = setcol, skip = 1,
                    na.strings = c('?'), stringsAsFactors = F)

# convert data.frame to data.table
setDT(train)
setDT(test)

# check missing values
table(is.na(train))
sapply(train, function(x) sum(is.na(x))/length(x))
table(is.na(test))
sapply(test, function(x) sum(is.na(x))/length(x))

# quick data cleaning
# remove extra character from target variable
library(stringr)
test [,target := substr(target, 1, nchar(target)-1)] # (?)

# remove leading whitespaces
char_col <- colnames(train)[sapply(test, is.character)]
for (i in char_col) set(train, j = i, value = str_trim(train[[i]], side = 'left'))
for (i in char_col) set(test, j = i, value = str_trim(test[[i]], side = 'left'))

# set all missing value as 'Missing'
train[is.na(train)] <- 'Missing'
test[is.na(test)] <- 'Missing'

# using one hot encoding
label <- train$target
label_test <- test$target
new_tr <- model.matrix(~ . + 0, data = train[,-c('target'), with = F])
new_ts <- model.matrix(~ . + 0, data = test[,-c('target'), with = F])

# convert factor to numeric
label <- as.numeric(label) - 1  # ?????
label_test <- as.numeric(label_test) - 1 # ?????
