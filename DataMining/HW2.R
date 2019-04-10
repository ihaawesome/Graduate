setwd('C:/Users/HK/Desktop/2019-1/DataMining/HW2')
library(ISLR)
library(MASS)
library(tidyverse)

# --------------------------------------------------------------------------------------
# Exercise 2.9

# (a)
  str(Auto)

# (b)
  quan <- select(Auto, -mpg, -origin, -name)
  qual <- select(Auto, origin, name)
  
  t(apply(quan, 2, range))
  
# (c)
  cbind(mean = apply(quan, 2, mean), sd = apply(quan, 2, sd)) %>% round(2)

# (d)
  quan.s <- quan[-c(10:85),] ; qual.s <- qual[-c(10:85),]
  t(apply(quan.s, 2, range))  
  cbind(mean = apply(quan.s, 2, mean), sd = apply(quan.s, 2, sd)) %>% round(2)
  
# (e)-(f)
  g <- ggplot(Auto) + theme_light()
  g + geom_histogram(aes(displacement), col='white')
  g + geom_bar(aes(cylinders), col='white')
  g + geom_histogram(aes(horsepower), col='white', binwidth=10)
  g + geom_histogram(aes(weight), col='white', binwidth=100)
  g + geom_histogram(aes(acceleration), col='white')
  g + geom_histogram(aes(year))
  g + geom_bar(aes(origin))
  
  g + geom_point(aes(displacement, mpg))
  g + geom_boxplot(aes(x=cylinders, y=mpg, group=cylinders))
  g + geom_point(aes(horsepower, mpg))
  g + geom_point(aes(weight, mpg))
  g + geom_point(aes(acceleration, mpg))
  g + geom_boxplot(aes(x=year, y=mpg,group=year))
  g + geom_point(aes(year, mpg))
  g + geom_boxplot(aes(x=origin, y=mpg, group=origin))
  g + geom_histogram(aes(mpg))
  
  panel.cor <- function(x, y){
    usr <- par('usr') 
    on.exit(par(usr))
    par(usr=c(0,1,0,1))
    r <- round(cor(x,y), digits=4)
    text(0.5, 0.5, r, cex=1.2)
  }
  
  pairs(quan, upper.panel=function(x, y) points(x, y, pch=20), lower.panel=panel.cor) 
  
  Auto %>% 
    select(-mpg,-name) %>% 
    gather(var, value, -origin) %>% 
    mutate(var=factor(var, levels=c('cylinders','displacement','horsepower','weight','acceleration','year'))) %>% 
    ggplot() + theme_test() +
    geom_boxplot(aes(x=origin, y=value, group=origin)) + 
    facet_wrap(~var, scales='free', ncol=3)
  
  Auto %>% 
    select(-name) %>% 
    gather(var, value, -mpg) %>% 
    mutate(var=factor(var, levels=colnames(Auto)[2:8])) %>% 
    ggplot() + theme_test() +
    geom_point(aes(x=value, y=mpg)) + 
    facet_wrap(~var, scales='free', ncol=4)
  
  cormat1 <- cor(Auto[,1:7])
  ctmp <- cormat1[1,] %>% round(4)
  par(mfrow=c(2,4))
  for (i in 2:8) {
    if (i==8) boxplot(mpg~origin, Auto, xlab='origin', pch=20, cex.lab=1.5)
    else plot(Auto[,i], Auto$mpg, xlab=colnames(Auto)[i], ylab='', main=ctmp[i], pch=20, cex.lab=1.5)
  }
  
# --------------------------------------------------------------------------------------  
# Exercise 2.10
# (a) 
  str(Boston)
  dim(Boston)
  
# (b)
  png('Pairwise Scatterplot of Predictors.png', width=1000, height=800)
  pairs(select(Boston, -medv, -chas))
  dev.off()
  
  tmp <- select(Boston, -chas, -medv)
  colnames(tmp)
  
  png('predictors-chas.png', width=1000, height=400)
  par(mfrow=c(2,6))
  for(i in 1:12) {
  boxplot(tmp[,i]~chas, Boston,
          main=colnames(tmp)[i], pch=20, cex.main=2)
  }
  dev.off()
  
  cormat <- cor(Boston)
  ctmp <- cormat[14,-14] %>% round(4)
  
  png('predictors-medv.png', width=800, height=600)
  par(mfrow=c(3,5))
  for (i in 1:13) {
    if(i==4) {
      boxplot(medv~chas, Boston, main=ctmp[i], cex.main=1.5,
              xlab='chas', ylab='medv', pch=20, cex.lab=1.5)
    } else {
      plot(Boston[,i], Boston$medv, pch=20, main=ctmp[i], cex.main=1.5,
           xlab=colnames(Boston)[i], ylab='medv', cex.lab=1.5)
    } 
  }
  dev.off()
  
# (c)
  ctmp <- cormat[1,-1] %>% round(4)
  png('predictors-crim.png', width=800, height=600)
  par(mfrow = c(3, 4))
  for (i in 2:13) {
    if(i==4) {
      boxplot(crim~chas, Boston, main=ctmp[i-1], cex.main=1.5,
              xlab='chas', ylab='medv', pch=19, cex.lab=1.5)
    } else {
      plot(Boston[,i], Boston$crim, pch=19, main=ctmp[i-1], cex.main=1.5, 
           xlab=colnames(Boston)[i], ylab='crim', cex.lab=1.5)
    } 
  }
  dev.off()
  
# or
  tmp <- gather(select(Boston, -medv, -chas), var, value, -crim)
  tmp$var <- factor(tmp$var, levels=colnames(Boston)[-1])
  ggplot(tmp) + theme_test() +
    geom_point(aes(value, crim)) + 
    facet_wrap(~var, scales='free', ncol=4)
  
  ggplot(Boston) + theme_test() +
    geom_boxplot(aes(x = factor(chas), y = medv, group = chas)) + labs(x = 'chas') 
  
  cor(select(Boston, -medv))[1,-1] %>% round(4)

# (d)
  g + geom_histogram(aes(crim), binwidth=1, fill='lightblue', color='black')
  g + geom_histogram(aes(tax), binwidth=10, fill='lightblue', color='black')
  g + geom_histogram(aes(ptratio), binwidth=0.1, fill='lightblue', color='black')
  
  g + geom_point(aes(1:nrow(Boston), crim)) + labs(x='row.index')
  g + geom_point(aes(1:nrow(Boston), tax)) + labs(x='row.index')
  g + geom_point(aes(1:nrow(Boston), ptratio)) + labs(x='row.index')
  
  png('index-crim.png', 300, 300) ; plot(Boston$crim, pch=20, main='crim', ylab='') ; dev.off()
  png('index-tax.png', 300, 300) ; plot(Boston$tax, pch=20, main='tax', ylab='') ; dev.off()
  png('index-ptratio.png', 300, 300) ; plot(Boston$ptratio, pch=20, main='ptratio', ylab='') ; dev.off()
  
  t(apply(select(Boston, crim, tax, ptratio), 2, range))   
  
# (e)
  table(Boston$chas)
  
# (f)
  median(Boston$ptratio)
  g + geom_histogram(aes(ptratio), fill = 'gray', binwidth = 0.5) +
    geom_vline(aes(xintercept = median(ptratio)), color = 'red', linetype = 2, size = 1)
  
  hist(Boston$ptratio, xlab='ptratio', breaks=50, col='gray', main='')
  abline(v=median(Boston$ptratio), lty=2)
  
# (g)
  filter(Boston, medv==min(medv))
  apply(Boston, 2, summary) %>%round(4) 
  which(Boston$medv==min(Boston$medv))
  
# (h) 
  nrow(filter(Boston, rm > 7))
  nrow(filter(Boston, rm > 8))
  nrow(filter(Boston, rm <= 8))
  
  rbind(apply(Boston,2, summary),
        apply(filter(Boston, rm > 8), 2, mean)) %>% round(2)
  
  g <- ggplot(Boston) + theme_bw()
  g + geom_histogram(aes(rm), binwidth = 0.1, fill = 'gray') +
    geom_vline(aes(xintercept = 7), linetype = 2, color = 'red', size = 1) +
    geom_vline(aes(xintercept = 8), linetype = 2, color = 'red', size = 1)

  