library(ISLR)
library(MASS)
library(car)
library(tidyverse)


# Data Mining HW3 ----------
# CH3: Exercise 8, 9, 13, 14
# --------------------------


# plot options
  plotlm <- function(model) {
    par(mfrow = c(2, 2))
    plot(model, pch = 19, col = 'gray', lwd = 2)
  }
  myplotlm <- function(model) {
    par(mfrow = c(1,2))
    plot(model, pch = 19, col = 'gray', lwd = 2, which = c(1,5))
  }
  

# Exercise 8 (Simple LR)
  attach(Auto)
# (a)
  summary(lm1a <- lm(mpg ~ horsepower))
  cor(mpg, horsepower)
  coef(lm1a)[2] 
  
  predict(lm1a, newdata = data.frame(horsepower = 98), interval = 'confidence')
  predict(lm1a, newdata = data.frame(horsepower = 98), interval = 'prediction')
  
# (b)  
  plot(horsepower, mpg, pch = 19, col = 'gray')
  abline(lm1a, col = 'steelblue', lwd = 2)  

# (c)
# plot all  
  plotlm(lm1a)
  
# residuals vs fitted
  plot(predict(lm1a), residuals(lm1a), col = 'gray', lwd = 2)
  plot(predict(lm1a), rstudent(lm1a), col = 'gray', lwd = 2)
# leverage
  plot(hatvalues(lm1a), col = 'gray', lwd = 2)
  points(which.max(hatvalues(lm1a)), max(hatvalues(lm1a)), col = 'red', pch = 16)


# Exercise 9 (Multiple LR)
# (a)
  pairs(Auto[,-9], col = 'lightblue') 

# (b)
  cor(Auto[,-9])

# (c)
  summary(lm2c <- lm(mpg ~ . -name, Auto))
  coef(lm2c)['year'] 
  
  vif(lm2c)
  
# (d)
  plotlm(lm2c)
  
  plot(predict(lm2c), residuals(lm2c), col = 'gray', lwd = 2)
  plot(predict(lm2c), rstudent(lm2c), col = 'gray', lwd = 2)
  plot(hatvalues(lm2c), col = 'gray', lwd = 2)
  points(which.max(hatvalues(lm2c)), max(hatvalues(lm2c)), col = 'red', pch = 16) 

# (e)
  summary(lm2e <- lm(mpg ~ .*., Auto[,-9]))
  coef2e <- as.data.frame(summary(lm2e)$coefficients) %>% round(4)

# (f)  
  summary(lm2f <- lm(mpg ~ cylinders + poly(displacement, 2) +
                       log(horsepower) + log(weight) + poly(acceleration, 2) + year + origin))

  
# Exercise 13
  n <- 100
# (a)
  set.seed(1)
  x <- rnorm(n)
# (b)
  eps <- rnorm(n, 0, 0.5)
# (c)
  beta0 <- -1 ; beta1 <- 0.5
  y <- beta0 + beta1*x + eps
  length(y)  

# (d)
  par(mfrow = c(1, 1))
  plot(x, y, col = 'gray', pch = 19, main = 'Scatterplot (y, x)')

# (e)
  summary(lm3e <- lm(y ~ x))
  plotlm(lm3e)
  
# (f)
  par(mfrow = c(1, 1))
  plot(x, y, col = 'gray', pch = 19)
  abline(-1, 0.5, col = 'steelblue', lwd = 2)
  abline(lm3e, col = 'red', lwd = 2, lty = 2)
  legend('topleft', legend = c('population', 'fitted'), col = c('steelblue', 'red'), 
         lty = c(1, 2), lwd = 2, bty = 'n')

# (g)
  summary(lm3g <- lm(y ~ poly(x, 2)))
  plotlm(lm3g)
  
  anova(lm3e, lm3g, test = 'F')
  
  par(mfrow = c(1, 1))
  plot(x, y, col = 'gray', pch = 19)
  abline(-1, 0.5, col = 'steelblue', lwd = 2)
  abline(x = sort(x), y = sort(fitted(lm3g)), col = 'red', lwd = 2, lty = 2)
  legend('topleft', legend = c('population', 'fitted'), col = c('steelblue', 'red'), 
         lty = c(1, 2), bty = 'n', lwd = 2)
  
  
# (h)  
  eps_less <- rnorm(n, 0, 0.1)
  y_less <- beta0 + beta1*x + eps_less

  summary(lm3h <- lm(y_less ~ x))  
  
  plotlm(lm3h)
  
  par(mfrow = c(1, 1))
  plot(x, y_less, col = 'gray', pch = 19)
  abline(-1, 0.5, col = 'steelblue', lwd = 2)
  abline(lm3h, col = 'red', lwd = 2, lty = 2)
  legend('topleft', legend = c('population', 'fitted'), col = c('steelblue', 'red'), 
         lty = c(1, 2), lwd = 2, bty = 'n')
  
# (i)
  eps_more <- rnorm(n, 0, 1)
  y_more <- beta0 + beta1*x + eps_more  
  
  summary(lm3i <- lm(y_more ~ x))
  
  plotlm(lm3i) 
  
  par(mfrow = c(1, 1))
  plot(x, y_more, col = 'gray', pch = 19)
  abline(-1, 0.5, col = 'steelblue', lwd = 2)
  abline(lm3i, col = 'red', lwd = 2, lty = 2)
  legend('topleft', legend = c('population', 'fitted'), col = c('steelblue', 'red'), 
         lty = c(1, 2), lwd = 2, bty = 'n')

# (j)
  rbind(
    cbind(coef(lm3e), confint(lm3e, interval = 'confidence')), # original
    cbind(coef(lm3h), confint(lm3h, interval = 'confidence')), # less
    cbind(coef(lm3i), confint(lm3i, interval = 'confidence'))  # more
  ) %>% round(4)
  

# Exercise 14: collinearity
# (a)  
  set.seed(1)
  x1 <- runif(100)
  x2 <- 0.5*x1 + rnorm(100)/10

  beta <- c(2, 2, 0.3)
  y <- beta[1] + beta[2]*x1 + beta[3]*x2 + rnorm(100)
  
# (b)  
  cor(x1, x2)
  plot(x1, x2, col = 'gray', pch = 19)
  abline(lm(x2 ~ x1), col = 'pink', lwd = 2)

# (c)    
  summary(lm4c <- lm(y ~ x1 + x2))
  rbind(beta, coef(lm4c))
  confint(lm4c, interval = 'confidence') %>% cbind(estimate = coef(lm4c))
  
  vif(lm4c)

# (d)
  summary(lm4d <- lm(y ~ x1)) 

# (e)    
  summary(lm4e <- lm(y ~ x2)) 
  
# (f)
  anova(lm4c)
  
# (g)
  data4g <- rbind(data.frame(y, x1, x2), c(6, 0.1, 0.8))
  
  summary(lm4g12 <- lm(y ~ x1 + x2, data4g)) ; myplotlm(lm4g12)
  summary(lm4g1 <- lm(y ~ x1, data4g)) ; myplotlm(lm4g1)  
  summary(lm4g2 <- lm(y ~ x2, data4g)) ; myplotlm(lm4g2)    
  
  gdat <- gather(data4g, var, x, -y)
  ggplot(gdat, aes(x, y)) + theme_bw() +
    geom_point(col = 'gray') + geom_smooth(method = 'lm', se = F, col = 'pink', size = 2) + 
    geom_point(data = gdat[c(101,202),], aes(x, y), col = 'red', size = 2) +
    facet_grid(~var, scales = 'free')
  