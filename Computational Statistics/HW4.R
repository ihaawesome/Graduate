library(dplyr)
library(ggplot2)
library(gridExtra)
library(randomForest)

# 1) GMM
  mygmm <- function(x, E=10^(-10))
  { 
    niter = 0 ; maxiter = 1000 ; error = 1
    
  # initial values
    n = length(x)
    prop = 0.5
    set.seed(0)
    mu = sort(sample(x, 2)) ; mu1 = mu[1] ; mu2 = mu[2]
    sd1 = sd(x) ; sd2 = sd(x)
    
  # E-M
    while (error >= E & niter <= maxiter) 
    { prop0 <- prop ; mu10 <- mu1 ; mu20 <- mu2 ; sd10 <- sd1 ; sd20 <- sd2
      Ey = c() ; y = c()
      for (i in 1:n)
      { Ey[i] <- prop0*dnorm(x[i], mu10, sd10) / ( prop0*dnorm(x[i], mu10, sd10) + (1-prop0)*dnorm(x[i], mu20, sd20) )
        y[i] <- round(Ey[i])
      }
      
      prop <- sum(y)/n
      mu1 <- sum(y*x)/sum(y) ; mu2 <- sum((1-y)*x)/sum(1-y)
      sd1 <- sqrt(sum(y*(x-mu1)^2)/sum(y)) ; sd2 <- sqrt(sum((1-y)*(x-mu2)^2)/sum(1-y))
      
      error <- abs(prop - prop0)
      niter <- niter + 1
      
      print(paste("error =", error, "niter = ", niter, sep = " "))
    }
    
    output <- list(y=y, prop=prop, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2)
    return(output)
  }
  
  set.seed(100)
  x1 <- c(rnorm(300, mean = 1, sd = 0.5), rnorm(700, mean = 3, sd = 0.3)) # (true value) prop=0.3 ; mu1=1 ; sd1=0.5 ; mu2=3 ; sd2=0.3
  ex1 <- mygmm(x1) ; ex1
  
  system.time( for (i in 1:10) mygmm(x1, E=10^-10) )

  ggplot() + theme_light() + labs(x = "x", y = "density", title = "GMM", subtitle = "N(1, 0.5) & N(3, 0.3)") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "") +
    scale_color_manual("", values = c(real = "gray", n1 = "darkgreen", n2 = "orange")) +
    scale_fill_manual("", values = c(real = "gray", n1 = "darkgreen", n2 = "orange")) +
    geom_area(aes(x1[1:300], dnorm(x1[1:300], 1, 0.5), fill = "real", alpha = 0.5, color = "real")) + 
    geom_area(aes(x1[301:1000], dnorm(x1[301:1000], 3, 0.3), fill = "real", alpha = 0.5, color = "real")) +
    geom_line(aes(x1[1:300], dnorm(x1[1:300], ex1$mu1, ex1$sd1), color = "n1"), size = 1) + 
    geom_line(aes(x1[301:1000], dnorm(x1[301:1000], ex1$mu2, ex1$sd2), color = "n2"), size = 1)
  
  true.y = c(rep(1, 300), rep(0, 700))
  table(true.y, ex1$y)

  
# 2) bivariate normal w/ missing values
  data2 <-data.frame(w1 = c(8, 11, 16, 18, 6, 4, 20, 25, 9, 13), 
                     w2 = c(10, 14, 16, 15, 20, 4, 18, 22, NA, NA))
  
  myem2 <- function(E=10^(-10))
  { 
    w1 = data2$w1 ; w2 = data2$w2 ; n = nrow(data2)
      
  # initial values
    T1 = sum(w1) ; T2 = sum(w2[1:8]) 
    T11 = sum(w1*w1) ; T12 = sum(w1[1:8]*w2[1:8]) ; T22 = sum(w2[1:8]*w2[1:8]) 
    
    mu1 = T1/n ; mu2 = T2/n
    sig11 = (T11-T1*T1/n) / n ; sig12 = (T12-T1*T2/n) / n ; sig22 = (T22-T2*T2/n) / n
    rho = sig12 / sqrt(sig11*sig22)
    
    niter = 0
    maxiter = 1000
    error = 1
    
    while (error >= E & niter <= maxiter)
    { T12_0 <- T12 ; T22_0 <- T22
      mu2_0 <- mu2 ; sig12_0 <- sig12 ; sig22_0 <- sig22 ; rho_0 <- rho
    # E-step  
      Ew29 <- mu2 + sig12/sig11*(w1[9]-mu1)
      Ew29sq <- Ew29^2 + sig22*(1-rho^2)
      Ew210 <- mu2 + sig12/sig11*(w1[10]-mu1)
      Ew210sq <- Ew210^2 + sig22*(1-rho^2)
      
      T2 <- sum(w2[1:8]) + Ew29 + Ew210
      T12 <- sum(w1[1:8]*w2[1:8]) + w1[9]*Ew29 + w1[10]*Ew210
      T22 <- sum(w2[1:8]*w2[1:8]) + Ew29sq + Ew210sq
      
    # M-step  
      mu2 = T2/n
      sig12 = (T12-T1*T2/n) / n ; sig22 = (T22-T2*T2/n) / n
      rho = sig12 / sqrt(sig11*sig22)
      
      error <- abs(mu2 - mu2_0)
      niter <- niter + 1
      
      print(paste("error =", error, "niter = ", niter, sep = " "))
    }
    
    output = data.frame(Ew29, Ew210, mu1, mu2, sig11, sig12, sig22)
    return(output)
  }

  myem2(E=10^-10)
  system.time( for (i in 1:100) myem2() )
  
# 3) Reg. model w/ missing Y's
  red <- read.csv("D:/STAT/CS/HW4/winequality-red.csv", sep = ";")
  white <- read.csv("D:/STAT/CS/HW4/winequality-white.csv", sep = ";")
  
  step.true <- step(lm(quality ~ 1, red), direction = "both",
                    scope = list(lower=lm(quality ~ 1, red), upper=lm(quality ~ ., red)))
  summary(step.true)
  rf.true <- randomForest(quality ~ ., red) ; rf.true
  
  set.seed(100)
  empty <- sort(sample(nrow(red), nrow(red)*0.1))
  true.y <- red[empty, 12]
  red[empty, 12] <- NA

# stepwise reg.    
  myem3.1 <- function(E=10^(-10))
  {
  # initialize  
    red[empty, 12] <- mean(as.numeric(red[, 12]), na.rm = T)
    
    niter = 0
    maxiter = 1000
    error = 1
    
    while(error >= E & niter <= maxiter)
    { yhat_0 <- red[empty, 12]
    # model fitting
      fit <- step(lm(quality ~ 1, red), direction = "both", trace = 0,
                  scope = list(lower=lm(quality ~ 1, red), upper=lm(quality ~ ., red)))
      yhat <- fitted(fit)[empty]
      red[empty, 12] <- yhat
      
      error <- mean( abs(yhat-yhat_0) / yhat_0 )
      niter <- niter + 1
      
      print(paste("error =", error, "niter = ", niter, sep = " "))
    }
    
    return(list(fit = fit, yhat = yhat))
  }  
    
  ex3.1 <- myem3.1()
  summary(ex3.1$fit)
  mean((ex3.1$yhat - true.y)^2)

 
  red[empty, 12] <- NA
  
  # random forest       
  myem3.2 <- function(E=10^(-10))
  {
  # initial
    red[empty, 12] <- mean(as.numeric(red[, 12]), na.rm = T)

    niter = 0
    maxiter = 1000
    error = 1
    
    while(error >= E & niter <= maxiter)
    { yhat_0 <- red[empty, 12]
 
     # model fitting
      set.seed(100)
      fit <- randomForest(quality ~ ., red)
      yhat <- fit$predicted[empty]
      mse <- min(fit$mse)
      red[empty, 12] <- yhat

      error <- mean( abs(yhat-yhat_0) / (yhat_0) )
      niter <- niter + 1
      
      print(paste("error =", error, "niter = ", niter, sep = " "))
    }
    
    return(list(fit = fit, yhat = yhat))
  }  
  
  ex3.2 <-myem3.2()  
  ex3.2 <- myem3.2(0.007)
  mean((ex3.2$yhat-true.y)^2)
  
  par(mfrow = c(1, 2))
  plot(rf.true, lwd = 2, main = "MSE in complete model")
  abline(h = min(rf.true$mse), lty = 2)
  plot(ex3.2$fit, lwd = 2, main = "MSE in E-M final model")
  abline(h = min(ex3.2$fit$mse), lty = 2)

  importance <- data.frame(variable = rownames(ex3.2$fit$importance), IncNodePurity = ex3.2$fit$importance)
  grf2 <- ggplot(importance, aes(reorder(variable, IncNodePurity), y = IncNodePurity, fill = IncNodePurity)) + 
          theme_test() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "") + 
          scale_fill_gradient(low = "lightblue", high = "darkblue") + 
          coord_flip() + labs(x = "", title = "Variable Importance", subtitle = "in E-M final model") + 
          geom_bar(stat = "identity", position = "dodge") ; grf2
  
  importance <- data.frame(variable = rownames(rf.true$importance), IncNodePurity = rf.true$importance)
  grf1 <- ggplot(importance, aes(reorder(variable, IncNodePurity), y = IncNodePurity, fill = IncNodePurity)) + 
    theme_test() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "") + 
    scale_fill_gradient(low = "lightblue", high = "darkblue") + 
    coord_flip() + labs(x = "", title = "Variable Importance", subtitle = "in complete model") + 
    geom_bar(stat = "identity", position = "dodge") ; grf1
  
  grid.arrange(grf1, grf2, ncol = 2)
  
  
# 4) multinomial with complex cell structure
  data4 <- data.frame(type = c("O", "A", "B", "AB"),
                      freq = c(176, 182, 60, 17))
  
  myem4 <- function(E=10^(-10))
  {
    n = 435
    no = 176 ; na = 182 ; nb = 60 ; nab = 17
    
    # initial values
    naa = 0.5*na ; nao = 0.5*na ; nbb = 0.5*nb ; nbo = 0.5*nb
    p = (naa + 0.5*nao + 0.5*nab) / n ; q = (nbb + 0.5*nbo + 0.5*nab) / n ; r = 1-p-q
    
    niter = 0
    maxiter = 1000
    error = 1
    
    while (error >= E & niter <= maxiter)
    { p_0 <- p ; q_0 <- q
    # E-step  
      naa <- na*(p^2/(p^2+2*p*r)) ; nao <- na*(2*p*r/(p^2+2*p*r))
      nbb <- nb*(q^2/(q^2+2*q*r)) ; nbo <- nb*(2*q*r/(q^2+2*q*r))
    # M-step  
      p = (naa + 0.5*nao + 0.5*nab) / n ; q = (nbb + 0.5*nbo + 0.5*nab) / n 
      
      error <- abs(p - p_0)
      niter <- niter + 1
      
      print(paste("error =", error, "niter = ", niter, sep = " "))
    }
    
    output <- data.frame(p, q, r = 1-p-q)
    return(output)
  }

  myem4()  
  system.time( for (i in 1:100 ) myem4() )
  