library(dplyr)
library(ggplot2)
library(GenSA)
library(GA) 

# (1) kmeans
  set.seed(1234)
  x1 = matrix(rnorm(100, sd = 0.5), ncol = 2)
  x2 = matrix(rnorm(100, 1, sd = 0.5), ncol = 2)
  x = rbind(x1, x2) ; colnames(x) = c("x1", "x2") ; x = as.data.frame(x)

  kmean = kmeans(x, 2, nstart = 10) ; kmean

  ggplot() + theme_test() + theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + 
    geom_point(data = data.frame(x, cluster = factor(kmean$cluster)), aes(x1, x2, color = cluster), size = 2) + scale_color_viridis_d() +
    geom_point(data = data.frame(kmean$centers), aes(x1, x2), size = 3, shape = 21, color = "black", fill = "white", stroke = 3) +
    labs(title = "Result of R kmeans function")


# function
  mykmeans <- function(x, k, maxiter=10, nstart=1) {
    
    n = nrow(x)
    tss = sum(dist(x)^2) / n
    niter = 0
    
    init.mean = list() ; output = list()
    for (ns in 1:nstart) {
  
      init.mean[[ns]] <- x[sample(n,k),]
    
      cl.mean = init.mean[[ns]]
      cl = numeric(k)
      ss = matrix(numeric(n*k), nrow = n)
      change = F
      
      while (change == F & niter <= maxiter) {
        cl.mean0 <- cl.mean
        cl0 <- cl
  
        for (j in 1:k) {
          for (i in 1:n) { 
            ss[i,j] = sum((x[i,]-cl.mean[j,])^2)
            cl[i] = which.min(ss[i,])
          } 
        }
        spl = split(x, cl)
        for (j in 1:k) { cl.mean[j,] = apply(spl[[j]], mean, MARGIN = 2) }
        
        change = identical(cl, cl0)
        niter <- niter + 1
      }
      
      size = c() ; wss = c()
      for (j in 1:k) { 
        size[j] = nrow(spl[[j]])
        wss[j] = sum(dist(spl[[j]])^2) / size[j] 
      }
      tot.wss = sum(wss)
      bss = tss - tot.wss
      
      output[[ns]] <- list(cluster=cl, centers=cl.mean, tss=tss, wss=wss, tot.wss=tot.wss, bss=bss, size=size, niter=niter)
    }
    
    final <- output[[which.min(tot.wss)]]
    return(final)
  }
  
  mykmean = mykmeans(x, k=2, nstart=10) ; mykmean
  
  data = data.frame(x, cluster = factor(mykmean$cluster))
  ggplot() + theme_test() + 
    geom_point(data = data, aes(x1, x2, color = cluster), size = 2) + scale_color_viridis_d() +
    geom_point(data = mykmean$centers, aes(x1, x2), size = 3, shape = 21, color = "black", fill = "white", stroke = 3) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + labs(title = "Result of mykmeans function") 
  
  
# (2) baseball
  baseball = read.csv("C:/Users/HG/Documents/R/data/baseball.txt", header = T, sep = " ")

  myaic <- function(par) {
    index = as.logical(par >= 0.5)
    y = baseball[1]
    x = baseball[-1][index]
    
    model <- lm(log(salary) ~., data.frame(y, x))
    aic <- extractAIC(model)[2]
    
    return(aic)
  }
  
  
# (2-1) GenSA
  gensa.out <- GenSA(fn = myaic, 
                     lower = rep(0, 27), upper = rep(1, 27), 
                     control = list(max.time = 10, seed = 1))
  gensa.out$value
  as.logical(gensa.out$par >= 0.5) %>% which()
  
  gensa.info <- as.data.frame(gensa.out$trace.mat)
  
  g1 <- ggplot(gensa.info) + 
    geom_line(aes(1:length(function.value), current.minimum), color = "blue", linetype = "dashed") + 
    geom_line(aes(1:length(function.value), function.value)) + 
    labs(x = "iterations", y = "AIC", title = "GenSA Output") + ylim(-420, -360) + theme_test() +
    theme(plot.title = element_text(hjust = 0.5))
  
  g2 <- ggplot(gensa.info) + geom_line(aes(1:length(temperature), temperature), linetype = "dashed") + 
    labs(x = "iterations", title = "")  + theme_test() +
    theme(plot.subtitle = element_text(hjust = 0.5))
  
  library(gridExtra)
  grid.arrange(g1, g2, nrow = 2)
  
  
  
# (2-2) GA  
  ga.out <- ga(type = "real-valued",
               fitness = function(x) -myaic(x),
               lower = rep(0, 27), upper = rep(1, 27),
               popSize = 20, maxiter = 100, seed = 1)
  
  ga.info <- summary(ga.out)
  ga.info$fitness
  which(ga.info$solution[1,] >= 0.5)
  
  plot(ga.out, main = "GA Output")


  
  
  
  
  
  
  
    
# try ---------------------------------------------
  
  myaic <- function(par) {
    
    y = baseball[1]
    x = sample(baseball[-1], as.integer(par))
    current <- lm(log(salary) ~ ., data.frame(y, x))
    
    new = sample(baseball[-1], 1)
    
    if (!(names(new) %in% names(x))) {
      
      newdata = data.frame(y, x, new)
      add <- lm(log(salary) ~ ., newdata)
      
      if (extractAIC(current)[2] > extractAIC(add)[2]) { current <- add }
    }
    else if (names(new) %in% names(x)) {
      
      x = x[names(x) != names(new)]
      newdata = data.frame(y, x)
      subst <- lm(log(salary) ~ ., newdata)
      
      if (extractAIC(current)[2] > extractAIC(subst)[2]) { current <- subst }
    }
    
    myaic = extractAIC(current)[2]
    return(myaic)
  }

  
  gensa.out <- GenSA(fn = myaic, 
                     lower = 1, upper = 26, 
                     control = list(maxit = 1000, seed = 123))
  
  gensa.info = as.data.frame(gensa.out$trace.mat)
  ggplot(gensa.info) + geom_line(aes(1:length(function.value), function.value)) + 
    labs(x = "iterations", y = "AIC", title = "GenSA Output") +  ylim(-420,-390) + theme_test()
  ggplot(gensa.info) + geom_line(aes(1:length(temperature), temperature), linetype = "dashed") + 
    labs(x = "index", title = "Temperature") + ylim(0,1000) + theme_test()
  
  plot(gensa.info$function.value, ylim = c(-420,-380), type = "l",
       ylab = "AIC", main = "GenSA Output")
  

  
  ga.out <- ga(type = "real-valued",
               fitness = function(x) -myaic(x),
               lower = 1, upper = 26,
               popSize = 20, maxiter = 300, seed = 123) # keepBest = T, optim = T
  ga.info = summary(ga.out)
  plot(ga.out, main = "GA Output")
  
      