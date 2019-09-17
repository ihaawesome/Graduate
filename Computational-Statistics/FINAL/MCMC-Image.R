library(dplyr)
library(ggplot2)
setwd("C:/Users/HK/Desktop/18-2/CS/FINAL")

############################## Image Restoration ##############################
  
  myimage <- function(value, main = "") {
    image(axes = F, matrix(value, ncol = n2), main = main, 
          col = grey(seq(1, 0, length = length(value))))
  }
  
  berry <- read.table("DATA/serviceberry.dat", header = T)$x
  n1 = 54 ; n2 = 46
  
  img <- berry
  y <- as.vector(img)
  
  x <- y
  p.swit = 0.5
  set.seed(1) ; swit <- sample(length(y), size = round(p.swit*length(y)))
  x[swit] <- ifelse(x[swit]==1, 0, 1) # make degraded image
  
  par(mfrow = c(1, 2))
  myimage(y, "TRUE")
  myimage(x, "DEGRADED")
  
  sum(y==x) / length(y)
  
  x <- rep(0, length(y))
  
  mynb2 <- function(i) {
    if (i%%n1==1) {
      idx <- c(ifelse(i<=n1, NA, i-n1),        ifelse(i>n1*(n2-1), NA, i+n1),
               ifelse(i<=n1, NA, i-n1+1), i+1, ifelse(i>n1*(n2-1), NA, i+n1+1))
    }
    else if (i%%n1==0) {
      idx <- c(ifelse(i<=n1, NA, i-n1-1), i-1, ifelse(i>n1*(n2-1), NA, i+n1-1),
               ifelse(i<=n1, NA, i-n1),        ifelse(i>n1*(n2-1), NA, i+n1))
    }
    else {
      idx <- c(ifelse(i<=n1, NA, i-n1-1), i-1, ifelse(i>n1*(n2-1), NA, i+n1-1),
               ifelse(i<=n1, NA, i-n1),        ifelse(i>n1*(n2-1), NA, i+n1),
               ifelse(i<=n1, NA, i-n1+1), i+1, ifelse(i>n1*(n2-1), NA, i+n1+1))
    }
    return(sort(idx[!is.na(idx)]))  
  }
  
  # 1st order neighborhood
  mynb1 <- function(i) {
    if (i%%n1==1) {      idx <- c(ifelse(i<=n1, NA, i-n1), i+1, ifelse(i>n1*(n2-1), NA, i+n1)) }
    else if (i%%n1==0) { idx <- c(ifelse(i<=n1, NA, i-n1), i-1, ifelse(i>n1*(n2-1), NA, i+n1)) }
    else {               idx <- c(ifelse(i<=n1, NA, i-n1), i-1, i+1, ifelse(i>n1*(n2-1), NA, i+n1)) }
    return(sort(idx[!is.na(idx)]))  
  }
  
  mynb <- mynb2
  
  mychain <- function(niter, alpha=1, beta=0.8) {
    out <- matrix(nrow = length(y), ncol = niter)
    out[,1] <- x
    for (t in 2:niter) {
      xt <- out[,(t-1)]
      for (i in 1:length(y)) {
        p <- 1 / (1 + exp(alpha*(((y[i]==0)-(y[i]==1))) + 
                            beta*sum((xt[mynb(i)]==0)-(xt[mynb(i)]==1))))
        xt[i] <- rbinom(1, 1, p)
      }
      out[,t] <- xt
      
      if (t%%10 == 0) print(paste("iter =", t))
    }
    return(out)
  }
  
  niter = 1000
  try1 <- mychain(niter)
  
  par(mfrow = c(2, 2))
  myimage(y, "TRUE")
  myimage(rowMeans(try1), "MEAN")
  myimage(try1[,2], "ITER 1")
  myimage(try1[,niter], paste("ITER", niter))
  
  sum(y==try1[,1]) / length(y)
  sum(y==try1[,niter]) / length(y)
