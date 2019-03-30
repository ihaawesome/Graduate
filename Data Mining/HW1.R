library(MASS)
library(tidyverse)
library(scatterplot3d)
library(rgl)

# ------------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------------

# generate p-dim random points 
  mysample <- function(n, p) {
    s <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
    d <- runif(n, 0, 1)
    S <- apply(s, 1, function(x) sum(x^2))
    
    s <- d^(1/p) * s / sqrt(S)
    return(s)
  }
  
  distance <- function(sample) {
    apply(sample, 1, function(x) sqrt(sum(x^2)))
  }
  
# distance to the closest point in the sample
  closest <- function(sample) {
    d <- distance(sample)
    return(c(m = min(d), sample[which.min(d),]))
  }
  
# true median distance
  closest_m <- function(n, p) (1-(0.5)^(1/n))^(1/p)

  n <- c(100, 1000, 10000)
  p <- c(2, 3, 5, 10)
  
  d_true <- matrix(nrow = 3, ncol = 4)
  rownames(d_true) <- n ; colnames(d_true) <- p
  for (i in 1:3) for (j in 1:4) d_true[i, j] <- closest_m(n[i], p[j])
  

# ------------------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------------------

# surface & inside (p = 2)
  surface2d <- mvrnorm(100, mu = rep(0, 2), Sigma = diag(2))
  S <- apply(surface2d, 1, function(x) sum(x^2))
  surface2d <- surface2d / sqrt(S)
  
  inside2d <- mvrnorm(100, mu = rep(0, 2), Sigma = diag(2))
  S <- apply(inside2d, 1, function(x) sum(x^2))
  inside2d <- inside2d / sqrt(S) *runif(100)
  
  
  plot(surface2d[,1], surface2d[,2], pch = 19, col = 'gray',
       xlab = 'x1', ylab = 'x2',  main = 'Sample on(inside) the Surface\np = 2')
  points(inside2d[,1], inside2d[,2], pch = 19)
  
  
# scatterplot (p = 2)  
  sample2d <- mysample(1000, 2) ; colnames(sample2d) <- paste0('x', 1:2)
  mininfo <- closest(sample2d)
  mind <- round(mininfo[1], 5) 
  minpoint <- mininfo[-1]
  
  par(mfrow = c(1, 1))
  plot(sample2d[,1], sample2d[,2], pch = 19, col = 'gray',
       xlim = c(-1, 1), ylim = c(-1, 1),  
       xlab = 'x1', ylab = 'x2', main = 'p = 2', sub = paste('Min Distance =', mind))
  points(0, 0, pch = 19)
  points(minpoint[1], minpoint[2], col = 'red', pch = 19)
  
  ggplot(as.data.frame(sample2d)) + theme_light() + 
    geom_point(aes(0, 0), size = 2) +
    geom_point(aes(x1, x2), color = 'gray', size = 2) +
    geom_point(aes(minpoint[1], minpoint[2]), color = 'red', size = 2) +
    labs(x = 'x1', y = 'x1', title = 'p = 2', subtitle = paste('Min Distance =', mind)) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

# scatterplot (p = 3)
  sample3d <- mysample(1000, 3) ; colnames(sample3d) <- paste0('x', 1:3)
  mininfo <- closest(sample3d)
  mind <- round(mininfo[1], 5) 
  minpoint <- mininfo[-1]
  
  plt3d <- scatterplot3d(
    sample3d[,1], sample3d[,2], sample3d[,3], color = 'gray', pch = 19, 
    xlab = 'x1', ylab = 'x2', zlab = 'x3',    
    main = 'p = 3', sub = paste('Min Distance =', mind)
  )
  plt3d$points3d(0, 0, 0, pch = 19)
  plt3d$points3d(minpoint[1], minpoint[2], minpoint[3], pch = 19, col = 'red')

  
# Histograms
  ggplot(data.frame(d = distance(mysample(10000, 2)))) + theme_light() +
    geom_histogram(aes(d), color = 'white') + labs(x = 'distance', title = 'p = 2') +
    theme(plot.title = element_text(hjust = 0.5))

  ggplot(data.frame(d = distance(mysample(10000, 3)))) + theme_light() +
    geom_histogram(aes(d), color = 'white', binwith = 0.01) + labs(x = 'distance', title = 'p = 3') +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot(data.frame(d = distance(mysample(10000, 5)))) + theme_light() +
    geom_histogram(aes(d), color = 'white') + labs(x = 'distance', title = 'p = 5') +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot(data.frame(d = distance(mysample(10000, 10)))) + theme_light() +
    geom_histogram(aes(d), color = 'white') + labs(x = 'distance', title = 'p = 10') +
    theme(plot.title = element_text(hjust = 0.5))
  
  
    
# ------------------------------------------------------------------------------------------
# 1000 simulations
# ------------------------------------------------------------------------------------------
  
  Nsim <- 1000
  
  d_100_1 <- matrix(nrow = Nsim, ncol = 4)
  for (i in 1:Nsim) d_100_1[i, 1] <- closest(mysample(100, 2))[1]
  for (i in 1:Nsim) d_100_1[i, 2] <- closest(mysample(100, 3))[1]
  for (i in 1:Nsim) d_100_1[i, 3] <- closest(mysample(100, 5))[1]
  for (i in 1:Nsim) d_100_1[i, 4] <- closest(mysample(100, 10))[1]
  
  d_1000_1 <- matrix(nrow = Nsim, ncol = 4)
  for (i in 1:Nsim) d_1000_1[i, 1] <- closest(mysample(1000, 2))[1]
  for (i in 1:Nsim) d_1000_1[i, 2] <- closest(mysample(1000, 3))[1]
  for (i in 1:Nsim) d_1000_1[i, 3] <- closest(mysample(1000, 5))[1]
  for (i in 1:Nsim) d_1000_1[i, 4] <- closest(mysample(1000, 10))[1]
  
  d_10000_1 <- matrix(nrow = Nsim, ncol = 4)
  for (i in 1:Nsim) d_10000_1[i, 1] <- closest(mysample(10000, 2))[1]
  for (i in 1:Nsim) d_10000_1[i, 2] <- closest(mysample(10000, 3))[1]
  for (i in 1:Nsim) d_10000_1[i, 3] <- closest(mysample(10000, 5))[1]
  for (i in 1:Nsim) d_10000_1[i, 4] <- closest(mysample(10000, 10))[1]
  
  d_sample_1 <- rbind(
    apply(d_100_1, 2, median),
    apply(d_1000_1, 2, median),
    apply(d_10000_1, 2, median)
  )
  rownames(d_sample_1) <- n
  colnames(d_sample_1) <- p

# ------------------------------------------------------------------------------------------
# 10000 simulations
# ------------------------------------------------------------------------------------------
  Nsim <- 10000
  
  d_100_2 <- matrix(nrow = Nsim, ncol = 4)
  for (i in 1:Nsim) d_100_2[i, 1] <- closest(mysample(100, 2))[1]
  for (i in 1:Nsim) d_100_2[i, 2] <- closest(mysample(100, 3))[1]
  for (i in 1:Nsim) d_100_2[i, 3] <- closest(mysample(100, 5))[1]
  for (i in 1:Nsim) d_100_2[i, 4] <- closest(mysample(100, 10))[1]
  
  d_1000_2 <- matrix(nrow = Nsim, ncol = 4)
  for (i in 1:Nsim) d_1000_2[i, 1] <- closest(mysample(1000, 2))[1]
  for (i in 1:Nsim) d_1000_2[i, 2] <- closest(mysample(1000, 3))[1]
  for (i in 1:Nsim) d_1000_2[i, 3] <- closest(mysample(1000, 5))[1]
  for (i in 1:Nsim) d_1000_2[i, 4] <- closest(mysample(1000, 10))[1]
  
  d_10000_2 <- matrix(nrow = Nsim, ncol = 4)
  for (i in 1:Nsim) d_10000_2[i, 1] <- closest(mysample(10000, 2))[1]
  for (i in 1:Nsim) d_10000_2[i, 2] <- closest(mysample(10000, 3))[1]
  for (i in 1:Nsim) d_10000_2[i, 3] <- closest(mysample(10000, 5))[1]
  for (i in 1:Nsim) d_10000_2[i, 4] <- closest(mysample(10000, 10))[1]
  
  d_sample_2 <- rbind(
    apply(d_100_2, 2, median),
    apply(d_1000_2, 2, median),
    apply(d_10000_2, 2, median)
  )
  rownames(d_sample_2) <- n
  colnames(d_sample_2) <- p
  

# ------------------------------------------------------------------------------------------
# Summary Plots
# ------------------------------------------------------------------------------------------
  P <- 1:20
  N <- 100
  Nsim <- 1000
  Samples <- list()
  MinD <- matrix(nrow = Nsim, ncol = 20)
  for (nn in 1:Nsim) {
    for (pp in P) Samples[[pp]] <- mysample(N, pp)
    MinD[nn,] <- lapply(Samples, function(s) closest(s)[1]) %>% unlist
  }
  
  data.frame(
    P, 
    Sample = apply(MinD, 2, median),
    True = closest_m(N, P)
  ) %>% ggplot() +
    geom_line(aes(P, True,  color = 'True', linetype = 'True'), size = 1) +
    geom_line(aes(P, Sample, color = 'Sample', linetype = 'Sample'), size = 2) + 
    theme_light() + theme(legend.position = 'top') +
    scale_color_manual(values = c(Sample = 'darkblue', True = 'orange')) +
    scale_linetype_manual(values = c(Sample = 3, True = 1)) +
    labs(y = 'Median Distance', color = '', linetype = '')
  

  
  
  