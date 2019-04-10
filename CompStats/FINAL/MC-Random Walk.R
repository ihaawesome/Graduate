library(dplyr)
library(ggplot2)
library(reshape2)
setwd("C:/Users/HK/Desktop/18-2/CS/FINAL")

############################## Self Avoiding Walks ##############################

# 현재 위치에서 4방향으로의 좌표
  neighbor <- function(i, j) {
    matrix(c(i+1, j,
             i-1, j,
             i, j-1,
             i, j+1), ncol = 2, byrow = T)
  }
  
# SAW 알고리즘
# wt (weight) = t-1시점에서 인접한 좌표 각각 방문한 적이 있으면 1, 없으면 0
# weight가 모두 0이 되면 trapped
  SAW <- function(t) {
    df <- matrix(0, ncol = 2)
    ti <- 1
    
    while(ti <= t) {
      findnext <- neighbor(df[ti,1], df[ti,2])
      avoid <- matrix(nrow = nrow(df), ncol = 4)
      for (fi in 1:4) {
        for (di in 1:nrow(df)) {
          avoid[di,fi] <- ifelse(all(findnext[fi,]==df[di,]), 0, 1)
        }
      }
      wt <- apply(avoid, 2, function(x) all(x==1))
      
      if (any(wt)) {
        choose <- sample(4, 1, prob = wt)
        df <- rbind(df, findnext[choose,])
        ti <- ti + 1
      } else { ti <- t + 1 }
    }
    
    df <- data.frame(df) ; colnames(df) <- c("x", "y")
    return(df)
  }
  
  # Path 그리기
  gpath <- function(df, t) {
    n <- nrow(df)
    title1 <- ifelse(n < t+1, "Trapped!", "Flying...")
    title2 <- ifelse(n < t+1, paste("stopped when t =", n-1), "")
    
    ggplot(df) + theme_gray() +
      theme(plot.title = element_text(hjust = 0.5), 
            plot.subtitle = element_text(hjust = 0.5), legend.position = "") + 
      labs(title = title1, subtitle = title2) +
      geom_hline(aes(yintercept = 0), linetype = "dotted", color = "lightgray") +
      geom_vline(aes(xintercept = 0), linetype = "dotted", color = "lightgray") +
      geom_path(aes(x, y), color = "red", size = 1) +
      geom_point(aes(0, 0, size = 3), color = "lightgray") +
      geom_text(aes(0, 0, label = "START")) +
      geom_point(aes(x[n], y[n], size = 3), color = "lightgray", shape = 23) +
      geom_text(aes(x[n], y[n], label = "NOW"))
  }

# 시뮬레이션
  SIM <- function(f, t, niter=10000) {
    out <- list()
    for (it in 1:niter) { out[[it]] <- f(t) ; if(it%%10==0) print(it) }
    
    nstep <- lapply(out, function(d) nrow(d)) %>% as.numeric
    Dt_sq <- lapply(out, function(d) (d[nrow(d),1])^2 + (d[nrow(d),2])^2) %>% as.numeric
    Dt_m <- lapply(out, function(d) abs(d[nrow(d),1]) + abs(d[nrow(d),2])) %>% as.numeric
    
    return(list(simulation = out, summary = data.frame(nstep, Dt_sq, Dt_m)))
  }

# Trial 

  saw30 <- SAW(30) ; gpath(saw30, 30)
  saw100 <- SAW(100) ; gpath(saw100, 100)
  saw500 <- SAW(500) ; gpath(saw500, 500)


  MC_SAW20 <- SIM(f = SAW, t = 20, niter = 10000)
  summary(MC_SAW20$summary) ; apply(MC_SAW20$summary[2:3], 2, sd) # Dt: mean = 5.784, sd = 4.0144
  sum(MC_SAW20$summary$nstep != 21) / 10000                       # 10000번 중 trapped된 비율 = 0.0877
  sum(MC_SAW20$summary$Dt_sq) / 897697164                         # mean-square distance = 49.5699
  
  ggplot(MC_SAW20[[2]]) + geom_bar(aes(x = Dt_m)) + theme_light() + 
    facet_wrap(~ifelse(nstep==21, "Flying", "Trapped"))
  
  MC_SAW30 <- SIM(f = SAW, t = 30, niter = 10000)
  summary(MC_SAW30$summary) ; apply(MC_SAW30$summary[2:3], 2, sd) # Dt: mean = 10.06, sd = 5.1088
  sum(MC_SAW30$summary$nstep != 31) / 10000                       # 10000번 중 trapped된 비율 = 0.1971
  
  ggplot(MC_SAW30[[2]]) + geom_bar(aes(x = Dt_e)) + theme_light() + 
    facet_wrap(~ifelse(nstep==31, "Flying", "Trapped"))
  
  MC_SAW100 <- SIM(f = SAW, t = 100, niter = 10000)
  summary(MC_SAW100$summary) ; apply(MC_SAW100$summary[2:3], 2, sd) # Dt: mean = 10.05, sd = 9.0388
  sum(MC_SAW100$summary$nstep != 101) / 10000                       # 10000번 중 trapped된 비율 = 0.7884
  mean(MC_SAW100$summary$Dt_e^2)
  
  ggplot(MC_SAW100[[2]]) + geom_bar(aes(x = Dt_e)) + theme_light() + 
    facet_wrap(~ifelse(nstep==101, "Flying", "Trapped"))
  
  MC_SAW3 <- SIM(f = SAW, t = 10000, niter = 1000)
  summary(MC_SAW3$summary$Dt) ; sd(MC_SAW3$summary$Dt)              # Dt: mean = 15.54, sd = 12.2606
  summary(MC_SAW3$summary$nstep)
  sum(MC_SAW3$summary$nstep <= 334) / 1000                          # 1000번 중 trapped된 비율 = 1
  
  ggplot(MC_SAW3[[2]]) + geom_bar(aes(x = Dt)) + theme_light()
  
  which.max(MC_SAW3$summary$Dt) ; which.max(MC_SAW3$summary$nstep)
  
  gpath(MC_SAW3$simulation[[795]], t = 10000)
  gpath(MC_SAW3$simulation[[722]], t = 10000)
  
  
  MC_SAW4 <- SIM(f = SAW, t = 10000, niter = 10000)
  which.max(MC_SAW4$summary$Dt) ; which.max(MC_SAW4$summary$nstep)
  summary(MC_SAW4$summary$Dt) ; sd(MC_SAW4$summary$Dt)              # Dt: mean = 15.02, sd = 11.71917
  summary(MC_SAW4$summary$nstep)
  
  ggplot(MC_SAW4[[2]]) + geom_bar(aes(x = Dt)) + theme_light()

  
############################## Random Walk ##############################
  
  RW <- function(t) {
    df <- matrix(0, ncol = 2)
    for (ti in 1:t) {
      findnext <- neighbor(df[ti,1], df[ti,2])
      choose <- sample(4, 1)
      df <- rbind(df, findnext[choose,])
    }
    df <- data.frame(df) ; colnames(df) <- c("x", "y")
    return(df)
  }
  
  deletepath <- function(df) { ifelse(nrow(filter(count(df, x, y), n > 1)) > 0, 1, 0) }
  
# Trial
  rw1 <- RW(t = 30) ; gpath(rw1, 30)
  rw2 <- RW(t = 1000) ; gpath(rw2, 30)
  
  MC_RW10 <- SIM(f = RW, t = 10, niter = 10000)
  MC_RW30 <- SIM(f = RW, t = 30, niter = 10000)
  MC_RW100 <- SIM(f = RW, t = 100, niter = 1000)
  MC_RW1000 <- SIM(f = RW, t = 1000, niter = 10000)
  
  summary(MC_RW1000[[2]])
  ggplot(MC_RW1000[[2]]) + geom_bar(aes(Dt)) + theme_light() 
  
  
  delete10 <- lapply(MC_RW10[[1]], function(d) deletepath(d)) %>% as.numeric
  MC_RW_SAW10 <- list(simulation = MC_RW10$simulation[which(delete10==0)],
                      summary = MC_RW10$summary[which(delete10==0),])
  
  summary(MC_RW_SAW10[[2]])
  ggplot(MC_RW_SAW10[[2]]) + geom_bar(aes(Dt)) + theme_light()
  
  delete1 <- lapply(MC_RW1[[1]], function(d) deletepath(d)) %>% as.numeric
  sum(delete1)
  
  delete3 <- lapply(MC_RW3[[1]], function(d) deletepath(d)) %>% as.numeric
  sum(delete3)
  
  
  log(mean(MC_SAW_D1$summary$Dt), base = 100)/2


############################## Pivot Algorithm (Try) ##############################
  
  pivot <- function(pivot, point) { # point as vector
    P <- t(pivot) ; A <- t(point)
    
    RX <- matrix(c(1,  0,
                   0, -1), ncol = 2, byrow = T)
    RY <- matrix(c(-1, 0,
                    0, 1), ncol = 2, byrow = T)
    R90 <- matrix(c(0, -1,
                    1,  0), ncol = 2, byrow = T)
    R180 <- matrix(c(-1,  0,
                      0, -1), ncol = 2, byrow = T)
    R270 <- matrix(c( 0, 1,
                     -1, 0), ncol = 2, byrow = T)
    
    r <- sample(5, 1)
    if (r==1) { A_ <- (A-Pm) %*% RX + Pm }
    if (r==2) { A_ <- (A-Pm) %*% RY + Pm }
    if (r==3) { A_ <- (A-Pm) %*% R90 + Pm }
    if (r==4) { A_ <- (A-Pm) %*% R180 + Pm }
    if (r==5) { A_ <- (A-Pm) %*% R270 + Pm }
    
    return(list(pivoted = A_, rotate = r))
  }
  
# 시뮬레이션  
  SIM_P <- function(f, t, niter = 10000) {
    out <- list()
    for (it in 1:niter) { 
      out[[it]] <- f(t)
      pivots <- pivot(as.matrix(out[[it]]))
      p <- pivots[[2]] ; n <- nrow(pivots[[1]])
  
      reject <- matrix(nrow = p-1, ncol = n-p)
      for (b in 1:(p-1)) { for (a in (p+1):n) { 
        reject[b, a-p] <- as.numeric(all(pivots[[1]][b,] == pivots[[1]][a,])) } }
      if(sum(reject) == 0) { out[[it]] <- pivots[[1]] }
      
      if(it%%10==0) { print(it) }
    }
  
    nstep <- lapply(out, function(d) nrow(d)) %>% as.numeric
    Dt_sq <- lapply(out, function(d) (d[nrow(d),1])^2 + (d[nrow(d),2])^2) %>% as.numeric
    Dt_m <- lapply(out, function(d) abs(d[nrow(d),1]) + abs(d[nrow(d),2])) %>% as.numeric
    
    return(list(simulation = out, summary = data.frame(nstep, Dt_sq, Dt_m)))
  }

# Trial   
  MC_PV20 <- SIM_P(SAW, 20, niter=10000)
  MC_PV30 <- SIM_P(SAW, 30, niter=10000)
  
  summary(MC_SAW20$summary) ; apply(MC_SAW20$summary[2:3], 2, sd)
  summary(MC_PV20$summary) ; apply(MC_PV20$summary[2:3], 2, sd)
  
  ggplot(MC_PV20$summary) + geom_bar(aes(Dt_sq))
  ggplot(MC_SAW20$summary) + geom_bar(aes(Dt_sq))
  
  summary(MC_PV30$summary) ; apply(MC_PV30$summary[2:3], 2, sd) 
  ggplot(MC_PV30$summary) + geom_bar(aes(Dt_sq)) + facet_wrap(~(nstep==31))
  
  gpath(as.data.frame(pivot(as.matrix(saw30))[[1]]), 31)
  
  
