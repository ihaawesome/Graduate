library(dplyr)
library(ggplot2)
library(extraDistr)
# HW7
# 1. slash dist
  m = 100000 ; n = 5000
  
# use slash as e, generate normal
# target f: normal, envelope g: slash  
  ym <- rslash(m)
  usiw <- dnorm(ym) / dslash(ym) ; siw <- usiw / sum(usiw)
  xn <- sample(ym, n, replace = T, prob = siw)
  
  hist(xn, freq = F, breaks = 50, col = "gray", border = "white", 
       xlab = "x", ylab = "normal density", main = "generate normal from slash")
  lines(sort(xn), dnorm(sort(xn)), col = "blue", lwd = 2)

# target f: slash, envelope g: normal
  ym <- rnorm(m)
  usiw <- dslash(ym) / dnorm(ym) ; siw <- usiw / sum(usiw)
  xn <- sample(ym, n, replace = T, prob = siw)
  
  hist(xn, freq = F, breaks = 50, col = "gray", border = "white",
       xlab = "y", ylab = "slash density", main = "generate slash from normal")
  lines(sort(xn), dslash(sort(xn)), col = "blue", lwd = 2)

  gg <- data.frame(x = seq(-5, 5, length = 100)) %>% mutate(slash = dslash(x), normal = dnorm(x))
  ggplot(gg) + theme_test() +
    geom_area(aes(x, slash, fill = "slash"), color = "black", alpha = 0.7) + 
    geom_area(aes(x, normal, fill = "normal"), color = "black", alpha = 0.5) +
    scale_fill_brewer("", palette = "Blues") + labs(y = "density") + theme(legend.position = "top")
  
  
# 2. Bayesian Inference
  data2 <- c(8, 3, 4, 3, 1, 7, 2, 6, 2, 7)
  Lfunc <- function(lambda) { 
    lik = c() ;for (i in 1:length(lambda)) lik[i] <- prod(dpois(data2, lambda[i]))
    return(lik) 
  }
  
  myfunc1 <- function(m, n) { 
    l <- rlnorm(m, log(4), 0.5) 
    usiw <- Lfunc(l)
    l.re <- sample(l, n, replace = T, prob = usiw)
    return(l.re)
  }
  
  myfunc2 <- function(n) {
    l <- rlnorm(n, log(4), 0.5) ; u <- runif(n)
    l.keep <- l[(u <= Lfunc(l)/Lfunc(4.3))]
    return(l.keep)
  }

  m = 100000 ; n = 10000
  l.sir <- myfunc1(m, n)
  l.rs <- myfunc2(m)
  l.rs.n <- sample(l.rs, n)
  
  l.qq <- data.frame(l.rs = sort(l.rs.n), l.sir = sort(l.sir))
  ggplot(l.qq) + 
    theme_test() + theme(plot.title = eleyment_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    geom_smooth(aes(l.rs, l.sir), size = 2, color = "skyblue", method = "lm") + geom_point(aes(l.rs, l.sir)) +
    labs(title = "Q-Q Plot", 
         subtitle = paste("Rej. Sampling vs SIR (", round(summary(lm(l.sir ~ l.rs, l.qq))$r.squared, 4)*100, "%)", sep = ""), 
         x = "Rejection Sampling", y = "SIR")
  
  ggplot() + theme_test() + labs(x = "lambda") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "top") +
    geom_density(aes(l.sir, fill = "SIR"), size = 1, alpha = 0.7) +
    geom_density(aes(l.rs, fill = "Rej.Sampling"), size = 1, alpha = 0.3) + 
    scale_fill_brewer("", palette = "Blues")
  
  
  n10000 = myfunc1(m, 10000) ; n1000 = myfunc1(m, 1000) ; n50000 = myfunc1(m, 50000)
  ggplot() + theme_minimal() + scale_color_brewer("n", palette = "Blues") +
    geom_density(aes(n1000, color = "1000"), size = 3) + 
    geom_density(aes(n10000, color = "10000"), size = 2) +
    geom_density(aes(n50000, color = "50000"), size = 1) +
    labs(x = "lambda", title = "Density Plot") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
  
# 3. Network Failure Prob.
  
  myroute <- function(x) {
    m <- matrix(c(1, x[1], x[2], x[3], x[4], 0, 0, 0, 0, 0,
                  x[1], 1, x[5], 0, 0, x[8], x[9], 0, 0, 0,
                  x[2], x[5], 1, x[6], 0, 0, x[10], 0, 0, 0,
                  x[3], 0, x[6], 1, x[7], 0, 0, x[11], 0, 0,
                  x[4], 0, 0, x[7], 1, 0, 0, x[12], x[13], 0,
                  0, x[8], 0, 0, 0, 1, x[14], 0, 0, x[17],
                  0, x[9], x[10], 0, 0, x[14], 1, x[15], 0, x[18],
                  0, 0, 0, x[11], x[12], 0, x[15], 1, x[16], x[19],
                  0, 0, 0, 0, x[13], 0, 0, x[16], 1, x[20],
                  0, 0, 0, 0, 0, x[17], x[18], x[19], x[20], 1), 
                ncol = 10, byrow = T)
    res <- ifelse((m %*% m %*% m %*% m %*% m %*% m %*% m %*% m %*% m)[1, 10] > 0, 0, 1)
    return(res)
  }
  
  n = 100000
  p = 0.05 ; p.star = 0.2
  
  x.mc <- matrix(sample(c(0, 1), 20*n, replace = T, prob = c(p, 1-p)), ncol = n) # 0: broken
  hx.mc <- apply(x.mc, myroute, MARGIN = 2)
  
  x.is <- matrix(sample(c(0, 1), 20*n, replace = T, prob = c(p.star, 1-p.star)), ncol = n)  
  hx.is <- apply(x.is, myroute, MARGIN = 2)   
  bx.is <- 20 - apply(x.is, sum, MARGIN = 2)
  usiw <- ((1-p)/(1-p.star))^20 * ((p*(1-p.star))/(p.star*(1-p)))^bx.is
  siw <- usiw / sum(usiw)
  
  sum(hx.mc) ; sum(hx.is) # of broken edges
  
  mu.mc <- mean(hx.mc)
  mu.is.star <- mean(hx.is*usiw)
  mu.is <- sum(hx.is*siw)
  c(mu.mc, mu.is.star, mu.is) %>% round(8)
  
  se.mc <- sqrt( var(hx.mc) / n )
  se.is.star <- sqrt( var(hx.is*usiw) / n )
  se.is <- sqrt( (var(hx.is*usiw) + mu.is^2*var(usiw) - 2*mu.is*cov(hx.is*usiw, usiw)) / n )
  c(se.mc, se.is.star, se.is) %>% round(8)
  
  myplot <- function(x) {
    xpt = c(0, 1, 1, 1, 1, 2, 2, 2, 2, 3) ; ypt = c(0, 1.5, 0.5, -0.5, -1.5, 1.5, 0.5, -0.5, -1.5, 0)
    g <- ggplot() + theme_minimal() + labs(x = "", y = "") + 
          theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
          geom_point(data = data.frame(xpt, ypt), aes(xpt, ypt), size = 5)
      
    if (x[1]==1) { g <- g + geom_line(aes(c(0, 1), c(0, 1.5))) }
    if (x[2]==1) { g <- g + geom_line(aes(c(0, 1), c(0, 0.5))) }
    if (x[3]==1) { g <- g + geom_line(aes(c(0, 1), c(0, -0.5))) }
    if (x[4]==1) { g <- g + geom_line(aes(c(0, 1), c(0, -1.5))) }
    if (x[5]==1) { g <- g + geom_line(aes(c(1, 1), c(1.5, 0.5))) }
    if (x[6]==1) { g <- g + geom_line(aes(c(1, 1), c(0.5, -0.5))) }
    if (x[7]==1) { g <- g + geom_line(aes(c(1, 1), c(-0.5, -1.5))) }
    if (x[8]==1) { g <- g + geom_line(aes(c(1, 2), c(1.5, 1.5))) }
    if (x[9]==1) { g <- g + geom_line(aes(c(1, 2), c(1.5, 0.5))) }
    if (x[10]==1) { g <- g + geom_line(aes(c(1, 2), c(0.5, 0.5))) }
    if (x[11]==1) { g <- g + geom_line(aes(c(1, 2), c(-0.5, -0.5))) }
    if (x[12]==1) { g <- g + geom_line(aes(c(1, 2), c(-1.5, -0.5))) }
    if (x[13]==1) { g <- g + geom_line(aes(c(1, 2), c(-1.5, -1.5))) }
    if (x[14]==1) { g <- g + geom_line(aes(c(2, 2), c(1.5, 0.5))) }
    if (x[15]==1) { g <- g + geom_line(aes(c(2, 2), c(0.5, -0.5))) }
    if (x[16]==1) { g <- g + geom_line(aes(c(2, 2), c(-0.5, -1.5))) }
    if (x[17]==1) { g <- g + geom_line(aes(c(2, 3), c(1.5, 0))) }
    if (x[18]==1) { g <- g + geom_line(aes(c(2, 3), c(0.5, 0))) }
    if (x[19]==1) { g <- g + geom_line(aes(c(2, 3), c(-0.5, 0))) }
    if (x[20]==1) { g <- g + geom_line(aes(c(2, 3), c(-1.5, 0))) }
    
    g <- g + geom_point(aes(c(0, 3), c(0, 0)), color = "orange", size = 5)
    return(g)
  }

  myplot(rep(1, 20))
  
# change p.star  
  p.star <- 0.1 ; x.is.new <- matrix(sample(c(0, 1), 20*n, replace = T, prob = c(p.star, 1-p.star)), ncol = n) 
  hx.is.new <- apply(x.is.new, myroute, MARGIN = 2)   
  bx.is.new <- 20 - apply(x.is.new, sum, MARGIN = 2)
  usiw <- ((1-p)/(1-p.star))^20 * ((p*(1-p.star))/(p.star*(1-p)))^bx.is.new
  siw <- usiw / sum(usiw)
  
  mu.is.star.new <- sum(hx.is.new*usiw) / n
  se.is.star.new <- sqrt( var(hx.is.new*usiw) / n )
  c(mu.is.star.new, se.is.star.new)
  
  
  
  
  
  