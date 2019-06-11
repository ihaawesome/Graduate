###### CH10 Exercises
library(ggplot2)
library(dplyr)

##### Exercise 10.7 #####
data(USArrests)
data.scaled <- scale(USArrests)
d1 <- dist(data.scaled)^2
d2 <- as.dist(1-cor(t(data.scaled)))
summary(d2/d1)

plot(d2/d1, col = 'darkgray')
plot(density(d2/d1), main = '', xlab = 'd2/d1')
boxplot(d2/d1, horizontal = T)


##### Exercise 10.10 #####
# (a)
n <- 20 ; p <- 50
set.seed(1010)
X <- rbind(
  matrix(rnorm(n*p, mean = 1*0.5), n, p),
  matrix(rnorm(n*p, mean = 2*0.5), n, p),
  matrix(rnorm(n*p, mean = 3*0.5), n, p)
)
y <- as.factor(rep(1:3, each = n))

plot(density(rnorm(100000,1*0.5)), main = 'True Density of each 3 Class', xlab = '', col = 'red1', lwd = 2)
lines(density(rnorm(100000,2*0.5)), col = 'green3', lwd = 2)
lines(density(rnorm(100000,3*0.5)), col = 'blue1', lwd = 2)

# (b)
summary(pc.out <- prcomp(X))
plot(pc.out, main = 'Variance Explained')

pc.two <- as.data.frame(pc.out$x[,1:2])
ggplot(pc.two) + theme_bw() + theme(legend.position = 'top') +
  geom_point(aes(PC1, PC2, shape = y, color = y), size = 2)

# (c)
set.seed(10)
(km.3class <- kmeans(X, 3, nstart = 100))
cl.3class <- factor(c(3:1)[km.3class$cluster], levels = c(1:3))

table(true = y, kmeans = cl.3class)
ggplot(pc.two) + theme_bw() + theme(legend.position = 'top') +
  geom_point(aes(PC1, PC2, shape = cl.3class, color = cl.3class), size = 2)
  
# (d)
set.seed(10)
(km.2class <- kmeans(X, 2, nstart = 100))
cl.2class <- factor(c(3,1)[km.2class$cluster], levels = c(1:3)) 

ggplot(pc.two) + theme_bw() + theme(legend.position = 'top') +
  geom_point(aes(PC1, PC2, shape = cl.2class, color = cl.2class), size = 2)
table(true = y, kmeans = cl.2class)

# (e)
set.seed(10)
(km.4class <- kmeans(X, 4, nstart = 100))
cl.4class <- factor(c(2,1,3,4)[km.4class$cluster], levels = c(1:4)) 

ggplot(pc.two) + theme_bw() + theme(legend.position = 'top') +
  geom_point(aes(PC1, PC2, shape = cl.4class, color = cl.4class), size = 2)
table(true = y, kmeans = cl.4class)

# (f)
set.seed(10)
(km.pc <- kmeans(pc.two, 3, nstart = 100))
cl.pc <- as.factor(c(3,1,2)[km.pc$cluster])

ggplot(pc.two) + theme_bw() + theme(legend.position = 'top') +
  geom_point(aes(PC1, PC2, shape = cl.pc, color = cl.pc), size = 2)
table(true = y, kmeans = cl.pc)

# (g) / compare to (b)
Xsc <- scale(X)
set.seed(10)
(km.scaled <- kmeans(Xsc, 3, nstart = 100))
cl.scaled <- factor(c(1,3,2)[km.scaled$cluster], levels = c(1:3))
ggplot(pc.two) + theme_bw() + theme(legend.position = 'top') +
  geom_point(aes(PC1, PC2, shape = cl.scaled, color = cl.scaled), size = 2)
table(true = y, kmeans = cl.pc)

