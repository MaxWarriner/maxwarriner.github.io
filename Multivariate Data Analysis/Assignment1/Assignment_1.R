# Script for Assignment 1: Max Warriner
# Question 1

load("C:/Users/12697/Documents/Multivariate Data Analysis/Assignment-1-2026/DigitsF.Rdata")

#Merge data together
dat <- cbind(rep(1, 1000), t(digitsF[,,1]))
for (digit in 2:10){
  dat <- rbind(dat, cbind(rep(digit, 1000), t(digitsF[,,digit])))
}

#Compute a Mahalanobis distance criterion

# a 

X = dat[,2:31]
labs=dat[,1]

mh <- function(X,labs){
  p = ncol(X) #dimension of the data
  
  J = length(table(labs)) #number of categories
  
  w_cov = array(0,c(p,p)) # weighted covariance placeholder
  
  N = nrow(X) # overall sample size
  
  xbar = apply(X,2,mean) #xbar vector
  
  for(j in 1:J) {
    
    Nj = length(labs[labs==j]) #sample size in each group
    
    w_cov = w_cov +(Nj-1)*var(X[labs==j,])/N #weight the overall cov matrix by the sample size in each group
    
  }
  
  mh=0
  
  for(j in 1:J) {
    
    xbar_j = apply(X[labs==j,],2,mean) #xbar vector for the j_th group
    
    mh = mh + t((xbar_j - xbar))%*%solve(w_cov)%*%(xbar_j - xbar) 
    
  }
  
  return(list(mh, w_cov))
  
}

#b
(statistic <- mh(X, labs)[[1]][1,1])
w_cov <- mh(X, labs)[[2]]

#c

NS = 1000

library(MASS)

# bootstrapping

stats_boot <- rep(NA, NS)

J = 10

boot_X <- X
for (j in 1:J){
  
  mu_j = apply(X[labs == j,], 2, mean)
  boot_X[labs == j,] = sweep(X[labs == j,], 2, STATS = mu_j, FUN = "-")
  
}

boot_X <- boot_X + mean(X)
set.seed(101)
for (i in 1:NS){
  random_X <- array(0, dim = c(10000,30))
  for (j in 1:J){
  X_j = boot_X[labs == j,]
  random_X_j <- X_j[sample(1000, replace=T),]
  random_X[(1+(j-1)*1000):(j*1000),] <- random_X_j
  }
  
  stats_boot[i] <- mh(random_X, labs = labs)[[1]][1,1]
}

#d
stats_boot[1:3]

(p <- mean(stats_boot > statistic))


#e

# Permutation Test simulation

stats_permute <- rep(NA, NS)

set.seed(101)
for (i in 1:NS){
  random_labs <- sample(labs, replace = F)
  stats_permute[i] <- mh(X, random_labs)[[1]][1,1]
}


#f 
stats_permute[1:3]
(p <- mean(stats_permute > statistic))


#g
ks.test(stats_boot, stats_permute)

qqplot(stats_boot, stats_permute, 
       xlab = "Boostrapping Distribution", 
       ylab = "Permutation Distribution", 
       main = "qq-plot Comparison", 
       sub = "K-S Test: D = 0.042, p = 0.341")
abline(0, 1, col = "red")




# Question 2
library(tidyverse)

dat <- USArrests

#a
cor <- cor(dat)
(lambda = -log(det(cor)))


#b
boots <- 5000
boot_lambdas <- rep(NA, boots)

set.seed(101)
for (i in (1:boots)){
  boot_X <- dat[sample(1:50, replace = T),]
  boot_cor <- cor(boot_X)
  boot_lambdas[i] <- -log(det(boot_cor))
}

(se_lambda <- sd(boot_lambdas) / sqrt(boots))

fitdistr(boot_lambdas, "gamma")

plot_dat <- tibble(lambdas = boot_lambdas) 

plot_dat <- plot_dat |>
  mutate(density = dgamma(lambdas, shape = 32.7092516, rate = 15.9143465))

ggplot(data = plot_dat) + 
  geom_histogram(aes(x = lambdas, y = after_stat(density), fill = "Bootstrapped Distribution"), 
                 bins = 50, alpha = 0.75) + 
  geom_line(aes(x = lambdas, y = density, color = "Fitted Gamma")) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlab(bquote(hat(lambda))) + 
  ylab('Density') + 
  labs(fill = "",
       colour = "") + 
  scale_fill_manual(values = c("Bootstrapped Distribution" = "lightblue2")) + 
  geom_point(aes(x = lambda, y = dgamma(lambda, shape = 32.7092516, rate = 15.9143465), color = "Test Statistic"), 
             size = 3)


#c
sims <- 5000
norm_lambdas <- rep(NA, sims)
xbar <- apply(dat, 2, mean)
sigma_hat <- diag(diag(cov(dat)))

library(MASS)
set.seed(101)
for(i in 1:sims){
  random_X <- mvrnorm(n = 50, mu = xbar, Sigma = sigma_hat)
  random_cor <- cor(random_X)
  norm_lambdas[i] <- -log(det(random_cor))
}


fitdistr(norm_lambdas, "gamma")

null_plot_dat <- tibble(lambdas = norm_lambdas) 

null_plot_dat <- null_plot_dat |>
  mutate(density = dgamma(lambdas, shape = 2.94872143, rate = 23.09649930))

ggplot(data = null_plot_dat) + 
  geom_histogram(aes(x = lambdas, y = after_stat(density), fill = "Null Distribution"), 
                 bins = 100, alpha = 0.75) + 
  geom_line(aes(x = lambdas, y = density, color = "Fitted Gamma")) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlab(bquote(hat(lambda))) + 
  ylab('Density') + 
  labs(fill = "",
       colour = "") + 
  scale_fill_manual(values = c("Null Distribution" = "lightblue2"))
  # geom_point(aes(x = lambda, y = dgamma(lambda, shape = 3.11104268, rate = 24.21261095), color = "Test Statistic"), 
  #            size = 3)

(null_p <- mean(norm_lambdas > lambda))

#d
boots <- 5000
scaled_boots <- rep(NA, 5000)
scaled_dat <- scale(dat)
Z_stats <- as.vector(scaled_dat)

set.seed(101)
for(i in 1:boots){
  random_Zs <- sample(Z_stats, size = 200, replace = T)
  random_X <- matrix(random_Zs, nrow = 50)
  random_cor <- cor(random_X)
  scaled_boots[i] <- -log(det(random_cor))
}

fitdistr(scaled_boots, "gamma")

scaled_plot_dat <- tibble(lambdas = scaled_boots) 

scaled_plot_dat <- scaled_plot_dat |>
  mutate(density = dgamma(lambdas, shape = 3.00364986, rate = 23.37034828))

ggplot(data = scaled_plot_dat) + 
  geom_histogram(aes(x = lambdas, y = after_stat(density), fill = "Scaled Bootstrap Distribution"), 
                 bins = 100, alpha = 0.75) + 
  geom_line(aes(x = lambdas, y = density, color = "Fitted Gamma")) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlab(bquote(hat(lambda))) + 
  ylab('Density') + 
  labs(fill = "",
       colour = "") +
  scale_fill_manual(values = c("Scaled Bootstrap Distribution" = "lightblue2")) + 
  scale_color_manual(values = c("Fitted Gamma" = "red"))

(scaled_p <- mean(scaled_boots > lambda))


#e
ks.test(scaled_boots, norm_lambdas)



# Question 3

#d

rm(list = ls())


nih.palette <- colorRampPalette(c(1,4,5,"lawngreen",7,2),space = "rgb")
gray.palette <- colorRampPalette(c(1,"white"),space = "rgb")


# Binary Image Read

fdata='cald0814Attn.F01.floatImage'
att=readBin(fdata,what='numeric',size=4,n=128*128*35,endian='big')
summary(att)
att=array(att,c(128,128,35))
layout(matrix(c(1,2,3,4,5,5), nrow = 3, byrow = T), heights = c(1,1,0.3))
image(att[,rev(1:128),17],main="PET Attenuation",col=gray.palette(100),axes=F,asp=1)
text(32/128,110/128,'Transverse',col='yellow')
brain_vals <- att[att>0]
hist(c(brain_vals),main="Attenuation",breaks=50)
image(att[,64,],col=gray.palette(100),axes=F,asp=2/4)
text(5/35,110/128,'Sagittal',col='yellow')
image(att[64,,],col=gray.palette(100),axes=F,asp=2/4)
text(5/35,110/128,'Coronal',col='yellow')

zlim <- range(brain_vals)

par(mar = c(2,1,2,1))
image(z = matrix(seq(zlim[1], zlim[2], length = 100), ncol = 1), col = gray.palette(100), axes = T, main = "Tissue density")
axis(4)


