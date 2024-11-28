# 2024-11-25
# 4 samples with normal conditioning on normal populations
# basis function 5
# estimating mean
# Doubly robust kernel imputation
# New formulation of DEL, idea 4, hh's algorithm


rm(list=ls())
setwd("/home/hhe/2024/mean/nn_5")

source("algorithms.R")
source("drm_algorithms.R")

library(doRNG)
library(doParallel)
library(dplyr)

set.seed(100)

for(n_obs in c(100,200,500)){
  
# n_obs <- 1000
n_samples <- rep(x=n_obs,times=4)  # sample sizes

theta <- c(5*5+2, 4.5*4+1.5, 1.2-4.2*5.2, 4.5*5.5+0.8)
m=3
d=2


# Define the basis functions
basis_func <- function(x) {
  return(c(1,x, x^2))
}



cl <- makeCluster(22)
registerDoParallel(cl)



results <- foreach(mc = 1:5000, .packages=c("mvtnorm","drmdel","EnvStats","dplyr")) %dorng% {
  
  x0.t <- rnorm(n_samples[1], mean = 5, sd = 1.3)
  x1.t <- rnorm(n_samples[2], mean = 4.5, sd = 1)
  x2.t <- rnorm(n_samples[3], mean = 5.2, sd = 1.5)
  x3.t <- rnorm(n_samples[4], mean = 5.5, sd = 1.8)
  mu0 <- 2 + 5*x0.t
  mu1 <- 1.5 + 4*x1.t
  mu2 <- 1.2 - 4.2*x2.t
  mu3 <- 0.8 + 4.5*x3.t
  y0 <- mapply(rnorm, mean = mu0, MoreArgs = list(n = 1, sd = 8))
  y1 <- mapply(rnorm, mean = mu1, MoreArgs = list(n = 1, sd = 9))
  y2 <- mapply(rnorm, mean = mu2, MoreArgs = list(n = 1, sd = 8.5))
  y3 <- mapply(rnorm, mean = mu3, MoreArgs = list(n = 1, sd = 8.7))
  
  y <- c(y0, y1, y2, y3)
  
  #### Full data estimator ####
  drmfit.full <- drmdel(x=y, n_samples=n_samples, basis_func=5)
  
  p0 <- 0.1+0.9*exp(x0.t-5)/(1+exp(x0.t-5))
  p0 <- ifelse(p0>1,1,p0)
  p1 <- 0.1+0.9*exp(x1.t-5)/(1+exp(x1.t-5))
  p1 <- ifelse(p1>1,1,p1)
  p2 <- 0.1+0.9*exp(x2.t-5)/(1+exp(x2.t-5))
  p2 <- ifelse(p2>1,1,p2)
  p3 <- 0.1+0.9*exp(x3.t-5)/(1+exp(x3.t-5))
  p3 <- ifelse(p3>1,1,p3)
  
  delta_0 <- rbinom(n_samples[1],1,p0)
  # sum(delta_0)/n_samples[1]
  delta_1 <- rbinom(n_samples[2],1,p1)
  # sum(delta_1)/n_samples[2]
  delta_2 <- rbinom(n_samples[3],1,p2)
  # sum(delta_2)/n_samples[3]
  delta_3 <- rbinom(n_samples[4],1,p3)
  # sum(delta_3)/n_samples[4]
  
  delta <- c(delta_0,delta_1,delta_2,delta_3)
  
  # par(mfrow=c(2,2))
  # plot(x0.t,y0)
  # points(x0.t[delta_0==0],y0[delta_0==0],col="red")
  # plot(x1.t,y1)
  # points(x1.t[delta_1==0],y1[delta_1==0],col="red")
  # plot(x2.t,y2)
  # points(x2.t[delta_2==0],y2[delta_2==0],col="red")
  # plot(x3.t,y3)
  # points(x3.t[delta_3==0],y3[delta_3==0],col="red")
  
  
  y0.obs <- y0[delta_0==1]
  y1.obs <- y1[delta_1==1]
  y2.obs <- y2[delta_2==1]
  y3.obs <- y3[delta_3==1]
  
  y0.na <- y0
  y0.na[delta_0==0] <- NA
  y1.na <- y1
  y1.na[delta_1==0] <- NA
  y2.na <- y2
  y2.na[delta_2==0] <- NA
  y3.na <- y3
  y3.na[delta_3==0] <- NA
  y.na <- c(y0.na,y1.na,y2.na,y3.na)
  
  x0.f <- rnorm(n_samples[1], mean = 5, sd = 1.5)
  x1.f <- rnorm(n_samples[2], mean = 5, sd = 1.5)
  x2.f <- rnorm(n_samples[3], mean = 5, sd = 1.5)
  x3.f <- rnorm(n_samples[4], mean = 5, sd = 1.5)
  
  
  #### Estimate propensity score ####
  
  #####
  #### kernel estimate propensity
  # p0.t <- Phat.kern_func(x = x0.t, delta = delta_0, h = n_samples[1]^(-1/5)*sqrt(var(x0.t)))
  # p1.t <- Phat.kern_func(x = x1.t, delta = delta_1, h = n_samples[1]^(-1/5)*sqrt(var(x1.t)))
  # p2.t <- Phat.kern_func(x = x2.t, delta = delta_2, h = n_samples[1]^(-1/5)*sqrt(var(x2.t)))
  # p3.t <- Phat.kern_func(x = x3.t, delta = delta_3, h = n_samples[1]^(-1/5)*sqrt(var(x3.t)))
  #### score estimate propensity
  # p0.t <- Phat.score_func(x = x0.t, delta = delta_0, c = 10)
  # p1.t <- Phat.score_func(x = x1.t, delta = delta_1, c = 10)
  # p2.t <- Phat.score_func(x = x2.t, delta = delta_2, c = 10)
  # p3.t <- Phat.score_func(x = x3.t, delta = delta_3, c = 10)
  #### true propensity
  p0.t <- p0
  p1.t <- p1
  p2.t <- p2
  p3.t <- p3
  
  p.t <- c(p0.t,p1.t,p2.t,p3.t)
  
  
  #####
  p0.f <- Phat.kern_func(x = x0.f, delta = delta_0, h = n_samples[1]^(-1/5)*sqrt(var(x0.f)))
  p1.f <- Phat.kern_func(x = x1.f, delta = delta_1, h = n_samples[1]^(-1/5)*sqrt(var(x1.f)))
  p2.f <- Phat.kern_func(x = x2.f, delta = delta_2, h = n_samples[1]^(-1/5)*sqrt(var(x2.f)))
  p3.f <- Phat.kern_func(x = x3.f, delta = delta_3, h = n_samples[1]^(-1/5)*sqrt(var(x3.f)))
  
  p.f <- c(p0.f,p1.f,p2.f,p3.f)
  
  #### Naive estimator ####
  y.obs <- c(y0.obs, y1.obs, y2.obs, y3.obs)
  n_samples.obs <- c(sum(delta_0==1), sum(delta_1==1), sum(delta_2==1), sum(delta_3==1))
  drmfit.cc <- drmdel(x=y.obs, n_samples=n_samples.obs, basis_func=5)
  
  
  #### Random Kernel Regression Imputation #### 
  
  #####
  y0_hat_t <- imp.randKern_func(h = n_samples[1]^(-1/3)*sqrt(var(x0.t)), x = x0.t, y = y0.na, k = 1, delta = delta_0)
  y1_hat_t <- imp.randKern_func(h = n_samples[2]^(-1/3)*sqrt(var(x1.t)), x = x1.t, y = y1.na, k = 1, delta = delta_1)
  y2_hat_t <- imp.randKern_func(h = n_samples[3]^(-1/3)*sqrt(var(x2.t)), x = x2.t, y = y2.na, k = 1, delta = delta_2)
  y3_hat_t <- imp.randKern_func(h = n_samples[4]^(-1/3)*sqrt(var(x3.t)), x = x3.t, y = y3.na, k = 1, delta = delta_3)
  y_hat_t <- c(y0_hat_t, y1_hat_t, y2_hat_t, y3_hat_t)

  
  # calculating MELE
  par_dim <- m*(d+1)
  par_init <- rep(0, par_dim)
  
  drm_opt <- do.call(optim, c(list(par=par_init, fn=negLDELfn_miss4,
                                   gr=negLDELgr_miss4,
                                   method="BFGS", control=list(maxit=10000),
                                   x=y.na, x_imp=y_hat_t, delta=delta, px_hat=p.t,
                                   n_total=sum(n_samples), n_samples=n_samples, m=m,
                                   basis_func=basis_func,
                                   d=d)))
  mele.tt <- drm_opt$par

  #####
  y0.tt <- rep(-1,n_obs)
  y0.tt.imp <- imp.randKern_func(h = n_samples[1]^(-1/3)*sqrt(var(x0.t)), x = x0.t, y = y0.na, k = 1, delta = delta_0)
  y0.tt[delta_0==1] <- 1/p0.t[delta_0==1]*y0[delta_0==1]+(1-1/p0.t[delta_0==1])*y0.tt.imp[delta_0==1]
  y0.tt[delta_0==0] <- y0.tt.imp[delta_0==0]

  y1.tt <- rep(-1,n_obs)
  y1.tt.imp <- imp.randKern_func(h = n_samples[2]^(-1/3)*sqrt(var(x1.t)), x = x1.t, y = y1.na, k = 1, delta = delta_1)
  y1.tt[delta_1==1] <- 1/p1.t[delta_1==1]*y1[delta_1==1]+(1-1/p1.t[delta_1==1])*y1.tt.imp[delta_1==1]
  y1.tt[delta_1==0] <- y1.tt.imp[delta_1==0]
  
  y2.tt <- rep(-1,n_obs)
  y2.tt.imp <- imp.randKern_func(h = n_samples[3]^(-1/3)*sqrt(var(x2.t)), x = x2.t, y = y2.na, k = 1, delta = delta_2)
  y2.tt[delta_2==1] <- 1/p2.t[delta_2==1]*y2[delta_2==1]+(1-1/p2.t[delta_2==1])*y2.tt.imp[delta_2==1]
  y2.tt[delta_2==0] <- y2.tt.imp[delta_2==0]
  
  y3.tt <- rep(-1,n_obs)
  y3.tt.imp <- imp.randKern_func(h = n_samples[4]^(-1/3)*sqrt(var(x3.t)), x = x3.t, y = y3.na, k = 1, delta = delta_3)
  y3.tt[delta_3==1] <- 1/p3.t[delta_3==1]*y3[delta_3==1]+(1-1/p3.t[delta_3==1])*y3.tt.imp[delta_3==1]
  y3.tt[delta_3==0] <- y3.tt.imp[delta_3==0]
  
  y.tt <- c(y0.tt, y1.tt, y2.tt, y3.tt)
  drmfit.tt <- drmdel(x=y.tt, n_samples=n_samples, basis_func=5)
  
  
  #####
  y0.tf <- rep(-1,n_obs)
  y0.tf.imp <- imp.randKern_func(h = n_samples[1]^(-1/3)*sqrt(var(x0.t)), x = x0.t, y = y0.na, k = 1, delta = delta_0)
  y0.tf[delta_0==1] <- 1/p0.f[delta_0==1]*y0[delta_0==1]+(1-1/p0.f[delta_0==1])*y0.tf.imp[delta_0==1]
  y0.tf[delta_0==0] <- y0.tf.imp[delta_0==0]
  
  y1.tf <- rep(-1,n_obs)
  y1.tf.imp <- imp.randKern_func(h = n_samples[2]^(-1/3)*sqrt(var(x1.t)), x = x1.t, y = y1.na, k = 1, delta = delta_1)
  y1.tf[delta_1==1] <- 1/p1.f[delta_1==1]*y1[delta_1==1]+(1-1/p1.f[delta_1==1])*y1.tf.imp[delta_1==1]
  y1.tf[delta_1==0] <- y1.tf.imp[delta_1==0]
  
  y2.tf <- rep(-1,n_obs)
  y2.tf.imp <- imp.randKern_func(h = n_samples[3]^(-1/3)*sqrt(var(x2.t)), x = x2.t, y = y2.na, k = 1, delta = delta_2)
  y2.tf[delta_2==1] <- 1/p2.f[delta_2==1]*y2[delta_2==1]+(1-1/p2.f[delta_2==1])*y2.tf.imp[delta_2==1]
  y2.tf[delta_2==0] <- y2.tf.imp[delta_2==0]
  
  y3.tf <- rep(-1,n_obs)
  y3.tf.imp <- imp.randKern_func(h = n_samples[4]^(-1/3)*sqrt(var(x3.t)), x = x3.t, y = y3.na, k = 1, delta = delta_3)
  y3.tf[delta_3==1] <- 1/p3.f[delta_3==1]*y3[delta_3==1]+(1-1/p3.f[delta_3==1])*y3.tf.imp[delta_3==1]
  y3.tf[delta_3==0] <- y3.tf.imp[delta_3==0]
  
  y.tf <- c(y0.tf, y1.tf, y2.tf, y3.tf)
  drmfit.tf <- drmdel(x=y.tf, n_samples=n_samples, basis_func=5)
  
  
  #####
  y0.ft <- rep(-1,n_obs)
  y0.ft.imp <- imp.randKern_func(h = n_samples[1]^(-1/3)*sqrt(var(x0.f)), x = x0.f, y = y0.na, k = 1, delta = delta_0)
  y0.ft[delta_0==1] <- 1/p0.t[delta_0==1]*y0[delta_0==1]+(1-1/p0.t[delta_0==1])*y0.ft.imp[delta_0==1]
  y0.ft[delta_0==0] <- y0.ft.imp[delta_0==0]
  
  y1.ft <- rep(-1,n_obs)
  y1.ft.imp <- imp.randKern_func(h = n_samples[2]^(-1/3)*sqrt(var(x1.f)), x = x1.f, y = y1.na, k = 1, delta = delta_1)
  y1.ft[delta_1==1] <- 1/p1.t[delta_1==1]*y1[delta_1==1]+(1-1/p1.t[delta_1==1])*y1.ft.imp[delta_1==1]
  y1.ft[delta_1==0] <- y1.ft.imp[delta_1==0]
  
  y2.ft <- rep(-1,n_obs)
  y2.ft.imp <- imp.randKern_func(h = n_samples[3]^(-1/3)*sqrt(var(x2.f)), x = x2.f, y = y2.na, k = 1, delta = delta_2)
  y2.ft[delta_2==1] <- 1/p2.t[delta_2==1]*y2[delta_2==1]+(1-1/p2.t[delta_2==1])*y2.ft.imp[delta_2==1]
  y2.ft[delta_2==0] <- y2.ft.imp[delta_2==0]
  
  y3.ft <- rep(-1,n_obs)
  y3.ft.imp <- imp.randKern_func(h = n_samples[4]^(-1/3)*sqrt(var(x3.f)), x = x3.f, y = y3.na, k = 1, delta = delta_3)
  y3.ft[delta_3==1] <- 1/p3.t[delta_3==1]*y3[delta_3==1]+(1-1/p3.t[delta_3==1])*y3.ft.imp[delta_3==1]
  y3.ft[delta_3==0] <- y3.ft.imp[delta_3==0]
  
  y.ft <- c(y0.ft, y1.ft, y2.ft, y3.ft)
  drmfit.ft <- drmdel(x=y.ft, n_samples=n_samples, basis_func=5)
  
  
  #####
  y0.ff <- rep(-1,n_obs)
  y0.ff.imp <- imp.randKern_func(h = n_samples[1]^(-1/3)*sqrt(var(x0.f)), x = x0.f, y = y0.na, k = 1, delta = delta_0)
  y0.ff[delta_0==1] <- 1/p0.f[delta_0==1]*y0[delta_0==1]+(1-1/p0.f[delta_0==1])*y0.ff.imp[delta_0==1]
  y0.ff[delta_0==0] <- y0.ff.imp[delta_0==0]
  
  y1.ff <- rep(-1,n_obs)
  y1.ff.imp <- imp.randKern_func(h = n_samples[2]^(-1/3)*sqrt(var(x1.f)), x = x1.f, y = y1.na, k = 1, delta = delta_1)
  y1.ff[delta_1==1] <- 1/p1.f[delta_1==1]*y1[delta_1==1]+(1-1/p1.f[delta_1==1])*y1.ff.imp[delta_1==1]
  y1.ff[delta_1==0] <- y1.ff.imp[delta_1==0]
  
  y2.ff <- rep(-1,n_obs)
  y2.ff.imp <- imp.randKern_func(h = n_samples[3]^(-1/3)*sqrt(var(x2.f)), x = x2.f, y = y2.na, k = 1, delta = delta_2)
  y2.ff[delta_2==1] <- 1/p2.f[delta_2==1]*y2[delta_2==1]+(1-1/p2.f[delta_2==1])*y2.ff.imp[delta_2==1]
  y2.ff[delta_2==0] <- y2.ff.imp[delta_2==0]
  
  y3.ff <- rep(-1,n_obs)
  y3.ff.imp <- imp.randKern_func(h = n_samples[4]^(-1/3)*sqrt(var(x3.f)), x = x3.f, y = y3.na, k = 1, delta = delta_3)
  y3.ff[delta_3==1] <- 1/p3.f[delta_3==1]*y3[delta_3==1]+(1-1/p3.f[delta_3==1])*y3.ff.imp[delta_3==1]
  y3.ff[delta_3==0] <- y3.ff.imp[delta_3==0]
  
  y.ff <- c(y0.ff, y1.ff, y2.ff, y3.ff)
  drmfit.ff <- drmdel(x=y.ff, n_samples=n_samples, basis_func=5)
  
  
  output <- rbind(c(mean(y0.obs),mean(y1.obs),mean(y2.obs),mean(y3.obs)),
                  c(mean(y0),mean(y1),mean(y2),mean(y3)),
                  c(mean(y0.k),mean(y1.k),mean(y2.k),mean(y3.k)),
                  c(mean(y0.tt),mean(y1.tt),mean(y2.tt),mean(y3.tt)),
                  c(mean(y0.tf),mean(y1.tf),mean(y2.tf),mean(y3.tf)),
                  c(mean(y0.ft),mean(y1.ft),mean(y2.ft),mean(y3.ft)),
                  c(mean(y0.ff),mean(y1.ff),mean(y2.ff),mean(y3.ff)),
                  c(sum(drmfit.cc$p_est$x[which(drmfit.cc$p_est$k==0)]*drmfit.cc$p_est$p_est[which(drmfit.cc$p_est$k==0)]), sum(drmfit.cc$p_est$x[which(drmfit.cc$p_est$k==1)]*drmfit.cc$p_est$p_est[which(drmfit.cc$p_est$k==1)]), sum(drmfit.cc$p_est$x[which(drmfit.cc$p_est$k==2)]*drmfit.cc$p_est$p_est[which(drmfit.cc$p_est$k==2)]), sum(drmfit.cc$p_est$x[which(drmfit.cc$p_est$k==3)]*drmfit.cc$p_est$p_est[which(drmfit.cc$p_est$k==3)])),
                  c(sum(drmfit.full$p_est$x[which(drmfit.full$p_est$k==0)]*drmfit.full$p_est$p_est[which(drmfit.full$p_est$k==0)]), sum(drmfit.full$p_est$x[which(drmfit.full$p_est$k==1)]*drmfit.full$p_est$p_est[which(drmfit.full$p_est$k==1)]), sum(drmfit.full$p_est$x[which(drmfit.full$p_est$k==2)]*drmfit.full$p_est$p_est[which(drmfit.full$p_est$k==2)]), sum(drmfit.full$p_est$x[which(drmfit.full$p_est$k==3)]*drmfit.full$p_est$p_est[which(drmfit.full$p_est$k==3)])),
                  c(sum(drmfit.k$p_est$x[which(drmfit.k$p_est$k==0)]*drmfit.k$p_est$p_est[which(drmfit.k$p_est$k==0)]), sum(drmfit.k$p_est$x[which(drmfit.k$p_est$k==1)]*drmfit.k$p_est$p_est[which(drmfit.k$p_est$k==1)]), sum(drmfit.k$p_est$x[which(drmfit.k$p_est$k==2)]*drmfit.k$p_est$p_est[which(drmfit.k$p_est$k==2)]), sum(drmfit.k$p_est$x[which(drmfit.k$p_est$k==3)]*drmfit.k$p_est$p_est[which(drmfit.k$p_est$k==3)])),
                  c(sum(drmfit.tt$p_est$x[which(drmfit.tt$p_est$k==0)]*drmfit.tt$p_est$p_est[which(drmfit.tt$p_est$k==0)]), sum(drmfit.tt$p_est$x[which(drmfit.tt$p_est$k==1)]*drmfit.tt$p_est$p_est[which(drmfit.tt$p_est$k==1)]), sum(drmfit.tt$p_est$x[which(drmfit.tt$p_est$k==2)]*drmfit.tt$p_est$p_est[which(drmfit.tt$p_est$k==2)]), sum(drmfit.tt$p_est$x[which(drmfit.tt$p_est$k==3)]*drmfit.tt$p_est$p_est[which(drmfit.tt$p_est$k==3)])),
                  c(sum(drmfit.tf$p_est$x[which(drmfit.tf$p_est$k==0)]*drmfit.tf$p_est$p_est[which(drmfit.tf$p_est$k==0)]), sum(drmfit.tf$p_est$x[which(drmfit.tf$p_est$k==1)]*drmfit.tf$p_est$p_est[which(drmfit.tf$p_est$k==1)]), sum(drmfit.tf$p_est$x[which(drmfit.tf$p_est$k==2)]*drmfit.tf$p_est$p_est[which(drmfit.tf$p_est$k==2)]), sum(drmfit.tf$p_est$x[which(drmfit.tf$p_est$k==3)]*drmfit.tf$p_est$p_est[which(drmfit.tf$p_est$k==3)])),
                  c(sum(drmfit.ft$p_est$x[which(drmfit.ft$p_est$k==0)]*drmfit.ft$p_est$p_est[which(drmfit.ft$p_est$k==0)]), sum(drmfit.ft$p_est$x[which(drmfit.ft$p_est$k==1)]*drmfit.ft$p_est$p_est[which(drmfit.ft$p_est$k==1)]), sum(drmfit.ft$p_est$x[which(drmfit.ft$p_est$k==2)]*drmfit.ft$p_est$p_est[which(drmfit.ft$p_est$k==2)]), sum(drmfit.ft$p_est$x[which(drmfit.ft$p_est$k==3)]*drmfit.ft$p_est$p_est[which(drmfit.ft$p_est$k==3)])),
                  c(sum(drmfit.ff$p_est$x[which(drmfit.ff$p_est$k==0)]*drmfit.ff$p_est$p_est[which(drmfit.ff$p_est$k==0)]), sum(drmfit.ff$p_est$x[which(drmfit.ff$p_est$k==1)]*drmfit.ff$p_est$p_est[which(drmfit.ff$p_est$k==1)]), sum(drmfit.ff$p_est$x[which(drmfit.ff$p_est$k==2)]*drmfit.ff$p_est$p_est[which(drmfit.ff$p_est$k==2)]), sum(drmfit.ff$p_est$x[which(drmfit.ff$p_est$k==3)]*drmfit.ff$p_est$p_est[which(drmfit.ff$p_est$k==3)]))
  )
  
}

stopCluster(cl)

Mean <- Reduce("+",results)/5000

bias.func <- function(x){
  (x-matrix(theta, byrow = TRUE, ncol = 4, nrow = 14))
}
rb.func <- function(x){
  (x-matrix(theta, byrow = TRUE, ncol = 4, nrow = 14))/matrix(theta, byrow = TRUE, ncol = 4, nrow = 14)
}
var.func <- function(x){
  (x-Mean)^2
}
mse.func <- function(x){
  (x-matrix(theta, byrow = TRUE, ncol = 4, nrow = 14))^2
}


Bias <- Reduce("+",lapply(results, bias.func))/5000
RB <- Reduce("+",lapply(results, rb.func))/5000*100
MSE <- Reduce("+",lapply(results, mse.func))/5000
Var <- Reduce("+",lapply(results, var.func))/5000

RB <- as.data.frame(RB,
                    row.names = c("sample.cc","sample.full","em.k","em.tt","em.tf","em.ft","em.ff",
                                  "drm.cc","drm.full","drm.k","drm.tt","drm.tf","drm.ft","drm.ff"), 
                    col.names = c("0","1","2","3"))
MSE <- as.data.frame(MSE,
                     row.names = c("sample.cc","sample.full","em.k","em.tt","em.tf","em.ft","em.ff",
                                   "drm.cc","drm.full","drm.k","drm.tt","drm.tf","drm.ft","drm.ff"), 
                     col.names = c("0","1","2","3"))
Bias <- as.data.frame(Bias,
                      row.names = c("sample.cc","sample.full","em.k","em.tt","em.tf","em.ft","em.ff",
                                    "drm.cc","drm.full","drm.k","drm.tt","drm.tf","drm.ft","drm.ff"), 
                      col.names = c("0","1","2","3"))
Var <- as.data.frame(Var,
                     row.names = c("sample.cc","sample.full","em.k","em.tt","em.tf","em.ft","em.ff",
                                   "drm.cc","drm.full","drm.k","drm.tt","drm.tf","drm.ft","drm.ff"), 
                     col.names = c("0","1","2","3"))


write.table(RB, paste("rb_x=",n_obs,".txt",sep = ""))
write.table(MSE, paste("mse_x=",n_obs,".txt",sep = ""))
write.table(Bias, paste("bias_x=",n_obs,".txt",sep = ""))
write.table(Var, paste("var_x=",n_obs,".txt",sep = ""))
saveRDS(results,file = paste("mu_",n_obs,".RData",sep = ""))

}

print("Done!")
warnings()
