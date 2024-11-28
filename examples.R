# 2024-11-25
# 4 samples with normal conditioning on normal populations
# basis function 5
# estimating mean
# Doubly robust kernel imputation
# New formulation of DEL, idea 4, hh's algorithm


rm(list=ls())

source("algorithms.R")
source("del.R")

library(drmdel)

set.seed(100)

n_samples <- c(150,125,200,185)  # sample sizes
n_total <- sum(n_samples)

mu <- c(5*5+2, 4.5*4+1.5, 1.2 - 4.2*5.2, 4.5*5.5+0.8)
var <- c(1.3^2*5^2+8^2,1^2*4^2+9^2,1.5^2*4.2^2+8.5^2,1.8^2*4.5^2+8.7^2)

theta_func.n <- function(par_0,par_1,par_2,par_3){
  # model parameters
  return(c(log(sqrt(par_0[2]/par_1[2]))+(par_0[1]^2)/(2*par_0[2])-(par_1[1]^2)/(2*par_1[2]), par_1[1]/par_1[2]-par_0[1]/par_0[2], 1/(2*par_0[2])-1/(2*par_1[2]),
           log(sqrt(par_0[2]/par_2[2]))+(par_0[1]^2)/(2*par_0[2])-(par_2[1]^2)/(2*par_2[2]), par_2[1]/par_2[2]-par_0[1]/par_0[2], 1/(2*par_0[2])-1/(2*par_2[2]),
           log(sqrt(par_0[2]/par_3[2]))+(par_0[1]^2)/(2*par_0[2])-(par_3[1]^2)/(2*par_3[2]), par_3[1]/par_3[2]-par_0[1]/par_0[2], 1/(2*par_0[2])-1/(2*par_3[2])
           
  )
  )
}

theta <- theta_func.n(c(27,106.25),c(19.5,97),c(-20.64,111.94),c(25.55,141.30))

m=3
d=2


# Define the basis functions
basis_func <- function(x) {
  return(c(1,x, x^2))
}



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
  
  p0 <- 0.1+0.9*exp(x0.t-5)/(1+exp(x0.t-5))
  p0 <- ifelse(p0>1,1,p0)
  p1 <- 0.1+0.9*exp(x1.t-5)/(1+exp(x1.t-5))
  p1 <- ifelse(p1>1,1,p1)
  p2 <- 0.1+0.9*exp(x2.t-5)/(1+exp(x2.t-5))
  p2 <- ifelse(p2>1,1,p2)
  p3 <- 0.1+0.9*exp(x3.t-5)/(1+exp(x3.t-5))
  p3 <- ifelse(p3>1,1,p3)
  
  delta_0 <- rbinom(n_samples[1],1,p0)
  delta_1 <- rbinom(n_samples[2],1,p1)
  delta_2 <- rbinom(n_samples[3],1,p2)
  delta_3 <- rbinom(n_samples[4],1,p3)
  delta <- c(delta_0,delta_1,delta_2,delta_3)
  
  y0.obs <- y0
  y0.obs[delta_0==0] <- NA
  y1.obs <- y1
  y1.obs[delta_1==0] <- NA
  y2.obs <- y2
  y2.obs[delta_2==0] <- NA
  y3.obs <- y3
  y3.obs[delta_3==0] <- NA
  y.obs <- c(y0.obs,y1.obs,y2.obs,y3.obs)
  
  x0.f <- rnorm(n_samples[1], mean = 5, sd = 1.5)
  x1.f <- rnorm(n_samples[2], mean = 5, sd = 1.5)
  x2.f <- rnorm(n_samples[3], mean = 5, sd = 1.5)
  x3.f <- rnorm(n_samples[4], mean = 5, sd = 1.5)
  
  
  #### Estimate propensity score ####
  
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
  
  # False Propensity
  p0.f <- Phat.kern_func(x = x0.f, delta = delta_0, h = n_samples[1]^(-1/5)*sqrt(var(x0.f)))
  p1.f <- Phat.kern_func(x = x1.f, delta = delta_1, h = n_samples[1]^(-1/5)*sqrt(var(x1.f)))
  p2.f <- Phat.kern_func(x = x2.f, delta = delta_2, h = n_samples[1]^(-1/5)*sqrt(var(x2.f)))
  p3.f <- Phat.kern_func(x = x3.f, delta = delta_3, h = n_samples[1]^(-1/5)*sqrt(var(x3.f)))
  
  p.f <- c(p0.f,p1.f,p2.f,p3.f)
  
  # calculating MELE
  par_dim <- m*(d+1)
  par_init <- rep(0, par_dim)

  
  #### Full data estimator ####
  # drmdel
  drmfit.full <- drmdel(x=y, n_samples=n_samples, basis_func=5)
  
  # new del
  drm_opt.full <- do.call(optim, c(list(par=par_init, fn=negLDELfn,
                                      gr=negLDELgr,
                                      method="BFGS", control=list(maxit=10000),
                                      x=y,
                                      n_total=n_total, n_samples=n_samples, m=m,
                                      basis_func=basis_func,
                                      d=d)))
  mele.full <- drm_opt.full$par
  
  
  #### Naive estimator ####
  # drmdel
  n_samples.obs <- c(sum(delta_0==1), sum(delta_1==1), sum(delta_2==1), sum(delta_3==1))
  drmfit.cc <- drmdel(x = y.obs[delta==1], 
                      n_samples=n_samples.obs, basis_func=5)
  
  # new del
  drm_opt.cc <- do.call(optim, c(list(par=par_init, fn=negLDELfn,
                                        gr=negLDELgr,
                                        method="BFGS", control=list(maxit=10000),
                                        x=y.obs[delta==1],
                                        n_total=sum(n_samples.obs), n_samples=n_samples.obs, m=m,
                                        basis_func=basis_func,
                                        d=d)))
  mele.cc <- drm_opt.cc$par
  
  
  #### Random Kernel Regression Imputation #### 
  
  # imputation
  y0_hat_t <- imp.randKern_func(h = n_samples[1]^(-1/3)*sqrt(var(x0.t)), x = x0.t, y = y0.obs, k = 1, delta = delta_0)
  y1_hat_t <- imp.randKern_func(h = n_samples[2]^(-1/3)*sqrt(var(x1.t)), x = x1.t, y = y1.obs, k = 1, delta = delta_1)
  y2_hat_t <- imp.randKern_func(h = n_samples[3]^(-1/3)*sqrt(var(x2.t)), x = x2.t, y = y2.obs, k = 1, delta = delta_2)
  y3_hat_t <- imp.randKern_func(h = n_samples[4]^(-1/3)*sqrt(var(x3.t)), x = x3.t, y = y3.obs, k = 1, delta = delta_3)
  y_hat_t <- c(y0_hat_t, y1_hat_t, y2_hat_t, y3_hat_t)

  # drmdel
  y0.k <- y0.obs
  y0.k[delta_0==0] <- y0_hat_t[delta_0==0]
  y1.k <- y1.obs
  y1.k[delta_1==0] <- y1_hat_t[delta_1==0]
  y2.k <- y2.obs
  y2.k[delta_2==0] <- y2_hat_t[delta_2==0]
  y3.k <- y3.obs
  y3.k[delta_3==0] <- y3_hat_t[delta_3==0]
  y.k <- c(y0.k, y1.k, y2.k, y3.k)
  
  drmfit.k <- drmdel(x=y.k, n_samples=n_samples, basis_func=5)

  # new del
  drm_opt.k <- do.call(optim, c(list(par=par_init, fn=negLDELfn_miss4,
                                   gr=negLDELgr_miss4,
                                   method="BFGS", control=list(maxit=10000),
                                   x=y.obs, x_imp=y_hat_t, delta=delta, px_hat=rep(1,n_total),
                                   n_total=n_total, n_samples=n_samples, m=m,
                                   basis_func=basis_func,
                                   d=d)))
  mele.k <- drm_opt.k$par
  
  
  #### doubly robust ####

  # imputation
  y0_hat_tt <- imp.randKern.dr_func(h = n_samples[1]^(-1/3)*sqrt(var(x0.t)), x = x0.t, y = y0.obs, k = 1, delta = delta_0, p_hat = p0.t)
  y1_hat_tt <- imp.randKern.dr_func(h = n_samples[2]^(-1/3)*sqrt(var(x1.t)), x = x1.t, y = y1.obs, k = 1, delta = delta_1, p_hat = p1.t)
  y2_hat_tt <- imp.randKern.dr_func(h = n_samples[3]^(-1/3)*sqrt(var(x2.t)), x = x2.t, y = y2.obs, k = 1, delta = delta_2, p_hat = p2.t)
  y3_hat_tt <- imp.randKern.dr_func(h = n_samples[4]^(-1/3)*sqrt(var(x3.t)), x = x3.t, y = y3.obs, k = 1, delta = delta_3, p_hat = p3.t)
  y0.tt <- y0.obs
  y0.tt[delta_0==0] <- y0_hat_tt[delta_0==0]
  y1.tt <- y1.obs
  y1.tt[delta_1==0] <- y1_hat_tt[delta_1==0]
  y2.tt <- y2.obs
  y2.tt[delta_2==0] <- y2_hat_tt[delta_2==0]
  y3.tt <- y3.obs
  y3.tt[delta_3==0] <- y3_hat_tt[delta_3==0]
  y.tt <- c(y0.tt, y1.tt, y2.tt, y3.tt)
  
  # drmdel
  drmfit.tt <- drmdel(x=y.tt, n_samples=n_samples, basis_func=5)

  # new del
  drm_opt.tt <- do.call(optim, c(list(par=par_init, fn=negLDELfn_miss4,
                                     gr=negLDELgr_miss4,
                                     method="BFGS", control=list(maxit=10000),
                                     x=y.obs, x_imp=y_hat_t, delta=delta, px_hat=p.t,
                                     n_total=n_total, n_samples=n_samples, m=m,
                                     basis_func=basis_func,
                                     d=d)))
  mele.tt <- drm_opt.tt$par
  
  output <- rbind(drmfit.cc$mele,mele.cc,drmfit.full$mele,mele.full,drmfit.k$mele,mele.k,drmfit.tt$mele,mele.tt)
  

