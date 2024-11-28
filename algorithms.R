# 2024-09-27
# This file includes algorithms for imputation and propensity score estimation

# require package "mvtnorm","dplyr"
require(mvtnorm)
# require(dplyr)




#### Estimating Propensity Score ####


Phat.kern_func <- function(x,delta,h = 1){
  # function to get the Phat using kernel method
  # x: data matrix (can be multivariate)
  # delta: response indicator
  # h: bandwidth (a value or a vector with same value of dim(x)[2])
  
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  n_sample <- dim(x.mat)[1]
  
  Phat <- rep(-1, n_sample)   # estimated propensity score
  
  if(length(h)==1){
    H_inv <- diag(1/h, nrow = x.dim, ncol = x.dim)
  } else if(length(h)==x.dim) {
    H_inv <- diag(1/h)
  } else {
    stop("different length for bandwidth and dimension of X")
  }
  
  for(i in 1:n_sample){
    # num.vec <- det(H_inv)^(-1)*dmvnorm(sweep(x.mat, 2, x.mat[i,]), sigma = diag(h^2,x.dim))
    num.vec <- dmvnorm(sweep(x.mat, 2, x.mat[i,])%*%H_inv)
    Phat[i] <- sum(delta*num.vec)/sum(num.vec)
  }
  
  return(Phat)
} # end function


Phat.score_func <- function(x,delta,c = 5){
  # function to get the Phat using score method
  # two methods for spliting data into different classes
  # x: data matrix (can be multivariate)
  # delta: response indicator
  # c: number of class (must be integer)
  
  
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  n_sample <- dim(x.mat)[1]
  df <- data.frame(delta,x)
  
  Phat <- rep(-1, n_sample)   # estimated propensity score
  
  # Logistic Regression
  log.mod <- glm(delta~., family = binomial(link = "logit"), data = df)
  
  # c groups
  # quantile
  row.order <- order(log.mod$fitted.values)
  for(g in 1:c){
    Phat[row.order[(1+(n_sample/c)*(g-1)):((n_sample/c)*g)]] <- sum(delta[row.order[(1+(n_sample/c)*(g-1)):((n_sample/c)*g)]])/(n_sample/c)
  }
  
  
  # # clustering
  # cl.obs <- kmeans(x = log.mod$fitted.values[delta==1], centers = c, iter.max = 1000, nstart = 50)      # based on observed units
  # 
  # # distance matrix to each center
  # dist.mat <- matrix(-1, n_sample, c)
  # for(j in 1:c){
  #   dist.mat[,j] <- sqrt((predict(log.mod,df)-cl.obs$centers[j])^2)
  # }
  # 
  # cl <- rep(-1, n_sample)    # clustering
  # for(i in 1:n_sample){
  #   cl[i] <- which(rank(dist.mat[i,]) == 1)
  # }
  # 
  # for(i in 1:n_sample){
  #   cl.id <- cl[i]
  #   Phat[i] <- sum(delta[which(cl==cl.id)])/length(delta[which(cl==cl.id)])
  # }
  
  return(Phat)
} # end function




#### Imputation ####


imp.randKern_func <- function(h = 1, k = 5, x, y, delta) {
  # random hot-deck imputation from an estimating condition distribution using kernel estimator
  
  # This algorithm imputes every Y (missing or observed).
  # The output is a n by k matrix with imputed values of Y.
  
  # h: bandwidth of the kernel estimator (default 1)
  # k: number imputations (default 5)
  # x: fully observed variables
  # y: missing variable
  # delta: response vector, where 0 means missing and 1 means observed
  # y and delta must be vectors, x can be a vector or a matrix
  
  
  
  # if (k%%1 != 0) {
  #   stop(print("k is not an integer."))
  # }
  # 
  # if (length(delta) != length(y)){
  #   stop(print("Y and delta have different length."))
  # }
  
  # if (length(delta) != dim(x)[1]){
  #   stop(print("The number of row for X is different from the length of delta."))
  # }
  
  
  
  n_sample <- length(y)       # sample size
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  
  # support of the imputation values 
  support.y <- y[delta==1]
  
  # Y_tilde_kernel is the matrix of the imputed value for all units (missing or not)
  Y_tilde <- matrix(NA, nrow = n_sample, ncol = k)
  
  # inverse of bandwidth matrix
  H_inv <- diag(1/h, nrow = dim(as.matrix(x))[2], ncol = dim(as.matrix(x))[2])
  
  # Estimated conditional distribution F_hat of Y_i = y_l given X_i(xi)
  # [F_hat]_{ij} = P(Y_i = Y_j|X_i)
  
  # get_Fhat <- function(xi){dmvnorm(sweep(x.mat,2,xi) %*% H_inv) / sum(delta * dmvnorm(sweep(x.mat,2,xi) %*% H_inv))}
  # F_hat <- as.matrix(apply(x.mat, 1, get_Fhat))
  # F_hat <- t(F_hat[delta == 1, ])
  # 
  # if(any(is.na(F_hat))){
  #   print("Can not sample from a kernel distribution")
  # }
  # stop
  # 
  # # Sampling function
  # sample_fun <- function(prob)(
  #   return(sample(support.y, n_sample=k, replace=TRUE, prob=prob))
  # )
  # 
  # if(k==1){
  #   # Sample from the discrete distribution
  #   Y_tilde <- as.numeric(apply(F_hat, 1, sample_fun))
  # } else {
  #   # Sample from the discrete distribution
  #   Y_tilde <- t(apply(F_hat, 1, sample_fun))
  # }
  # Y_tilde <- as.matrix(Y_tilde)
  
  for(i in 1:n_sample){
    # Estimated conditional distribution F_hat of Y_i = y_l given X_i
    F_hat <- dmvnorm(sweep(x.mat,2,x.mat[i,]) %*% H_inv) / sum(delta * dmvnorm(sweep(x.mat,2,x.mat[i,]) %*% H_inv))
    F_hat <- F_hat[delta == 1]
    
    if(any(is.na(F_hat))){
      stop("Can not sample from a kernel distribution")
    }
    
    Y_tilde[i,] <- sample(support.y, size = k, replace = TRUE, prob = F_hat)
  }
  
  return(Y_tilde)
} # end function


imp.detKern_func <- function(h = 1, x, y, delta) {
  # single deterministic kernel imputation
  
  # This algorithm imputes every Y (missing or observed).
  # The output is a n by 1 matrix with imputed values of Y.
  
  # h: bandwidth of the kernel estimator (default 1)
  # x: fully observed variables
  # y: missing variable
  # delta: response vector, where 0 means missing and 1 means observed
  # y and delta must be vectors, x can be a vector or a matrix
  
  
  
  # if (k%%1 != 0) {
  #   stop(print("k is not an integer."))
  # }
  # 
  # if (length(delta) != length(y)){
  #   stop(print("Y and delta have different length."))
  # }
  
  # if (length(delta) != dim(x)[1]){
  #   stop(print("The number of row for X is different from the length of delta."))
  # }
  
  
  
  n_sample <- length(y)       # sample size
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  
  # Y_tilde_kernel is the matrix of the imputed value for all units (missing or not)
  Y_tilde <- matrix(NA, nrow = n_sample, ncol = 1)
  
  weights.k <- getWeight.kern_func(x,delta,h)
  
  Y_tilde <- weights.k[,delta==1] %*% y[delta==1]
  
  return(Y_tilde)
} # end function


imp.score_func <- function(c = 10, k = 5, x, y, delta) {
  # random hot deck score method imputation.
  
  # This algorithm imputes every Y (missing or observed).
  # The output is a n by k matrix with imputed values of Y.
  
  # c: number of classes (default 10)
  # k: number imputations (default 5)
  # x: fully observed variables
  # y: missing variable
  # delta: response vector, where 0 means missing and 1 means observed
  # y and delta must be vectors, can be a vector or a matrix
  
  # The output is a n by k matrix with fully imputed values of Y.
  
  
  # if (k%%1 != 0) {
  #   stop(print("k is not an integer."))
  # }
  # 
  # if (c%%1 != 0) {
  #   stop(print("c is not an integer."))
  # }
  # 
  # if (length(delta) != length(y)){
  #   stop(print("Y and delta have different length."))
  # }
  # 
  # if (length(delta) != dim(x)[1]){
  #   stop(print("The number of row for X is different from the length of delta."))
  # }
  
  n_sample <- length(y)       # sample size
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  df <- data.frame(y,x)
  df.r <- df[delta==1,]
  g_levels <- 1:c
  
  # support of the imputation values 
  support.y <- y[delta==1]
  
  # Y_tilde_kernel is the matrix of the imputed value for all units (missing or not)
  Y_tilde <- matrix(NA, nrow = n_sample, ncol = k)
  
  # preliminary linear Regression
  lr.mod <- lm(y~., data = df.r)
  
  # class id
  g_id <- rep(-1, n_sample)
  
  # c groups
  
  # equal-quantile method (Eltinge & Yansaneh, 1997)
  row.order <- order(predict(lr.mod, df))
  # row.order <- order(x)
  
  for(g in g_levels){
    g_id[row.order[(1+(n_sample/length(g_levels))*(g-1)):((n_sample/length(g_levels))*g)]] <- g
  }
  
  # # classification algorithm (Seber, 1984, Chapter7)
  # cl.obs <- kmeans(x = lr.mod$fitted.values, centers = c, iter.max = 1000, nstart = 50)      # based on observed units
  # 
  # # distance matrix to each center
  # dist.mat <- matrix(-1, n_sample, c)
  # for(j in 1:c){
  #   dist.mat[,j] <- sqrt((predict(lr.mod,df)-cl.obs$centers[j])^2)
  # }
  # 
  # for(i in 1:n_sample){
  #   g_id[i] <- which(rank(dist.mat[i,]) == 1)
  # }
  for(g in g_levels[-c]){
    # print(sum(delta[which(g_id==g)]))
    if(sum(delta[which(g_id==g)]) == 0){
      g_id[which(g_id==g)] <- g+1
      g_levels <- g_levels[-which(g_levels==g)]
    }
  }

  for(g in rev(g_levels[-1])){
    # print(sum(delta[which(g_id==g)]))
    if(sum(delta[which(g_id==g)]) == 0){
      g_id[which(g_id==g)] <- g-1
      g_levels <- g_levels[-which(g_levels==g)]
    }
  }
  
  g_id.support <- g_id[which(delta==1)]
  
  for(i in 1:n_sample){
    temp_id <- g_id[i]
    support.y.temp_id <- support.y[which(g_id.support == temp_id)]
    Y_tilde[i,] <- sample(support.y.temp_id, size = k, replace = TRUE)
  }
  
  
  return(Y_tilde)
} # end function


imp.detLinear_func <- function(x, y, delta) {
  # deterministic linear regression imputation
  
  # This algorithm imputes every Y (missing or observed).
  # The output is a n by 1 matrix with imputed values of Y.
  
  # x: fully observed variables (x can be polynomial)
  # y: missing variable
  # delta: response vector, where 0 means missing and 1 means observed
  # y and delta must be vectors, can be a vector or a matrix
  
  
  if (length(delta) != length(y)){
    stop(print("Y and delta have different length."))
  }
  
  n_sample <- length(y)       # sample size
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  df <- data.frame(y,x)
  df.r <- df[delta==1,]
  
  # Linear Regression
  lr.mod <- lm(y~., data = df.r)
  
  # estimate Y_tilde through the linear regression model
  Y_tilde <- predict(lr.mod, newdata = df)
  
  
  return(Y_tilde)
  
} # end function


imp.randLinear_func <- function(k = 5, x, y, delta) {
  # random linear regression imputation
  
  # This algorithm imputes every Y (missing or observed).
  # The output is a n by 1 matrix with imputed values of Y.
  
  # x: fully observed variables (x can be polynomial)
  # y: missing variable
  # delta: response vector, where 0 means missing and 1 means observed
  # y and delta must be vectors, can be a vector or a matrix
  
  # The output is a n by k matrix with fully imputed values of Y.
  
  if (k%%1 != 0) {
    stop(print("k is not an integer."))
  }
  
  if (length(delta) != length(y)){
    stop(print("Y and delta have different length."))
  }
  
  n_sample <- length(y)       # sample size
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  # df <- data.frame(y,x)
  df <- data.frame(x)
  df.r <- data.frame(y,x)[delta==1,]
  
  # Y_tilde_kernel is the matrix of the imputed value for all units (missing or not)
  Y_tilde <- matrix(NA, nrow = n_sample, ncol = k)
  
  # Linear Regression
  lr.mod <- lm(y~., data = df.r)
  res.obs <- lr.mod$residuals
  
  # estimate Y_tilde through the linear regression model
  Y_hat <- predict(lr.mod, newdata = df)
  
  for(i in 1:n_sample){
    Y_tilde[i,] <- Y_hat[i] + sample(res.obs, k, replace=TRUE)
  }
  
  # if(k==1){
  #   Y_tilde <- as.numeric(sapply(Y_hat, function(x){x+sample(res, k, replace=TRUE)}))
  # } else {
  #   Y_tilde <- t(sapply(Y_hat, function(x){x+sample(res, k, replace=TRUE)}))
  #   colnames(Y_tilde) <- 1:k
  # }
  # Y_tilde <- as.matrix(Y_tilde)
  
  return(Y_tilde)
  
} # end function


imp.cluster_func <- function(c = 10, k = 5, x, y, delta) {
  # random hot deck clustering method imputation.
  
  # This algorithm imputes every Y (missing or observed).
  # The output is a n by k matrix with imputed values of Y.
  
  # c: number of classes (default 10)
  # k: number imputations (default 5)
  # x: fully observed variables
  # y: missing variable
  # delta: response vector, where 0 means missing and 1 means observed
  # y and delta must be vectors, can be a vector or a matrix
  
  # The output is a n by k matrix with fully imputed values of Y.
  
  
  # if (k%%1 != 0) {
  #   stop(print("k is not an integer."))
  # }
  # 
  # if (c%%1 != 0) {
  #   stop(print("c is not an integer."))
  # }
  # 
  # if (length(delta) != length(y)){
  #   stop(print("Y and delta have different length."))
  # }
  # 
  # if (length(delta) != dim(x)[1]){
  #   stop(print("The number of row for X is different from the length of delta."))
  # }
  
  n_sample <- length(y)       # sample size
  x.mat <- as.matrix(x)
  
  # standardize data
  x.std <- scale(x.mat, center = TRUE, scale = apply(x.mat, 2, sd, na.rm = TRUE))
  x.std.obs <- x.std[delta==1,]
  
  # support of the imputation values 
  support.y <- y[delta==1]
  
  # Y_tilde_kernel is the matrix of the imputed value for all units (missing or not)
  Y_tilde <- matrix(NA, nrow = n_sample, ncol = k)
  
  # class id
  cl <- rep(-1, n_sample)
  
  # clustering
  cl.obs <- kmeans(x = x.std.obs, centers = c, iter.max = 1000, nstart = 50)      # based on observed units
  
  # distance matrix to each center
  dist.mat <- matrix(-1, n_sample, c)
  for(j in 1:c){
    dist.mat[,j] <- sqrt(rowSums(sweep(x.std,2,cl.obs$centers[j,])^2))
  }
  
  for(i in 1:n_sample){
    cl[i] <- which(rank(dist.mat[i,]) == 1)
  }
  
  cl.support <- cl[which(delta==1)]
  
  for(i in 1:n_sample){
    cl.id <- cl[i]
    support.y.id <- support.y[which(cl.support == cl.id)]
    Y_tilde[i,] <- sample(support.y.id, size = k, replace = TRUE)
  }
  
  
  return(Y_tilde)
} # end function



imp.randKern.dr_func <- function(h = 1, k = 5, x, y, delta, p_hat) {
  # random hot-deck imputation with doubly robust from an 
  # estimating condition distribution using kernel estimator
  
  # This algorithm imputes every Y (missing or observed).
  # The output is a n by k matrix with imputed values of Y.
  
  # h: bandwidth of the kernel estimator (default 1)
  # k: number imputations (default 5)
  # x: fully observed variables
  # y: missing variable
  # delta: response vector, where 0 means missing and 1 means observed
  # y and delta must be vectors, x can be a vector or a matrix
  
  
  
  # if (k%%1 != 0) {
  #   stop(print("k is not an integer."))
  # }
  # 
  # if (length(delta) != length(y)){
  #   stop(print("Y and delta have different length."))
  # }
  
  # if (length(delta) != dim(x)[1]){
  #   stop(print("The number of row for X is different from the length of delta."))
  # }
  
  
  
  n_sample <- length(y)       # sample size
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  
  # support of the imputation values 
  support.y <- y[delta==1]
  
  # Y_tilde_kernel is the matrix of the imputed value for all units (missing or not)
  Y_tilde <- matrix(NA, nrow = n_sample, ncol = k)
  
  # inverse of bandwidth matrix
  H_inv <- diag(1/h, nrow = dim(as.matrix(x))[2], ncol = dim(as.matrix(x))[2])
  
  # Estimated conditional distribution F_hat of Y_i = y_l given X_i(xi)
  # [F_hat]_{ij} = P(Y_i = Y_j|X_i)
  
  # get_Fhat <- function(xi){dmvnorm(sweep(x.mat,2,xi) %*% H_inv) / sum(delta * dmvnorm(sweep(x.mat,2,xi) %*% H_inv))}
  # F_hat <- as.matrix(apply(x.mat, 1, get_Fhat))
  # F_hat <- t(F_hat[delta == 1, ])
  # 
  # if(any(is.na(F_hat))){
  #   print("Can not sample from a kernel distribution")
  # }
  # stop
  # 
  # # Sampling function
  # sample_fun <- function(prob)(
  #   return(sample(support.y, n_sample=k, replace=TRUE, prob=prob))
  # )
  # 
  # if(k==1){
  #   # Sample from the discrete distribution
  #   Y_tilde <- as.numeric(apply(F_hat, 1, sample_fun))
  # } else {
  #   # Sample from the discrete distribution
  #   Y_tilde <- t(apply(F_hat, 1, sample_fun))
  # }
  # Y_tilde <- as.matrix(Y_tilde)
  
  for(i in 1:n_sample){
    # Estimated conditional distribution F_hat of Y_i = y_l given X_i
    F_hat <- dmvnorm(sweep(x.mat,2,x.mat[i,]) %*% H_inv)*(1-p_hat)/p_hat / sum(delta * dmvnorm(sweep(x.mat,2,x.mat[i,]) %*% H_inv) * (1-p_hat)/p_hat)
    F_hat <- F_hat[delta == 1]
    
    if(any(is.na(F_hat))){
      stop("Can not sample from a kernel distribution")
    }
    
    Y_tilde[i,] <- sample(support.y, size = k, replace = TRUE, prob = F_hat)
  }
  
  return(Y_tilde)
} # end function





#### Estimating g function in proposed method ####

# for proposed method, estimating the weight of mg function

getWeight.kern_func <- function(x,delta,h){
  # function to get the weights (W) of mg_hat using kernel method
  # x: data matrix (can be multivariate)
  # delta: response indicator
  # h: bandwidth (a value or a vector with same value of dim(x)[2])
  
  x.mat <- as.matrix(x)
  num.mat <- matrix(-1, n, n)   # numerator matrix
  x.dim <- dim(x.mat)[2]
  
  if(length(h)==1){
    H_inv <- diag(1/h, nrow = x.dim, ncol = x.dim)
  } else if(length(h)==x.dim) {
    H_inv <- diag(1/h)
  } else {
    stop("different length for bandwidth and dimension of X")
  }
  
  for(j in 1:n_sample){
    # num.mat[,j] <- delta[j]*dmvnorm(sweep(x.mat, 2, x[j,])%*%H_inv)
    num.mat[,j] <- delta[j]*det(H_inv)^(-1)*dmvnorm(sweep(x.mat, 2, x.mat[j,]), sigma = diag(h^2,x.dim))
  }
  
  weight.mat <- apply(num.mat, 2, function(y){y/rowSums(num.mat)})
  
  return(weight.mat)
} # end function


getWeight.score_func <- function(x,y,delta,c){
  # function to get the weights (W) of mg_hat using score method
  # x: data matrix (can be multivariate)
  # delta: response indicator
  # h: bandwidth (a value or a vector with same value of dim(x)[2])
  
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  n_sample <- dim(x.mat)[1]
  df <- data.frame(y,x)
  df.r <- df[delta==1,]
  
  # preliminary Logistic Regression
  lr.mod <- lm(y~., data = df.r)
  
  # class id
  cl <- rep(-1, n_sample)
  
  # c groups
  # quantile
  row.order <- order(predict(lr.mod, df))
  for(g in 1:c){
    cl[row.order[(1+(n_sample/c)*(g-1)):((n_sample/c)*g)]] <- g
  }
  
  # # clustering
  # cl.obs <- kmeans(x = lr.mod$fitted.values, centers = c, iter.max = 1000, nstart = 50)      # based on observed units
  # 
  # # distance matrix to each center
  # dist.mat <- matrix(-1, n_sample, c)
  # for(j in 1:c){
  #   dist.mat[,j] <- sqrt((predict(lr.mod,df)-cl.obs$centers[j])^2)
  # }
  # 
  # cl <- rep(-1, n_sample)    # clustering
  # for(i in 1:n_sample){
  #   cl[i] <- which(rank(dist.mat[i,]) == 1)
  # }
  
  dum.mat <- matrix(0, n, n)    # dummy variable matrix
  for(i in 1:n_sample){
    dum.mat[i, which(cl==cl[i])] <- 1
    dum.mat[i,delta==0] <- 0
  }
  
  weight.mat <- apply(dum.mat,2,function(y){y/rowSums(dum.mat)})
  
  return(weight.mat)
} # end function


getWeight.cluster_func <- function(x,delta,c){
  # function to get the weights (W) of mg_hat using clustering method
  # x: data matrix (can be multivariate)
  # delta: response indicator
  # h: bandwidth (a value or a vector with same value of dim(x)[2])
  
  x.mat <- as.matrix(x)
  
  # standardize data
  x.std <- scale(x.mat, center = TRUE, scale = apply(x.mat, 2, sd, na.rm = TRUE))
  x.std.obs <- x.std[delta==1,]
  
  n_sample <- dim(x.mat)[1]
  
  # clustering
  cl.obs <- kmeans(x = x.std.obs, centers = c, iter.max = 1000, nstart = 50)      # based on observed units
  
  # distance matrix to each center
  dist.mat <- matrix(-1, n_sample, c)
  for(j in 1:c){
    dist.mat[,j] <- sqrt(rowSums(sweep(x.std,2,cl.obs$centers[j,])^2))
  }
  
  cl <- rep(-1, n_sample)    # clustering
  for(i in 1:n_sample){
    cl[i] <- which(rank(dist.mat[i,]) == 1)
  }
  
  dum.mat <- matrix(0, n, n)    # dummy variable matrix
  for(i in 1:n_sample){
    dum.mat[i, which(cl==cl[i])] <- 1
    dum.mat[i,delta==0] <- 0
  }
  
  weight.mat <- apply(dum.mat,2,function(y){y/rowSums(dum.mat)})
  
  return(weight.mat)
} # end function




#### Random Imputation for EE ####


gimp.score_func <-  function(c = 10, k = 1, x, y, delta) {
  # random hot deck score method imputation for ee. 
  
  # This algorithm imputes every Y (missing or observed).
  # The output is a n by k matrix with the ID of imputed values of Y.
  
  # c: number of classes (default 10)
  # k: number imputations (default 5)
  # x: fully observed variables
  # y: missing variable
  # delta: response vector, where 0 means missing and 1 means observed
  # y and delta must be vectors, can be a vector or a matrix
  
  # The output is a n by k matrix with fully imputed values of Y.
  
  
  # if (k%%1 != 0) {
  #   stop(print("k is not an integer."))
  # }
  # 
  # if (c%%1 != 0) {
  #   stop(print("c is not an integer."))
  # }
  # 
  # if (length(delta) != length(y)){
  #   stop(print("Y and delta have different length."))
  # }
  # 
  # if (length(delta) != dim(x)[1]){
  #   stop(print("The number of row for X is different from the length of delta."))
  # }
  
  n_sample <- length(y)       # sample size
  x.mat <- as.matrix(x)
  x.dim <- dim(x.mat)[2]
  df <- data.frame(y,x)
  dr.r <- dr[delta==1,]
  
  # support of the imputation values 
  support.id <- which(delta==1)
  
  # g_fun.id is the matrix of the id for g function (missing or not)
  g_fun.id <- matrix(NA, nrow = n_sample, ncol = k)
  
  # preliminary Linear Regression
  lr.mod <- lm(y~., data = df.r)
  
  # class id
  cl <- rep(-1, n_sample)
  
  # c groups
  # quantile
  row.order <- order(predict(lr.mod, df))
  for(g in 1:c){
    cl[row.order[(1+(n_sample/c)*(g-1)):((n_sample/c)*g)]] <- g
  }
  
  # # clustering
  # cl.obs <- kmeans(x = lr.mod$fitted.values, centers = c, iter.max = 1000, nstart = 50)      # based on observed units
  # 
  # # distance matrix to each center
  # dist.mat <- matrix(-1, n_sample, c)
  # for(j in 1:c){
  #   dist.mat[,j] <- sqrt((predict(lr.mod,df)-cl.obs$centers[j])^2)
  # }
  # 
  # for(i in 1:n_sample){
  #   cl[i] <- which(rank(dist.mat[i,]) == 1)
  # }
  
  for(i in 1:n_sample){
    cl.id <- cl[i]
    support.id <- which((cl == cl.id)&(delta==1))
    g_fun.id[i,] <- sample(support.id, size = k, replace = TRUE)
  }
  
  
  return(g_fun.id)
} # end function
