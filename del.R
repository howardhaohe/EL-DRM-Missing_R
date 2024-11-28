# 2024-11-25
#


#### Classic DEL-EL formulation ####

negLDELfn <- function(par, x, n_total, n_samples, m, basis_func, d) {
  # Calculate negative log dual empirical likelihood for a
  # given value of parameter.
  # par: A vector of parameter vectors theta_r for each population r (2 to m+1)
  # x: A vector of all samples from all populations (numeric vector)
  # basis_func: A function to compute Q(x) given an observation x
  # n_samples: A vector of sample sizes for each population
  # n_total: Total sample size (sum of all n_samples)
  # m+1: number samples
  # d+1: dimension of basis function
  
  
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (sum(n_samples) != n_total) {
    stop("total sample size is not consistent with 'n_samples'!")
  }  # check n_samples lengths
  
  if (length(par) != m*(d+1)) {
    stop("parameters length is not m*(d+1)!")
  }  # check parameter length
  
  
  rho <- n_samples / n_total    # Calculate rho values for each population
  theta <- c(rep(0,d+1),par)    # add theta_0 to parameters
  Qx <- t(sapply(x, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  
  
  if (dim(Qx)[2] != m) { 
    stop("dimension of Q(x_kj) should be n*m!")
  }   # check Q(x_kj)
  
  
  # Split x into a list of samples based on n_samples
  x_list <- split(x = x, f = rep(1:(m+1), n_samples))
  # Split parameter into a list of parameters based on n_samples
  theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  
  # h(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  h_r_kj <- lapply(1:(m+1), function(r) {
    rho[r] * exp(Qx %*% theta_list[[r]])
  })
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj <- Reduce("+",h_r_kj)
  
  # linear term: A list of length m+1, each element is a length n_total vector with each row to be theta_r^T*Q(x_kj)
  tail_term <- lapply(1:(m+1), function(k) {
    t(sapply(x_list[[k]], basis_func)) %*% (theta_list[[k]])
  })
  
  
  # summation over all k,j
  # returns the negative value of DEL evaluate at a specific theta
  return(sum(log(s_kj)-do.call(rbind,tail_term)))
  
}


negLDELgr <- function(par, x, n_total, n_samples, m, basis_func, d) {
  # Calculate gradient of negative log dual empirical likelihood for a
  # given value of parameter.
  # par: A vector of parameter vectors theta_r for each population r (2 to m+1)
  # x: A vector of all samples from all populations (numeric vector)
  # basis_func: A function to compute Q(x) given an observation x
  # n_samples: A vector of sample sizes for each population
  # n_total: Total sample size (sum of all n_samples)
  # m+1: number samples
  # d+1: dimension of basis function
  
  
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (sum(n_samples) != n_total) {
    stop("total sample size is not consistent with 'n_samples'!")
  }  # check n_samples lengths
  
  if (length(par) != m*(d+1)) {
    stop("parameters length is not m*(d+1)!")
  }  # check parameter length
  
  
  rho <- n_samples / n_total    # Calculate rho values for each population
  theta <- c(rep(0,d+1),par)    # add theta_0 to parameters
  Qx <- t(sapply(x, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  
  
  if (dim(Qx)[2] != m) { 
    stop("dimension of Q(x_kj) should be n*m!")
  }   # check Q(x_kj)
  
  # Split x into a list of samples based on n_samples
  x_list <- split(x = x, f = rep(1:(m+1), n_samples))
  
  # Split parameter into a list of parameters based on n_samples
  theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  
  # Split Q(x_kj) into a list based on n_samples
  Qx_list <- lapply(1:(m+1), function(k) {
    t(sapply(x_list[[k]], basis_func))
  })
  
  
  # h(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  h_r_kj <- lapply(1:(m+1), function(r) {
    rho[r] * exp(Qx %*% theta_list[[r]])
  })
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj <- Reduce("+",h_r_kj)
  
  # Find start and end indices for each samples
  start_indices <- cumsum(c(1, n_samples[-length(n_samples)])) # Start indices for each segment
  end_indices <- cumsum(n_samples) # End indices for each segment
  
  # e_k as a list of length m+1
  e_k <- lapply(seq_along(n_samples), function(k) {
    .vec <- rep(0, sum(n_samples))  # Initialize a vector of zeros of the required length
    .vec[start_indices[k]:end_indices[k]] <- 1  # Set elements from start to end index to 1
    return(.vec)
  })
  
  # parLk_alpha: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. alpha
  parLk_alpha <- mapply(function(mat1,mat2,e_k){as.vector(mat1/mat2)-e_k}, mat1 = h_r_kj, e_k = e_k, MoreArgs = list(mat2 = s_kj), SIMPLIFY = FALSE)
  
  # parLk_theta: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. theta
  parLk_theta <-  mapply(function(mat1,mat2){mat1*mat2}, mat1 = parLk_alpha, MoreArgs = list(mat2 = Qx), SIMPLIFY = FALSE)
  
  
  # summation over all k,j
  # returns a m*(d+1)-vector of negative value of the partial derivatives of DEL evaluate at a specific theta 
  return(c(t(do.call(rbind, lapply(parLk_theta, colSums))[-1,])))
  
}




#### DEL-EL with missing value Idea 4 ####

negLDELfn_miss4 <- function(par, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d) {
  # Calculate negative log dual empirical likelihood for a
  # given value of parameter.
  # par: A vector of parameter vectors theta_r for each population r (2 to m+1)
  # x: A vector of all samples from all populations (numeric vector)
  # x_imp: A vector of all imputed values from all populations (numeric vector)
  # delta: A vector of all response indicator from all populations (0,1)
  # px_hat: A vector of all estimated propensity score from all populations (numeric vector)
  # basis_func: A function to compute Q(x) given an observation x
  # n_samples: A vector of sample sizes for each population
  # n_total: Total sample size (sum of all n_samples)
  # m+1: number samples
  # d+1: dimension of basis function
  
  
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(x_imp) != n_total) {
    stop("Length of imputed data vector 'x_imp' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(delta) != n_total) {
    stop("Length of response indicator data vector 'delta' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (any(is.na(x[delta==1]))) {
    stop("There are 'NA's in the response!")
  }  # check data lengths
  
  if (length(px_hat) != n_total) {
    stop("Length of estimated response probability data vector 'px_hat' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (sum(n_samples) != n_total) {
    stop("total sample size is not consistent with 'n_samples'!")
  }  # check n_samples lengths
  
  if (length(par) != m*(d+1)) {
    stop("parameters length is not m*(d+1)!")
  }  # check parameter length
  
  
  rho <- n_samples / n_total    # Calculate rho values for each population
  theta <- c(rep(0,d+1),par)    # add theta_0 to parameters
  Qx <- t(sapply(x, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  Qx_imp <- t(sapply(x_imp, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  
  
  if (dim(Qx)[2] != m) { 
    stop("dimension of Q(x_kj) should be n*m!")
  }   # check Q(x_kj)
  
  
  # Split x and x_imp into a list of samples based on n_samples
  x_list <- split(x = x, f = rep(1:(m+1), n_samples))
  x_imp_list <- split(x = x_imp, f = rep(1:(m+1), n_samples))
  
  # Split parameter into a list of parameters based on n_samples
  theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  
  # h(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  h_r_kj <- lapply(1:(m+1), function(r) {
    rho[r] * exp(Qx %*% theta_list[[r]])
  })
  h_r_kj_imp <- lapply(1:(m+1), function(r) {
    rho[r] * exp(Qx_imp %*% theta_list[[r]])
  })
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj <- Reduce("+",h_r_kj)
  s_kj_imp <- Reduce("+",h_r_kj_imp)
  
  # linear term: A list of length m+1, each element is a length n_total vector with each row to be theta_r^T*Q(x_kj)
  tail_term <- lapply(1:(m+1), function(k) {
    t(sapply(x_list[[k]], basis_func)) %*% (theta_list[[k]])
  })
  tail_term_imp <- lapply(1:(m+1), function(k) {
    t(sapply(x_imp_list[[k]], basis_func)) %*% (theta_list[[k]])
  })
  
  # tilde L_kj
  L_kj <- 1/px_hat*(log(s_kj) - do.call(rbind,tail_term)) + (1-1/px_hat)*(log(s_kj_imp) - do.call(rbind,tail_term_imp))
  # L_kj <- (log(s_kj) - do.call(rbind,tail_term))
  L_kj[delta==0] <- (log(s_kj_imp) - do.call(rbind,tail_term_imp))[delta==0]
  
  # summation over all k,j
  # returns the negative value of DEL evaluate at a specific theta
  return(sum(L_kj))
  
}


negLDELgr_miss4 <- function(par, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d) {
  # Calculate negative log dual empirical likelihood for a
  # given value of parameter.
  # par: A vector of parameter vectors theta_r for each population r (2 to m+1)
  # x: A vector of all samples from all populations (numeric vector)
  # x_imp: A vector of all imputed values from all populations (numeric vector)
  # delta: A vector of all response indicator from all populations (0,1)
  # px_hat: A vector of all estimated propensity score from all populations (numeric vector)
  # basis_func: A function to compute Q(x) given an observation x
  # n_samples: A vector of sample sizes for each population
  # n_total: Total sample size (sum of all n_samples)
  # m+1: number samples
  # d+1: dimension of basis function
  
  
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(x_imp) != n_total) {
    stop("Length of imputed data vector 'x_imp' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(delta) != n_total) {
    stop("Length of response indicator data vector 'delta' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (any(is.na(x[delta==1]))) {
    stop("There are 'NA's in the response!")
  }  # check data lengths
  
  if (length(px_hat) != n_total) {
    stop("Length of estimated response probability data vector 'px_hat' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (sum(n_samples) != n_total) {
    stop("total sample size is not consistent with 'n_samples'!")
  }  # check n_samples lengths
  
  if (length(par) != m*(d+1)) {
    stop("parameters length is not m*(d+1)!")
  }  # check parameter length
  
  
  rho <- n_samples / n_total    # Calculate rho values for each population
  theta <- c(rep(0,d+1),par)    # add theta_0 to parameters
  Qx <- t(sapply(x, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  Qx_imp <- t(sapply(x_imp, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  
  
  if (dim(Qx)[2] != m) { 
    stop("dimension of Q(x_kj) should be n*m!")
  }   # check Q(x_kj)
  
  
  # Split x and x_imp into a list of samples based on n_samples
  x_list <- split(x = x, f = rep(1:(m+1), n_samples))
  x_imp_list <- split(x = x_imp, f = rep(1:(m+1), n_samples))
  
  # Split parameter into a list of parameters based on n_samples
  theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  
  # h(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  h_r_kj <- lapply(1:(m+1), function(r) {
    rho[r] * exp(Qx %*% theta_list[[r]])
  })
  h_r_kj_imp <- lapply(1:(m+1), function(r) {
    rho[r] * exp(Qx_imp %*% theta_list[[r]])
  })
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj <- Reduce("+",h_r_kj)
  s_kj_imp <- Reduce("+",h_r_kj_imp)
  
  # Find start and end indices for each samples
  start_indices <- cumsum(c(1, n_samples[-length(n_samples)])) # Start indices for each segment
  end_indices <- cumsum(n_samples) # End indices for each segment
  
  # e_k as a list of length m+1
  e_k <- lapply(seq_along(n_samples), function(k) {
    .vec <- rep(0, sum(n_samples))  # Initialize a vector of zeros of the required length
    .vec[start_indices[k]:end_indices[k]] <- 1  # Set elements from start to end index to 1
    return(.vec)
  })
  
  # parLk_alpha: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. alpha
  parLk_alpha <- mapply(function(mat1,mat2,e_k){as.vector(mat1/mat2)-e_k}, mat1 = h_r_kj, e_k = e_k, MoreArgs = list(mat2 = s_kj), SIMPLIFY = FALSE)
  parLk_alpha_imp <- mapply(function(mat1,mat2,e_k){as.vector(mat1/mat2)-e_k}, mat1 = h_r_kj_imp, e_k = e_k, MoreArgs = list(mat2 = s_kj_imp), SIMPLIFY = FALSE)
  
  # parLk_theta: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. theta
  parLk_theta <-  mapply(function(mat1,mat2){mat1*mat2}, mat1 = parLk_alpha, MoreArgs = list(mat2 = Qx), SIMPLIFY = FALSE)
  parLk_theta_imp <-  mapply(function(mat1,mat2){mat1*mat2}, mat1 = parLk_alpha_imp, MoreArgs = list(mat2 = Qx_imp), SIMPLIFY = FALSE)
  
  # partial tilde L_kj/partial theta
  parLk_kj <- lapply(1:(m+1), function(r) {
    .mat <- 1/px_hat * parLk_theta[[r]] + (1-1/px_hat) * parLk_theta_imp[[r]]
    .mat[delta==0,] <- parLk_theta_imp[[r]][delta==0,]
    return(.mat)
  })
  
  
  # summation over all k,j
  return(c(t(do.call(rbind, lapply(parLk_kj, colSums))[-1,])))  
  
}



#### Estimate p_kj and d F_r(y_kj) for r ####

p_hat_miss4 <- function(mele, x, x_imp, delta, n_total, n_samples, m, basis_func, d) {
## d F_r(\tilde{Y}_{kj}) = \delta_{kj}Y_{kj} + (1-\delta_{kj})\tilde{Y}_{kj}
  
  # Calculate mean estimator based on DRM.
  # mele: A vector of mele of theta_r for each population r (2 to m+1)
  # x: A vector of all samples from all populations (numeric vector)
  # x_imp: A vector of all imputed values from all populations (numeric vector)
  # delta: A vector of all response indicator from all populations (0,1)
  # basis_func: A function to compute Q(x) given an observation x
  # n_samples: A vector of sample sizes for each population
  # n_total: Total sample size (sum of all n_samples)
  # m+1: number samples
  # d+1: dimension of basis function
  
  
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(x_imp) != n_total) {
    stop("Length of imputed data vector 'x_imp' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(delta) != n_total) {
    stop("Length of response indicator data vector 'delta' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (any(is.na(x[delta==1]))) {
    stop("There are 'NA's in the response!")
  }  # check data lengths
  
  if (length(px_hat) != n_total) {
    stop("Length of estimated response probability data vector 'px_hat' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (sum(n_samples) != n_total) {
    stop("total sample size is not consistent with 'n_samples'!")
  }  # check n_samples lengths
  
  if (length(par) != m*(d+1)) {
    stop("parameters length is not m*(d+1)!")
  }  # check parameter length
  
  
  rho <- n_samples / n_total    # Calculate rho values for each population
  theta <- c(rep(0,d+1),par)    # add theta_0 to parameters
  Qx <- t(sapply(x, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  Qx_imp <- t(sapply(x_imp, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  
  
  if (dim(Qx)[2] != m) { 
    stop("dimension of Q(x_kj) should be n*m!")
  }   # check Q(x_kj)
  
  
  # Split x and x_imp into a list of samples based on n_samples
  x_list <- split(x = x, f = rep(1:(m+1), n_samples))
  x_imp_list <- split(x = x_imp, f = rep(1:(m+1), n_samples))
  
  # Split parameter into a list of parameters based on n_samples
  theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  
  # omega_r(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  omega_r_kj <- lapply(1:(m+1), function(r) {
    exp(Qx %*% theta_list[[r]])
  })
  omega_r_kj_imp <- lapply(1:(m+1), function(r) {
    exp(Qx_imp %*% theta_list[[r]])
  })
  
  h_r_kj <- mapply("*",omega_r_kj, rho, SIMPLIFY = FALSE)
  h_r_kj_imp <- mapply("*",omega_r_kj_imp, rho, SIMPLIFY = FALSE)
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj <- Reduce("+",h_r_kj)*n_total
  s_kj_imp <- Reduce("+",h_r_kj_imp)*n_total
  
  
  # eta_kj: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. alpha
  eta_kj <- mapply(function(mat1,mat2){as.vector(mat1/mat2)}, mat1 = omega_r_kj, MoreArgs = list(mat2 = s_kj), SIMPLIFY = FALSE)
  eta_kj_imp <- mapply(function(mat1,mat2){as.vector(mat1/mat2)}, mat1 = omega_r_kj_imp, MoreArgs = list(mat2 = s_kj_imp), SIMPLIFY = FALSE)
  
  # partial tilde L_kj/partial theta
  eta_tilde_kj <- lapply(1:(m+1), function(r) {
    .vec <- eta_kj[[r]]
    .vec[delta==0] <- eta_kj_imp[[r]][delta==0]
    return(.vec)
  })
  
  
  # summation over all k,j
  return(c(do.call(rbind, lapply(eta_kj, sum))))  
  
}



#### Estimate Population Parameters ####


mu_hat_miss4 <- function(mele, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d) {
  # Calculate d F_r(y_kj) estimator based on DRM.
  # mele: A vector of mele of theta_r for each population r (2 to m+1)
  # x: A vector of all samples from all populations (numeric vector)
  # x_imp: A vector of all imputed values from all populations (numeric vector)
  # delta: A vector of all response indicator from all populations (0,1)
  # px_hat: A vector of all estimated propensity score from all populations (numeric vector)
  # basis_func: A function to compute Q(x) given an observation x
  # n_samples: A vector of sample sizes for each population
  # n_total: Total sample size (sum of all n_samples)
  # m+1: number samples
  # d+1: dimension of basis function
  
  
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(x_imp) != n_total) {
    stop("Length of imputed data vector 'x_imp' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(delta) != n_total) {
    stop("Length of response indicator data vector 'delta' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (any(is.na(x[delta==1]))) {
    stop("There are 'NA's in the response!")
  }  # check data lengths
  
  if (length(px_hat) != n_total) {
    stop("Length of estimated response probability data vector 'px_hat' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (sum(n_samples) != n_total) {
    stop("total sample size is not consistent with 'n_samples'!")
  }  # check n_samples lengths
  
  if (length(par) != m*(d+1)) {
    stop("parameters length is not m*(d+1)!")
  }  # check parameter length
  
  
  rho <- n_samples / n_total    # Calculate rho values for each population
  theta <- c(rep(0,d+1),par)    # add theta_0 to parameters
  Qx <- t(sapply(x, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  Qx_imp <- t(sapply(x_imp, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  
  
  if (dim(Qx)[2] != m) { 
    stop("dimension of Q(x_kj) should be n*m!")
  }   # check Q(x_kj)
  
  
  # Split x and x_imp into a list of samples based on n_samples
  x_list <- split(x = x, f = rep(1:(m+1), n_samples))
  x_imp_list <- split(x = x_imp, f = rep(1:(m+1), n_samples))
  
  # Split parameter into a list of parameters based on n_samples
  theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  
  # omega_r(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  omega_r_kj <- lapply(1:(m+1), function(r) {
    exp(Qx %*% theta_list[[r]])
  })
  omega_r_kj_imp <- lapply(1:(m+1), function(r) {
    exp(Qx_imp %*% theta_list[[r]])
  })
  
  h_r_kj <- mapply("*",omega_r_kj, rho, SIMPLIFY = FALSE)
  h_r_kj_imp <- mapply("*",omega_r_kj_imp, rho, SIMPLIFY = FALSE)
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj <- Reduce("+",h_r_kj)*n_total
  s_kj_imp <- Reduce("+",h_r_kj_imp)*n_total
  
  
  # eta_kj: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. alpha
  eta_kj <- mapply(function(mat1,mat2){as.vector(mat1/mat2)*x}, mat1 = omega_r_kj, MoreArgs = list(mat2 = s_kj), SIMPLIFY = FALSE)
  eta_kj_imp <- mapply(function(mat1,mat2){as.vector(mat1/mat2)*x_imp}, mat1 = omega_r_kj_imp, MoreArgs = list(mat2 = s_kj_imp), SIMPLIFY = FALSE)
  
  # partial tilde L_kj/partial theta
  eta_tilde_kj <- lapply(1:(m+1), function(r) {
    .vec <- 1/px_hat * eta_kj[[r]] + (1-1/px_hat) * eta_kj_imp[[r]]
    .vec[delta==0] <- eta_kj_imp[[r]][delta==0]
    return(.vec)
  })
  
  
  # summation over all k,j
  return(c(do.call(rbind, lapply(eta_kj, sum))))  
  
}
