# 2024-12-05
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




#### Estimate Population Parameters ####


# mu_hat <- function(mele, x, n_total, n_samples, m, basis_func, d) {
#   # Calculate d F_r(y_kj) estimator based on DRM.
#   # mele: A vector of mele of theta_r for each population r (2 to m+1)
#   # x: A vector of all samples from all populations (numeric vector)
#   # x_imp: A vector of all imputed values from all populations (numeric vector)
#   # delta: A vector of all response indicator from all populations (0,1)
#   # px_hat: A vector of all estimated propensity score from all populations (numeric vector)
#   # basis_func: A function to compute Q(x) given an observation x
#   # n_samples: A vector of sample sizes for each population
#   # n_total: Total sample size (sum of all n_samples)
#   # m+1: number samples
#   # d+1: dimension of basis function
#   
#   
#   if (length(x) != n_total) {
#     stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
#   }  # check data lengths
#   
#   if (sum(n_samples) != n_total) {
#     stop("total sample size is not consistent with 'n_samples'!")
#   }  # check n_samples lengths
#   
#   if (length(mele) != m*(d+1)) {
#     stop("parameters length is not m*(d+1)!")
#   }  # check parameter length
#   
#   
#   rho <- n_samples / n_total    # Calculate rho values for each population
#   theta <- c(rep(0,d+1),mele)    # add theta_0 to parameters
#   Qx <- t(sapply(x, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
#   
#   
#   if (dim(Qx)[2] != m) { 
#     stop("dimension of Q(x_kj) should be n*m!")
#   }   # check Q(x_kj)
#   
#   
#   # Split x and x_imp into a list of samples based on n_samples
#   x_list <- split(x = x, f = rep(1:(m+1), n_samples))
#   
#   # Split parameter into a list of parameters based on n_samples
#   theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
#   
#   # omega_r(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
#   omega_r_kj <- lapply(1:(m+1), function(r) {
#     exp(Qx %*% theta_list[[r]])
#   })
#   
#   h_r_kj <- mapply("*",omega_r_kj, rho, SIMPLIFY = FALSE)
#   
#   # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
#   s_kj <- Reduce("+",h_r_kj)*n_total
#   
#   
#   # eta_kj: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. alpha
#   eta_kj <- mapply(function(mat1,mat2){as.vector(mat1/mat2)*x}, mat1 = omega_r_kj, MoreArgs = list(mat2 = s_kj), SIMPLIFY = FALSE)
#   
#   # summation over all k,j
#   return(c(do.call(rbind, lapply(eta_kj, sum))))  
#   
# }



#### Idea 2: DEL-EL with missing value ####

negLDELfn_miss2 <- function(par, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d) {
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
  
  
  if(sum(delta) == n_total) {
    stop("No missing data. Use 'drmdel'!")
  }
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
  
  
  Qx_tilde <- 1/px_hat*Qx + (1-1/px_hat)*Qx_imp
  Qx_tilde[delta==0] <- Qx_imp[delta==0]
  
  # Ensure the sum of n_samples matches the number of rows in the Qx_tilde
  if (sum(n_samples) != nrow(Qx_tilde)) {
    stop("Sum of n_samples must equal the number of rows in the Qx_tilde")
  }
  
  # Use cumulative sums to define splitting indices
  split_indices <- cumsum(c(0, n_samples))
  
  # Split the Qx_tilde into a list of smaller matrices
  Qx_tilde_list <- lapply(seq_along(n_samples), function(i) {
    Qx_tilde[(split_indices[i] + 1):split_indices[i + 1], , drop = FALSE]
  })
  
  # Split parameter into a list of parameters based on n_samples
  theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  
  # h(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  h_r_kj_tilde <- lapply(1:(m+1), function(r) {
    rho[r] * exp(Qx_tilde %*% theta_list[[r]])
  })
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj_tilde <- Reduce("+",h_r_kj_tilde)
  
  # linear term: A list of length m+1, each element is a length n_total vector with each row to be theta_r^T*Q(x_kj)
  tail_term <- lapply(1:(m+1), function(k) {
    (Qx_tilde_list[[k]]) %*% (theta_list[[k]])
  })


  # summation over all k,j
  # returns the negative value of DEL evaluate at a specific theta
  return(sum(log(s_kj_tilde)-do.call(rbind,tail_term)))
  
}


negLDELgr_miss2 <- function(par, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d) {
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
  
  
  if(sum(delta) == n_total) {
    stop("No missing data. Use 'drmdel'!")
  }
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
  
  
  Qx_tilde <- 1/px_hat*Qx + (1-1/px_hat)*Qx_imp
  Qx_tilde[delta==0] <- Qx_imp[delta==0]
  
  # Ensure the sum of n_samples matches the number of rows in the Qx_tilde
  if (sum(n_samples) != nrow(Qx_tilde)) {
    stop("Sum of n_samples must equal the number of rows in the Qx_tilde")
  }
  
  # Use cumulative sums to define splitting indices
  split_indices <- cumsum(c(0, n_samples))
  
  # Split the Qx_tilde into a list of smaller matrices
  Qx_list <- lapply(seq_along(n_samples), function(i) {
    Qx_tilde[(split_indices[i] + 1):split_indices[i + 1], , drop = FALSE]
  })
  
  # Split parameter into a list of parameters based on n_samples
  theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  
  # tilde h(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  h_r_kj_tilde <- lapply(1:(m+1), function(r) {
    rho[r] * exp(Qx_tilde %*% theta_list[[r]])
  })
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj_tilde <- Reduce("+",h_r_kj)
  
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
  parLk_alpha <- mapply(function(mat1,mat2,e_k){as.vector(mat1/mat2)-e_k}, mat1 = h_r_kj_tilde, e_k = e_k, MoreArgs = list(mat2 = s_kj_tilde), SIMPLIFY = FALSE)

  # parLk_theta: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. theta
  parLk_theta <-  mapply(function(mat1,mat2){mat1*mat2}, mat1 = parLk_alpha, MoreArgs = list(mat2 = Qx_tilde), SIMPLIFY = FALSE)
  
  # summation over all k,j
  return(c(t(do.call(rbind, lapply(parLk_theta, colSums))[-1,])))  
  
}



#### Idea 3: DEL-EL with missing value ####

negLDELfn_miss3 <- function(par, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d) {
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
  
  
  if(sum(delta) == n_total) {
    stop("No missing data. Use 'drmdel'!")
  }
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
  
  # omega(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  omega_kj <- lapply(1:(m+1), function(r) {
    exp(Qx %*% theta_list[[r]])
  })
  omega_kj_imp <- lapply(1:(m+1), function(r) {
    exp(Qx_imp %*% theta_list[[r]])
  })
  
  omega_kj_tilde <- lapply(1:(m+1), function(r) {
    .vec <- 1/px_hat*omega_kj[[r]] + (1-1/px_hat)*omega_kj_imp[[r]]
    .vec[delta==0,] <- omega_kj_imp[[r]][delta==0,]
    return(.vec)
  })
    
  # h(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  h_r_kj_tilde <- lapply(1:(m+1), function(r) {
    rho[r] * omega_kj_tilde[[r]]
  })

  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj_tilde <- Reduce("+",h_r_kj_tilde)
  
  # linear term: A list of length m+1, each element is a length n_total vector with each row to be theta_r^T*Q(x_kj)
  tail_term <- lapply(1:(m+1), function(k) {
    t(sapply(x_list[[k]], basis_func)) %*% (theta_list[[k]])
  })
  tail_term_imp <- lapply(1:(m+1), function(k) {
    t(sapply(x_imp_list[[k]], basis_func)) %*% (theta_list[[k]])
  })

  
  # tilde L_kj
  L_kj <- 1/px_hat*(- do.call(rbind,tail_term)) + (1-1/px_hat)*(- do.call(rbind,tail_term_imp))
  # L_kj <- (log(s_kj) - do.call(rbind,tail_term))
  L_kj[delta==0] <- (- do.call(rbind,tail_term_imp))[delta==0]
  L_kj <- log(s_kj_tilde) + L_kj
  
  # summation over all k,j
  # returns the negative value of DEL evaluate at a specific theta
  return(sum(L_kj))
  
}


negLDELgr_miss3 <- function(par, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d) {
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
  
  
  if(sum(delta) == n_total) {
    stop("No missing data. Use 'drmdel'!")
  }
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
  
  # omega(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  omega_kj <- lapply(1:(m+1), function(r) {
    exp(Qx %*% theta_list[[r]])
  })
  omega_kj_imp <- lapply(1:(m+1), function(r) {
    exp(Qx_imp %*% theta_list[[r]])
  })
  
  omega_kj_tilde <- lapply(1:(m+1), function(r) {
    .vec <- 1/px_hat*omega_kj[[r]] + (1-1/px_hat)*omega_kj_imp[[r]]
    .vec[delta==0,] <- omega_kj_imp[[r]][delta==0,]
    return(.vec)
  })
  
  # h(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  h_r_kj_tilde <- lapply(1:(m+1), function(r) {
    rho[r] * omega_kj_tilde[[r]]
  })
  
  # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  s_kj_tilde <- Reduce("+",h_r_kj_tilde)
  
  
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
  parLk_alpha_tilde <- mapply(function(mat1,mat2,e_k){as.vector(mat1/mat2)-e_k}, mat1 = h_r_kj_tilde, e_k = e_k, MoreArgs = list(mat2 = s_kj_tilde), SIMPLIFY = FALSE)
  
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



#### Idea 4: DEL-EL with missing value ####

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
  
  
  # if(sum(delta) == n_total) {
  #   stop("No missing data. Use 'drmdel'!")
  # }
  # if (length(x) != n_total) {
  #   stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (length(x_imp) != n_total) {
  #   stop("Length of imputed data vector 'x_imp' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (length(delta) != n_total) {
  #   stop("Length of response indicator data vector 'delta' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (any(is.na(x[delta==1]))) {
  #   stop("There are 'NA's in the response!")
  # }  # check data lengths
  # 
  # if (length(px_hat) != n_total) {
  #   stop("Length of estimated response probability data vector 'px_hat' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (sum(n_samples) != n_total) {
  #   stop("total sample size is not consistent with 'n_samples'!")
  # }  # check n_samples lengths
  # 
  # if (length(par) != m*(d+1)) {
  #   stop("parameters length is not m*(d+1)!")
  # }  # check parameter length
  
  
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
  
  
  # if(sum(delta) == n_total) {
  #   stop("No missing data. Use 'drmdel'!")
  # }
  # if (length(x) != n_total) {
  #   stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (length(x_imp) != n_total) {
  #   stop("Length of imputed data vector 'x_imp' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (length(delta) != n_total) {
  #   stop("Length of response indicator data vector 'delta' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (any(is.na(x[delta==1]))) {
  #   stop("There are 'NA's in the response!")
  # }  # check data lengths
  # 
  # if (length(px_hat) != n_total) {
  #   stop("Length of estimated response probability data vector 'px_hat' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (sum(n_samples) != n_total) {
  #   stop("total sample size is not consistent with 'n_samples'!")
  # }  # check n_samples lengths
  # 
  # if (length(par) != m*(d+1)) {
  #   stop("parameters length is not m*(d+1)!")
  # }  # check parameter length
  
  
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



#### Estimate dF_r(y_kj) for each r ####

pEST_func <- function(mele, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d) {
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
  
  
  # if(sum(delta) == n_total) {
  #   stop("No missing data. Use 'drmdel'!")
  # }
  # if (length(x) != n_total) {
  #   stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (length(x_imp) != n_total) {
  #   stop("Length of imputed data vector 'x_imp' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (length(delta) != n_total) {
  #   stop("Length of response indicator data vector 'delta' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (any(is.na(x[delta==1]))) {
  #   stop("There are 'NA's in the response!")
  # }  # check data lengths
  # 
  # if (length(px_hat) != n_total) {
  #   stop("Length of estimated response probability data vector 'px_hat' is not consistent with sample sizes 'n_samples'!")
  # }  # check data lengths
  # 
  # if (sum(n_samples) != n_total) {
  #   stop("total sample size is not consistent with 'n_samples'!")
  # }  # check n_samples lengths
  # 
  # if (length(mele) != m*(d+1)) {
  #   stop("parameters length is not m*(d+1)!")
  # }  # check parameter length
  # 
  
  # rho <- n_samples / n_total    # Calculate rho values for each population
  # theta <- c(rep(0,d+1),mele)    # add theta_0 to parameters
  # x_tilde <- x
  # x_tilde[delta==0] <- x_imp[delta==0]
  # Qx_tilde <- t(sapply(x_tilde, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  # 
  # 
  # if (dim(Qx)[2] != m) { 
  #   stop("dimension of Q(x_kj) should be n*m!")
  # }   # check Q(x_kj)
  # 
  # # Split x and x_imp into a list of samples based on n_samples
  # x_list <- split(x = x_tilde, f = rep(1:(m+1), n_samples))
  # 
  # # Split parameter into a list of parameters based on n_samples
  # theta_list <- split(x = theta, f = rep(1:(m+1), rep(d+1,m+1)))
  # 
  # # omega_r(theta,x_kj): A list of length m+1, each element is a length n_total vector with each row to be rho_r*exp(theta_r^T*Q(x_kj))
  # omega_r_kj <- lapply(1:(m+1), function(r) {
  #   exp(Qx_tilde %*% theta_list[[r]])
  # })
  # 
  # h_r_kj <- mapply("*",omega_r_kj, rho, SIMPLIFY = FALSE)
  # 
  # # s(theta,x_kj): A vector of length n_total with each row to be s(theta,x_kj)
  # s_kj <- Reduce("+",h_r_kj)*n_total
  # 
  # 
  # # eta_kj: A list of length m+1, each element is the partial derivative of L_k(theta,x_kj) w.r.t. alpha
  # eta_kj <- mapply(function(mat1,mat2){as.vector(mat1/mat2)}, mat1 = omega_r_kj, MoreArgs = list(mat2 = s_kj), SIMPLIFY = FALSE)
  # 
  # # partial tilde L_kj/partial theta
  # eta_tilde_kj <- lapply(1:(m+1), function(r) {
  #   #1
  #   .vec <- eta_kj[[r]]
  #   .vec[delta==0] <- eta_kj_imp[[r]][delta==0]
  #   # #2
  #   # .vec <- 1/px_hat * eta_kj[[r]] + (1-1/px_hat) * eta_kj_imp[[r]]
  #   # .vec[delta==0] <- eta_kj_imp[[r]][delta==0]
  #   return(.vec)
  # })
  # 
  # pop_ind <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))
  # 
  # p_est <- data.frame(k=pop_ind, x=x_tilde, p_est=unlist(eta_tilde_kj))
  # return(p_est)
  
###########################
  
  rho <- n_samples / n_total    # Calculate rho values for each population
  theta <- c(rep(0,d+1),mele)    # add theta_0 to parameters
  Qx <- t(sapply(x, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)
  Qx_imp <- t(sapply(x_imp, basis_func)) # Use t() to transpose, so each row represents Q(x_kj)


  # if (dim(Qx)[2] != m) {
  #   stop("dimension of Q(x_kj) should be n*m!")
  # }   # check Q(x_kj)

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
    # #1
    # .vec <- eta_kj[[r]]
    # .vec[delta==0] <- eta_kj_imp[[r]][delta==0]
    #2
    .vec <- 1/px_hat * eta_kj[[r]] + (1-1/px_hat) * eta_kj_imp[[r]]
    .vec[delta==0] <- eta_kj_imp[[r]][delta==0]
    .vec <- .vec/sum(.vec)
    return(.vec)
  })

  pop_ind <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))
  x_tilde <- x
  x_tilde[delta==0] <- x_imp[delta==0]
  p_est <- data.frame(k=pop_ind, x=x_tilde, p_est=unlist(eta_tilde_kj))
  return(p_est)
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


  if(sum(delta) == n_total) {
    stop("No missing data. Use 'drmdel'!")
  }
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
    # #1
    # .vec <- eta_kj[[r]]
    # .vec[delta==0] <- eta_kj_imp[[r]][delta==0]
    #2
    .vec <- 1/px_hat * eta_kj[[r]] + (1-1/px_hat) * eta_kj_imp[[r]]
    .vec[delta==0] <- eta_kj_imp[[r]][delta==0]
    return(.vec)
  })


  # summation over all k,j
  return(c(do.call(rbind, lapply(eta_tilde_kj, sum))))

}




#### Estimate Population CDF ####


cdfEST_func <- function(p_est, n_samples) {
  # Calculate F_r(y) estimator based on DRM.
  # order of x
  x_tilde <- p_est[p_est$k==0,2]
  x_order <- order(x_tilde)
  x_sort <- x_tilde[x_order]
  
  # extract p_est for k^{th} sample
  m <- length(n_samples) - 1
  n_total <- sum(n_samples)
  cum_p <- NULL
  for (i in 1:(m+1)) {
    p_temp <- p_est[p_est$k==(i-1), 3]
    cum_p <- c(cum_p, cumsum(p_temp[x_order]))
  }
  
  pop_ind <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))
  x_data <- rep(x=x_sort, m+1)
  cdf_est <- data.frame(k=pop_ind, x=x_data, cdf_est=cum_p)
  
  return(cdf_est)
  
}



#### Estimate Population Quantile ####


quantileEST_func <- function(k, p, cdf_est) {
  # Calculate estimated quantile based on DRM.
  
  # check length
  if ((length(k) != length(p)) && (length(k) != 1 ) && (length(p) != 1)) {
    stop("The lengths of vector 'k' and 'p' must be the same, or one of the vector must have length one!")
  }
  nK <- max(length(k), length(p))
  if (length(k)==1) k <- rep(k, nK)
  if (length(p)==1) p <- rep(p, nK)
  
  x_sort <- cdf_est[cdf_est$k==0,2]
  
  q_est <- numeric(nK)
  
  for (i in 1:nK) {
    
    cdf_est_tmp <- cdf_est[cdf_est$k==k[i], 3]
    
    if (p[i] >= tail(cdf_est_tmp, 1)) {
      pos_tmp <- length(x_sort) - 1
    } else {
      pos_tmp <- length(which(cdf_est_tmp < p[i]))
    }
    
    # # linear interpolation
    # q_est[i] <- x_sort[pos_tmp] + (x_sort[pos_tmp + 1] - x_sort[pos_tmp]) *
    #   (p[i] - cdf_est_tmp[pos_tmp]) /
    #   (cdf_est_tmp[pos_tmp+1] - cdf_est_tmp[pos_tmp])
    
    # no interpolation
    q_est[i] <- x_sort[pos_tmp]
    
  }
  
  return(q_est)
  
}


#### Fit EL-DRM ####

eldrm_miss <- function(x, x_imp, delta, px_hat, n_samples, basis_func) {
  # This function fit a DRM with EL inference
  # x: A vector of all samples from all populations (numeric vector)
  # x_imp: A vector of all imputed values from all populations (numeric vector)
  # delta: A vector of all response indicator from all populations (0,1)
  # px_hat: A vector of all estimated propensity score from all populations (numeric vector)
  # basis_func: A function to compute Q(x) given an observation x
  # n_samples: A vector of sample sizes for each population
  
  # n_total: Total sample size (sum of all n_samples)
  # m+1: number samples
  # d+1: dimension of basis function
  
  m <- length(n_samples) - 1
  n_total <- sum(n_samples)
  d <- length(basis_func(1)) - 1
  
  
  if(sum(delta) == n_total) {
    stop("No missing data. Use 'drmdel'!")
  }
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(x_imp) != n_total) {
    stop("Length of imputed data vector 'x_imp' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (length(delta) != n_total) {
    stop("Length of response indicator data vector 'delta' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (anyNA(x[delta==1])) {
    stop("There are 'NA's in the response!")
  }  # check data lengths
  
  if (!anyNA(x)) {
    stop("Missing value need to be 'NA's in the response!")
  }  # check data lengths
  
  if (length(px_hat) != n_total) {
    stop("Length of estimated response probability data vector 'px_hat' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  
  if (sum(n_samples) != n_total) {
    stop("total sample size is not consistent with 'n_samples'!")
  }  # check n_samples lengths
  
  
  # calculating MELE
  par_dim <- m*(d+1)
  par_init <- rep(0, par_dim)
  
  # new del
  drm_opt <- do.call(optim, c(list(par=par_init, fn=negLDELfn_miss4,
                                   gr=negLDELgr_miss4,
                                   method="BFGS", control=list(maxit=10000),
                                   x=x, x_imp=x_imp, delta=delta, px_hat=px_hat,
                                   n_total=n_total, n_samples=n_samples, m=m,
                                   basis_func=basis_func,
                                   d=d)))
  par_mele <- drm_opt$par
  
  
  # Obtain the estimated p_kj 
  p_est <- pEST_func(par_mele, x, x_imp, delta, px_hat, n_total, n_samples, m, basis_func, d)
  
  # Obtain the estimated cdf F(x_kj) 
  cdf_est <- cdfEST_func(p_est, n_samples)
  
  return(list(mele=par_mele,
              negldl=drm_opt$value,
              p_est=p_est, cdf_est=cdf_est))
  
}


