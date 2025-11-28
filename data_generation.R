# scenario 1 ####
# data generating function for scenario 1
# N: number of generated curves (defalut 100)
# replications: number of datasets generated (default 100)
# current_seed: sets the seed (default 1)

# outputs a list of two elements:
# 1: "coefficients" is a list of length replications, each element of which is a matrix of coefficients
# 2: "true_labels" is a list of length replications, each element of which is a vector of labels assignments

library(fda)
scenario1 <- function(N = 100, replications = 100, current_seed = 1) {
  # number of clusters
  true_G <- 2
  # group 1 assignment probability
  p1 <- 0.5
  # noise
  sigma2 <- 1/12
  # functions
  h1 <- function(x){matrix(pmax(6-abs(x-7),0))}
  h2 <- function(x){matrix(pmax(6-abs(x-15),0))}
  
  u12 <- function(n) {matrix(rnorm(n, 0, sqrt(sigma2)))}
  et <- function(t, n) {
    replicate(n, rnorm(length(t), 0, sqrt(sigma2)))
  }
  
  # timepoints
  n_points <- 1001
  t <- seq(1,21,length.out=n_points)
  T_0 <- t[1]
  T_n <- t[length(t)]
  
  # parameters found by GCV
  best_basi <- 100
  best_lambda <- 1
  
  basi <- create.bspline.basis(rangeval = c(T_0,T_n), nbasis = best_basi)
  fd_Par <- fdPar(basi, lambda = best_lambda)
  
  set.seed(current_seed)
  
  n1_vec <- rbinom(replications, N, 1 - p1)
  n2_vec <- N - n1_vec
  clusters <- cbind(n1 = n1_vec, n2 = n2_vec)
  
  data <- lapply(1:replications, function(idx) {
    n1 <- clusters[idx,1]
    n2 <- clusters[idx,2]
    
    X1 <- h1(t) %*% t(u12(n1)) + h2(t) %*% t(u12(n1)) + et(t, n1)
    X2 <- h1(t) %*% t(u12(n2)) + et(t, n2)
    
    cbind(X1, X2)
  })
  
  coefficients <- lapply(1:replications, function(idx) {
    fd_data <- smooth.basis(y = as.matrix(data[[idx]]), fdParobj = fd_Par, argvals = t)
    
    t(fd_data$fd$coefs)
  })
  
  true_labels <- lapply(1:replications, function(idx) {
    rep(1:true_G, clusters[idx,])
  })
  
  return(list("coefficients"=coefficients, "true_labels"=true_labels))
}
