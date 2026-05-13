library(fda)
# generation of coefficients for the three scenarios in the paper

# scenario 1 ####
# data generating function for scenario 1
# N: number of generated curves (defalut 100)
# replications: number of datasets generated (default 100)
# current_seed: sets the seed (default 1)

# outputs a list of two elements:
# 1: "coefficients" is a list of length replications, each element of which is a matrix of coefficients
# 2: "true_labels" is a list of length replications, each element of which is a vector of labels assignments

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

# scenario 2 ####
# data generating function for scenario 2
# N: number of generated curves (defalut 150)
# replications: number of datasets generated (default 100)
# current_seed: sets the seed (default 1)

# outputs a list of two elements:
# 1: "coefficients" is a list of length replications, each element of which is a matrix of coefficients
# 2: "true_labels" is a list of length replications, each element of which is a vector of labels assignments

scenario2 <- function(N = 150, replications = 100, current_seed = 1) {
  # number of clusters
  true_G <- 3
  # group assignment probability
  p1 <- 1/3
  # noise
  sigma2 <- 1
  # functions
  h1 <- function(x){matrix(pmax(6-abs(x-11),0))}
  h2 <- function(x){matrix(pmax(6-abs(x-15),0))}
  h3 <- function(x){matrix(pmax(6-abs(x-7),0))}
  
  u <- function(n) {matrix(runif(n))}
  et <- function(t, n) {
    replicate(n, rnorm(length(t), 0, sqrt(sigma2)))
  }
  
  # timepoints
  n_points <- 1001
  t <- seq(1,21,length.out=n_points)
  T_0 <- t[1]
  T_n <- t[length(t)]
  
  # parameters found by GCV
  best_basi <- 200
  best_lambda <- 10
  
  basi <- create.bspline.basis(rangeval = c(T_0,T_n), nbasis = best_basi)
  fd_Par <- fdPar(basi, lambda = best_lambda)
  
  set.seed(current_seed)
  
  clusters <- t(sapply(1:replications, function(i) {
    table(sample(1:G_vero, N, replace = TRUE, prob = rep(p1, G_vero))) # prob can change
  }))
  
  data <- lapply(1:replications, function(idx) {
    n1 <- clusters[idx,]
    
    u1 <- u(n1[1])
    u2 <- u(n1[2])
    u3 <- u(n1[3])
    X1 <- h1(t)%*%t(u1) + h2(t)%*%t(1-u1) + et(t,n1[1])
    X2 <- h1(t)%*%t(u2) + h3(t)%*%t(1-u2) + et(t,n1[2])
    X3 <- h2(t)%*%t(u3) + h3(t)%*%t(1-u3) + et(t,n1[3])
    
    cbind(X1,X2,X3)
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


# scenario 3 ####
# data generating function for scenario 3
# N: number of generated curves (defalut 100)
# replications: number of datasets generated (default 100)
# current_seed: sets the seed (default 1)

# outputs a list of two elements:
# 1: "coefficients" is a list of length replications, each element of which is a matrix of coefficients
# 2: "true_labels" is a list of length replications, each element of which is a vector of labels assignments

scenario3 <- function(N = 100, replications = 100, current_seed = 1) {
  # number of clusters
  true_G <- 4
  # group assignment probability
  p1 <- 1/4
  # noise
  sigma2 <- 0.2^2
  # auxiliary functions
  psi1 <- function(x, wiggle=F) {
    if(!wiggle) {
      y <- as.vector(x+sin(pi*x)*exp(-x))
    } 
    else {
      y <- as.vector(cos(20*x))
    }
    y <- y-mean(y)
    s1 <- sum(y*y)
    y <- y*s1^-0.5
    return(matrix(y))
  }
  psi2 <- function(x, wiggle=F) {
    if(!wiggle) {
      y <- as.vector(cos(3+pi*x))
    } 
    else {
      y <- as.vector(sin(20*x))
    }
    y <- y-mean(y)
    s1 <- sum(y*y)
    y <- y*s1^-0.5
    return(matrix(y))
  }
  # functions
  m1 <- function(x,wiggle=F){psi1(x,wiggle)}
  m2 <- function(x,wiggle=F){psi2(x,wiggle)}
  m3 <- function(x,wiggle1=F,wiggle2=F){psi1(x,wiggle1)+psi2(x,wiggle2)}
  m4 <- function(x,wiggle1=F,wiggle2=F){-psi1(x,wiggle1)+psi2(x,wiggle2)}
  
  et <- function(t, n) {
    replicate(n, rnorm(length(t), 0, sqrt(sigma2)))
  }
  
  
  # timepoints
  n_points <- 101
  t <- seq(-1,1,length.out=n_points)
  T_0 <- t[1]
  T_n <- t[length(t)]
  
  # parameters
  best_basi <- 101
  best_lambda <- 1e-4
  
  basi <- create.bspline.basis(rangeval = c(T_0,T_n), nbasis = best_basi)
  fd_Par <- fdPar(basi, lambda = best_lambda)
  
  set.seed(current_seed)
  
  clusters <- t(sapply(1:replications, function(i) {
    table(sample(1:G_vero, N, replace = TRUE, prob = rep(p1, G_vero))) # prob can change
  }))
  
  data <- lapply(1:replications, function(idx) {
    n1 <- clusters[idx,]
    
    X1 <- m1(t,wiggle = T)%*%t(rep(1, n1[1])) + et(t,n1[1])
    X2 <- m2(t,wiggle = T)%*%t(rep(1, n1[2])) + et(t,n1[2])
    X3 <- m3(t,wiggle1 = T, wiggle2 = T)%*%t(rep(1, n1[3])) + et(t,n1[3])
    X4 <- m4(t,wiggle = T, wiggle2 = T)%*%t(rep(1, n1[4])) + et(t,n1[4])
    
    cbind(X1,X2,X3,X4)
  })
  
  coefficients <- lapply(1:ripetizioni, function(idx) {
    fd_data <- smooth.basis(y = as.matrix(data[[idx]]), fdParobj = fd_Par, argvals = t)
    
    t(fd_data$fd$coefs)
  })
  
  true_labels <- lapply(1:replications, function(idx) {
    rep(1:true_G, clusters[idx,])
  })
  
  return(list("coefficients"=coefficients, "true_labels"=true_labels))
}
