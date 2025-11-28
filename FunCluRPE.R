# Main function - for cycle ####
# coefficients: matrix of coefficients.
# alternatively: fd_object: object of class fd, from fda package
# rp_type: string, either "haar" (default) or "gaussian" for the type of RPs to use
# criterion: which sorting criterion to use, one of "kl" (Kullback-Leibler, default), "wass" (Wasserstein), "ent" (entropy)
# d: integer, dimension of the projected space
# alternatively: G, integer. Hypothesised number of groups, used in the formula ceiling(5*log(G))+1 to compute d
# B (default 1000) number of RPs
# Bstar (default 100) number of retained RPs
# try_models: one or a combination of the models for Mclust (see mclust.options("emModelNames")). Default is "all", meaning all models are tried out.
# G_gmm: integer vector specifying the numbers of mixture components used in all GMMs, if not specified 2:9
# crisp: TRUE for hard clustering, FALSE for fuzzy clustering
# current_seed: sets a seed. defalut is 1

# Returned output: a list
## final_cluster_output: either a vector of cluster labels or a matrix of posterior probabilities;
## proj_cluster: N times B matrix with the single projection partial clustering results
## criterion: which criterion was used in ordering the partial clustering results
## criterion_values: vector of length B with the criterion values for the corresponding projection

# libraries needed
library(mclust)
library(expm)
library(clue)

FunCluRPE <- function(coefficients, fd_object, 
                      rp_type = "haar", criterion = "kl", 
                      d, G, B = 1000, Bstar = 100, try_models = "all", G_gmm, crisp = TRUE, current_seed = 1) {
  if (!missing(coefficients) + !missing(fd_object) != 1) {
    stop("Supply exactly one of 'coefficients' or 'fd_object'.")
  }
  
  if (!missing(coefficients)) {
    if (!is.matrix(coefficients) |
        !is.numeric(coefficients) |
        any(!is.finite(coefficients))) {
      stop("'coefficients' must be a numeric matrix of finite real values.")
    }
  }
  
  if (!missing(fd_object)) {
    if (!inherits(fd_object, "fd")) {
      stop("'fd_object' must be an object of class 'fd'.")
    }
  }
  
  if (!is.character(rp_type) || length(rp_type) != 1) {
    stop("'rp_type' must be a single character value: 'haar' or 'gaussian'.")
  }
  
  if (!(rp_type %in% c("haar", "gaussian"))) {
    stop("'rp_type' must be either 'haar' or 'gaussian'.")
  }
  
  if (!is.character(criterion) ||
      length(criterion) != 1 ||
      !(criterion %in% c("kl", "wass", "ent"))) {
    stop("'criterion' must be one of: 'kl', 'wass', 'ent'.")
  }
  
  if (missing(d)) {
    if (!missing(G)) {
      d <- ceiling(5*log(G))+1
    } else {
      stop("Specify either d or G")
    }
  }
  
  if (!is.numeric(B) || length(B) != 1 ||
      !is.finite(B) || B %% 1 != 0) {
    stop("'B' must be a single finite integer.")
  }
  
  if (!is.numeric(Bstar) || length(Bstar) != 1 ||
      !is.finite(Bstar) || Bstar %% 1 != 0) {
    stop("'Bstar' must be a single finite integer.")
  }
  
  if (!(is.character(try_models) && length(try_models) >= 1)) {
    stop("'try_models' in invalid.")
  }
  
  if (length(try_models) == 1 && try_models == "all") {
    # valid
  } else {
    # Must be a subset of valid model names
    if (!all(try_models %in% mclust.options("emModelNames"))) {
      stop("'try_models' in invalid.")
    }
  }
  
  if (!missing(G_gmm)) {
    if (!is.numeric(G_gmm) || length(G_gmm) != 1 ||
        !is.finite(G_gmm) || G_gmm %% 1 != 0 ||
        G_gmm < 2) {
      stop("'G_gmm' must be an integer greater than or equal to 2.")
    }
  } else {
    G_gmm <- 2:9
  }
  
  if (!is.logical(crisp) || length(crisp) != 1) {
    stop("'crisp' must be a logical value (TRUE or FALSE).")
  }
  
  if (!is.numeric(current_seed) || length(current_seed) != 1 || !is.finite(current_seed)) {
    stop("'current_seed' must be a single finite numeric value.")
  }
  
  # random matrices generation
  dim_coef <- dim(coefficients)[2]
  set.seed(current_seed)
  if (rp_type == "gaussian") {
    Pmat_list <- replicate(B,
                           matrix(1 / sqrt(dim_coef) * rnorm(dim_coef * d, 0, 1), nrow = dim_coef),
                           simplify = FALSE)
    
  }
  if (rp_type == "haar") {
    Pmat_list <- replicate(B,
                           qr.Q(qr(matrix(
                             1 / sqrt(dim_coef) * rnorm(dim_coef * d, 0, 1), nrow = dim_coef
                           )[, 1:d]))[, 1:d],
                           simplify = FALSE)
  }
  
  if (try_models == "all") {
    try_models <- mclust.options("emModelNames")
  }
  
  # for loop
  
  proj_cluster <- matrix(nrow = nrow(coefficients), ncol = 0)
  proj_G <- vector()
  criterion_values <- vector()
  
  for (B_i in 1:B) {
    coefficients_pro <- coefficients%*%Pmat_list[[B_i]]
    
    # mclust
    set.seed(current_seed)
    normal_clustering <- Mclust(coefficients_pro, G, modelNames = try_models, verbose = FALSE)
    
    mclust_cl <- normal_clustering$classification
    mclust_posteriors <- normal_clustering$z
    mclust_G <- normal_clustering$G
    
    # split according to groups
    coefficients_pro_split <- split.data.frame(coefficients_pro, mclust_cl)
    
    # extract gmm parameters and make them lists
    mu_group <- as.list(data.frame(normal_clustering$parameters$mean))
    sigma_group <- normal_clustering$parameters$variance$sigma
    sigma_group <- lapply(seq(dim(sigma_group)[3]), function(y) sigma_group[,,y])
    pro_group <- normal_clustering$parameters$pro
    
    temp <- expand.grid(x = 1:mclust_G, y = 1:mclust_G)
    temp <- subset(temp, x != y)
    temp <- list(temp$x, temp$y)
    
    if (criterion == "kl") {
      kl_base <- mapply(
        function(mean_index, distr_index)
        {
          (log(determinant(sigma_group[[distr_index]], logarithm = F)$modulus/determinant(sigma_group[[mean_index]], logarithm = F)$modulus)
           - d
           + sum(diag(solve(sigma_group[[distr_index]], tol=1e-20)%*%sigma_group[[mean_index]])) 
           + t(mu_group[[distr_index]] - mu_group[[mean_index]])%*%solve(sigma_group[[distr_index]], tol=1e-20)%*%(mu_group[[distr_index]]-mu_group[[mean_index]]))/2
        },
        mean_index = temp[[1]],
        distr_index = temp[[2]]
      )
      
      criterion_values <- c(criterion_values, mean(kl_base, na.rm = T))
    }
    
    if (criterion == "wass") {
      temp <- asplit(combn(1:mclust_G, 2),1)
      Wasserstein <- mapply(
        function(v1_index, v2_index)
        {
          (sum((mu_group[[v1_index]] - mu_group[[v2_index]])^2) 
           + sum(diag(sigma_group[[v1_index]] + sigma_group[[v2_index]] 
                      - 2*sqrtm(sqrtm(sigma_group[[v1_index]])%*%sigma_group[[v2_index]]%*%sqrtm(sigma_group[[v1_index]])))))
        },
        v1_index = temp[[1]],
        v2_index = temp[[2]])
      
      criterion_values <- c(criterion_values, mean(Wasserstein, na.rm = T))
    }
    
    if (criterion == "ent") {
      Entropy <- -sum(mclust_posteriors * log(pmax(mclust_posteriors, 1e-12))) / (nrow(mclust_posteriors) * log(ncol(mclust_posteriors)))
      
      criterion_values <- c(criterion_values, mean(Entropy, na.rm = T))
    }
    
    proj_cluster <- cbind(proj_cluster, apply(mclust_posteriors, 1, which.max))
    proj_G <- c(proj_G, ncol(mclust_posteriors))
  }
  
  if (criterion%in%c("kl", "wass")) {
    ordered_cluster <- proj_cluster[,order(criterion_values, decreasing = T)]
  }
  if (criterion%in%c("ent")) {
    ordered_cluster <- proj_cluster[,order(criterion_values, decreasing = F)]
  }
  
  set.seed(current_seed)
  mclust_sorted <- as.data.frame(ordered_cluster[,1:Bstar])
  membership <- lapply(mclust_sorted, function(x) as.cl_membership(x))
  ensemble_fuzzy <- cl_consensus(membership, method = "SE")$.Data
  
  if (crisp) {
    final_cluster_output <- factor(apply(ensemble_fuzzy,1,which.max))
  } else {
    final_cluster_output <- ensemble_fuzzy
  }
  
  results <- list("final_cluster_output" = final_cluster_output,
                  "proj_cluster" = proj_cluster,
                  "criterion" = criterion,
                  "criterion_values" = criterion_values)
  
  return(results)
}
