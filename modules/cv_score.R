# code used to make the score function of interest
library(ape)
library(dendextend)
library(cluster)
library(tibble)
library(magrittr)
library(dplyr)
library(phytools)
library(limSolve)
library(mltools)
library(data.table)
library(factoextra)
source("convert_to_parenthesis.R")
library(foreach)
library(doParallel)
library(tictoc)

# defining one hot encoder for simap
onehotencoder = function(miss_data){
  nlvls = nlevels(miss_data)
  onehot = matrix(0, nrow = length(miss_data), ncol = nlvls)
  factors = as.numeric(miss_data)
  for (i in 1:length(factors)){
    if(is.na(factors[i]) == TRUE){
      onehot[i, ] = rep(1/nlvls, nlvls)
    }else{
      onehot[i, factors[i]] = 1
    }
  }
  row.names(onehot) = names(miss_data)
  colnames(onehot) = levels(miss_data)
  return(onehot)
}


continuous.missing = function(dend, miss_data, tol){
  max.tl = miss_data
  del.ind = match(NA, max.tl)
  new_max.tl = max.tl[-del.ind]
  fit = anc.ML(dend, new_max.tl, model = "BM", tol = tol)
  y.hat = fit$missing.x
  y.hat.var = fit$sig2
  return(list(y.hat, y.hat.var))
}

factorial.missing = function(dend, miss_data){
  onehot = onehotencoder(miss_data)
  fit_discrete = make.simmap(dend, onehot, model = "ER", pi = "estimated",
                             message = FALSE, nsim = 10)
  Q = as.matrix(fit_discrete$Q)
  dim = nrow(Q)
  pi = xranges(E = matrix(1, nrow = dim, ncol = dim), F = rep(1, dim),
          G = t(Q), H = rep(0, dim))[, 1]
  A.norm = dim
  return(list(pi, A.norm))
}

row_computing = function(types, dend, original_data, tol, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = row.names(original_data)
  j = col
  if (types[j] == "integer" | types[j] == "numeric"){
    current_data = original_data[, j]
    scaled = scale(current_data)
    for (i in (1:n)){
      saved_value = scaled[i]
      miss_data = scaled
      names(miss_data) = names
      miss_data[i] = NA
      results = continuous.missing(dend, miss_data, tol)
      y.pred = results[[1]]
      y.pred_var = results[[2]]
      mse = ((y.pred - saved_value)^2)/y.pred_var
      lines[i] =  mse
    }
  }else{
    current_data = original_data[, j]
    if (types[j] == "logical") current_data = as.factor(current_data)
    for (i in (1:n)){
      saved_value = numeric(nlevels(current_data))
      saved_value[as.numeric(current_data)[i]] = 1
      miss_data = current_data
      names(miss_data) = names
      miss_data[i] = NA
      results = factorial.missing(dend, miss_data)
      pi = results[[1]]
      A.norm = results[[2]]
      mse = (sum((saved_value - pi)^2))/A.norm
      lines[i] = mse
    }
  }
  return(lines)
}

# paralelized
L_score = function(dend, original_data, tol = 1e-18){
  types = sapply(original_data, class)
  p = length(original_data[1, ])
  n = length(original_data[, 1])
  total = p*n
  score.matrix = matrix(0, nrow = n, ncol = p)
  names = row.names(original_data)
  row.names(original_data) = c(1:nrow(original_data))
  # paralellizing
  cores = detectCores()
  cl = makeCluster(cores[1] - 1)
  clusterExport(cl, c("row_computing", "factorial.missing", "continuous.missing",
                      "onehotencoder"))
  registerDoParallel(cl)
  score.matrix = foreach(j = 1:p, .combine = cbind,
                         .export = c("row_computing", "factorial.missing", 
                          "continuous.missing","onehotencoder"),
                         .packages = c("ape", "phytools")) %dopar% {
    lines = row_computing(types, dend, original_data, tol, j)
    lines
}
  stopCluster(cl)
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}

# serialized option
L_score_2 = function(dend, original_data, tol = 1e-18){
  types = sapply(original_data, class)
  p = length(original_data[1, ])
  n = length(original_data[, 1])
  total = p*n
  score.matrix = matrix(0, nrow = n, ncol = p)
  names = row.names(original_data)
  row.names(original_data) = c(1:nrow(original_data))
  for (j in (1:p)){
    lines = numeric(n)
    if (types[j] == "integer" | types[j] == "numeric"){
      current_data = original_data[, j]
      scaled = scale(current_data)
      for (i in (1:n)){
        saved_value = scaled[i]
        miss_data = scaled
        names(miss_data) = names
        miss_data[i] = NA
        results = continuous.missing(dend, miss_data, tol)
        y.pred = results[[1]]
        y.pred_var = results[[2]]
        mse = ((y.pred - saved_value)^2)/y.pred_var
        lines[i] =  mse
      }
    }else{
      current_data = original_data[, j]
      if (types[j] == "logical") current_data = as.factor(current_data)
      for (i in (1:n)){
        saved_value = numeric(nlevels(current_data))
        saved_value[as.numeric(current_data)[i]] = 1
        miss_data = current_data
        names(miss_data) = names
        miss_data[i] = NA
        results = factorial.missing(dend, miss_data)
        pi = results[[1]]
        A.norm = results[[2]]
        mse = (sum((saved_value - pi)^2))/A.norm
        lines[i] = mse
      }
    }
    score.matrix[, j] = lines  
  }
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}

