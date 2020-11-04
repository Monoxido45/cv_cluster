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
source("C:/Users/lucru/Estatística_UFSCar/cv_cluster/modules/convert_to_parenthesis.R")
library(foreach)
library(doParallel)
library(tictoc)
library(mvMORPH)

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
  cname = colnames(original_data)[j]
  if (types[j] == "integer" | types[j] == "numeric"){
    current_data = as.matrix(original_data[, j])
    # not scaling
    colnames(current_data) = cname
    rownames(current_data) = names
    fit = mvBM(dend, current_data, model = "BMM", echo = T, method = "inverse")
    for (i in (1:n)){
      saved_value = current_data[i]
      miss_data = current_data
      miss_data[i, 1] = NA
      imp = estim(dend, miss_data, fit)
      y.pred = imp$estimates[imp$NA_index, 1]
      y.pred_var = imp$var[imp$NA_index]
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

row_computing2 = function(types, dend, original_data, tol, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = row.names(original_data)
  j = col
  cname = colnames(original_data)[j]
  if (types[j] == "integer" | types[j] == "numeric"){
    bool = (types == "numeric" | types == "integer")
    selected_data = as.data.frame(original_data[, bool])
    row.names(selected_data) = names
    colnames(selected_data) = colnames(original_data)[bool]
    fit = mvBM(dend, selected_data, model = "BMM", echo = T, method = "inverse")
    current_data = as.matrix(original_data[, j])
    # not scaling
    for (i in (1:n)){
      saved_value = current_data[i]
      miss_data = selected_data
      miss_data[i, j] = NA
      imp = estim(dend, miss_data, fit)
      y.pred = imp$estimates[imp$NA_index, 1]
      y.pred_var = imp$var[imp$NA_index]
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

row_computing3 = function(types, dend, original_data, tol, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = row.names(original_data)
  j = col
  cname = colnames(original_data)[j]
  if (types[j] == "integer" | types[j] == "numeric"){
    current_data = as.matrix(original_data[, j])
    # not scaling
    colnames(current_data) = cname
    rownames(current_data) = names
    for (i in (1:n)){
      saved_value = current_data[i]
      miss_data = current_data
      miss_data = as.data.frame(miss_data)
      colnames(miss_data) = colnames(current_data)
      rownames(miss_data) = rownames(current_data)
      miss_data[i, 1] = NA
      fit = mvBM(sim.test, miss_data, model = "BMM", method = "inverse")
      imp = estim(dend, miss_data, fit)
      y.pred = imp$estimates[imp$NA_index, 1]
      y.pred_var = imp$var[imp$NA_index]
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

row_computing4 = function(types, dend, original_data, tol, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = row.names(original_data)
  j = col
  cname = colnames(original_data)[j]
  if (types[j] == "integer" | types[j] == "numeric"){
    bool = (types == "numeric" | types == "integer")
    selected_data = as.data.frame(original_data[, bool])
    row.names(selected_data) = names
    colnames(selected_data) = colnames(original_data)[bool]
    current_data = original_data[, j]
    # not scaling
    for (i in (1:n)){
      saved_value = current_data[i]
      miss_data = selected_data
      miss_data[i, j] = NA
      fit = mvBM(dend, miss_data, model = "BMM", echo = T, method = "inverse")
      imp = estim(dend, miss_data, fit)
      y.pred = imp$estimates[imp$NA_index, 1]
      y.pred_var = imp$var[imp$NA_index]
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


# using mvMORPH
L_score = function(dend, original_data, tol  = 1e-20){
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
                         .export = c("row_computing", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                           lines = row_computing2(types, dend, original_data, tol, j)
                           lines
                         }
  stopCluster(cl)
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}


L_score_2 = function(dend, original_data, tol  = 1e-20){
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
                         .export = c("row_computing2", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                           lines = row_computing2(types, dend, original_data, tol, j)
                           lines
                         }
  stopCluster(cl)
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}

L_score_3 = function(dend, original_data, tol  = 1e-20){
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
                         .export = c("row_computing3", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                           lines = row_computing2(types, dend, original_data, tol, j)
                           lines
                         }
  stopCluster(cl)
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}

L_score_4 = function(dend, original_data, tol  = 1e-20){
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
                         .export = c("row_computing4", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                           lines = row_computing2(types, dend, original_data, tol, j)
                           lines
                         }
  stopCluster(cl)
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}



