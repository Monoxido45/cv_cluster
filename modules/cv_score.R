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


factorial.missing = function(dend, miss_data, row){
  onehot = onehotencoder(miss_data)
  fit_discrete = make.simmap(as.multiPhylo(dend), onehot, model = "ER",
                             message = FALSE)
  results = describe.simmap(fit_discrete, plot = F)
  which_row = which(rownames(results$tips) == row)
  pi = results$tips[which_row, ]
  return(pi)
}


row_computing = function(types, dend, original_data, tol,seed, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = row.names(original_data)
  j = col
  tree.order = dend$tip.label
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
      names(miss_data) = names
      miss_data[i, 1] = NA
      imp = estim(dend, miss_data, fit)
      y.pred = imp$estimates[i, 1]
      y.pred_var = imp$var[i, 1]
      mse = ((y.pred - saved_value)^2)/y.pred_var
      lines[i] =  mse
    }
  }else{
    current_data = original_data[, j]
    if (types[j] == "logical") current_data = as.factor(current_data)
    set.seed(seed)
    for (i in (1:n)){
      saved_value = numeric(nlevels(current_data))
      saved_value[as.numeric(current_data)[i]] = 1
      miss_data = current_data
      A = nlevels(miss_data)
      names(miss_data) = names
      miss_data[i] = NA
      pi = factorial.missing(dend, miss_data, i)
      mse = (sum((saved_value - pi)^2))
      lines[i] = mse/A
    }
  }
  return(lines)
}

row_computing2 = function(types, dend, original_data, tol, seed, col){
  n = nrow(original_data)
  lines = numeric(n)
  row.names(original_data) = c(1:nrow(original_data))
  j = col
  cname = colnames(original_data)[j]
  if (types[j] == "integer" | types[j] == "numeric"){
    bool = (types == "numeric" | types == "integer")
    selected_data = as.data.frame(original_data[, bool])
    row.names(selected_data) = row.names(original_data)
    colnames(selected_data) = colnames(original_data)[bool]
    fit = mvBM(dend, selected_data, model = "BMM", echo = T, method = "inverse")
    current_data = as.matrix(original_data[, j])
    # not scaling
    for (i in (1:n)){
      saved_value = current_data[i]
      miss_data = selected_data
      miss_data[i, j] = NA
      imp = estim(dend, miss_data, fit)
      print(imp$estimates)
      print(imp$var)
      print(saved_value)
      y.pred = imp$estimates[i, j]
      y.pred_var = imp$var[i, j]
      mse = ((y.pred - saved_value)^2)/y.pred_var
      lines[i] =  mse
    }
  }else{
    current_data = original_data[, j]
    if (types[j] == "logical") current_data = as.factor(current_data)
    set.seed(seed)
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

row_computing3 = function(types, dend, original_data, tol, seed, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = c(1:nrow(original_data))
  j = col
  cname = colnames(original_data)[j]
  tree.order = dend$tip.label
  if (types[j] == "integer" | types[j] == "numeric"){
    bool = (types == "numeric" | types == "integer")
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
      fit = mvBM(dend, miss_data, model = "BMM", method = "inverse")
      imp = estim(dend, miss_data, fit)
      y.pred = imp$estimates[i, 1]
      y.pred_var = imp$var[i ,1]
      mse = ((y.pred - saved_value)^2)/y.pred_var
      lines[i] =  mse
    }
  }else{
    current_data = original_data[, j]
    if (types[j] == "logical") current_data = as.factor(current_data)
    set.seed(seed)
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

row_computing4 = function(types, dend, original_data, tol, seed, col){
  n = nrow(original_data)
  lines = numeric(n)
  names = c(1:nrow(original_data))
  j = col
  cname = colnames(original_data)[j]
  if (types[j] == "integer" | types[j] == "numeric"){
    bool = (types == "numeric" | types == "integer")
    selected_data = as.data.frame(original_data[, bool])
    colnames(selected_data) = colnames(original_data)[bool]
    row.names(selected_data) = names
    current_data = as.matrix(original_data[, j])
    # not scaling
    for (i in (1:n)){
      saved_value = current_data[i]
      miss_data = selected_data
      miss_data[i, j] = NA
      fit = mvBM(dend, miss_data, model = "BMM", echo = T, method = "inverse")
      imp = estim(dend, miss_data, fit)
      print(imp$estimates)
      print(imp$var)
      print(saved_value)
      y.pred = imp$estimates[i, j]
      y.pred_var = imp$var[i, j]
      mse = ((y.pred - saved_value)^2)/y.pred_var
      lines[i] =  mse
    }
  }else{
    current_data = original_data[, j]
    if (types[j] == "logical") current_data = as.factor(current_data)
    set.seed(seed)
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
L_score = function(dend, original_data, tol  = 1e-20, seed = 99){
  types = sapply(original_data, class)
  p = length(original_data[1, ])
  n = length(original_data[, 1])
  total = p*n
  score.matrix = matrix(0, nrow = n, ncol = p)
  dend$edge.length[which(dend$edge.length %in% c(0))] = 10^(-3)
  names = row.names(original_data)
  row.names(original_data)  = c(1:nrow(original_data))
  if(is.null(row.names(original_data))) row.names(original_data) = c(1:nrow(original_data))
  # paralellizing
  cores = detectCores()
  cl = makeCluster(cores[1] - 1)
  clusterExport(cl, c("row_computing", "factorial.missing", "continuous.missing",
                      "onehotencoder"))
  registerDoParallel(cl)
  score.matrix = foreach(j = 1:p, .combine = cbind,
                         .export = c("row_computing", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                           lines = row_computing(types, dend, original_data, tol, seed, j)
                           lines
                         }
  stopCluster(cl)
  if(p == 1){
    score.matrix = as.matrix(score.matrix)
  }
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}


L_score_2 = function(dend, original_data, tol  = 1e-20, seed = 99){
  types = sapply(original_data, class)
  p = length(original_data[1, ])
  n = length(original_data[, 1])
  total = p*n
  score.matrix = matrix(0, nrow = n, ncol = p)
  names = row.names(original_data)
  row.names(original_data)  = c(1:nrow(original_data))
  if(is.null(row.names(original_data))) row.names(original_data) = c(1:nrow(original_data))
  # paralellizing
  cores = detectCores()
  cl = makeCluster(cores[1] - 1)
  clusterExport(cl, c("row_computing2", "factorial.missing", "continuous.missing",
                      "onehotencoder"))
  registerDoParallel(cl)
  score.matrix = foreach(j = 1:p, .combine = cbind,
                         .export = c("row_computing2", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                           lines = row_computing2(types, dend, original_data, tol, seed, j)
                           lines
                         }
  stopCluster(cl)
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}

L_score_3 = function(dend, original_data, tol  = 1e-20, seed = 99){
  types = sapply(original_data, class)
  p = length(original_data[1, ])
  n = length(original_data[, 1])
  total = p*n
  score.matrix = matrix(0, nrow = n, ncol = p)
  names = row.names(original_data)
  row.names(original_data)  = c(1:nrow(original_data))
  if(is.null(row.names(original_data))) row.names(original_data) = c(1:nrow(original_data))
  # paralellizing
  cores = detectCores()
  cl = makeCluster(cores[1] - 1)
  clusterExport(cl, c("row_computing", "factorial.missing", "continuous.missing",
                      "onehotencoder"))
  registerDoParallel(cl)
  score.matrix = foreach(j = 1:p, .combine = cbind,
                         .export = c("row_computing3", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                           lines = row_computing3(types, dend, original_data, tol, seed, j)
                           lines
                         }
  stopCluster(cl)
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}

L_score_4 = function(dend, original_data, tol  = 1e-20, seed = 99){
  types = sapply(original_data, class)
  p = length(original_data[1, ])
  n = length(original_data[, 1])
  total = p*n
  score.matrix = matrix(0, nrow = n, ncol = p)
  names = row.names(original_data)
  row.names(original_data)  = c(1:nrow(original_data))
  if(is.null(row.names(original_data))) row.names(original_data) = c(1:nrow(original_data))
  # paralellizing
  cores = detectCores()
  cl = makeCluster(cores[1] - 1)
  clusterExport(cl, c("row_computing", "factorial.missing", "continuous.missing",
                      "onehotencoder"))
  registerDoParallel(cl)
  score.matrix = foreach(j = 1:p, .combine = cbind,
                         .export = c("row_computing4", "factorial.missing","onehotencoder"),
                         .packages = c("ape", "phytools", "mvMORPH")) %dopar% {
                           lines = row_computing4(types, dend, original_data, tol, seed, j)
                           lines
                         }
  stopCluster(cl)
  partial_score = colSums(score.matrix)
  score = sum(partial_score)/total
  return(score)
}

L_cross_val = function(original_data, cl.method = "all", tol = 1e-20, seed = 99, method = NULL){
  p = ncol(original_data)
  if (cl.method == "all"){
    results = matrix(nrow = p, ncol = 10)
    for(k in (1:p)){
      test_data = as.data.frame(original_data[, k])
      colnames(test_data) = colnames(original_data)[k]
      rownames(test_data) = rownames(original_data)
      training_data = original_data[, -c(k)]
      types = sapply(training_data, class)
      bool =  (types == "integer" | types == "numeric")
      if (length(bool[bool != T]) == 0){ 
      d = dist(scale(training_data))
      hc.list = list(ward.hc = hclust(d, method = "ward.D"),
                     ward.d2.hc = hclust(d, method = "ward.D2"),
                     single.hc = hclust(d, method = "single"),
                     complete.hc = hclust(d, method = "complete"),
                     ave.hc = hclust(d, method = "average"),
                     med.hc = hclust(d, method = "median"),
                     centroid.hc = hclust(d, method = "centroid"),
                     mcquitty.hc = hclust(d, method = "mcquitty"),
                     agnes.hc = agnes(d, diss = TRUE),
                     diana.hc = diana(d, diss = TRUE))
      tops = lapply(hc.list, convert_to_phylo)
      results[k, ] = sapply(tops, FUN = L_score, original_data = test_data)
      }else{
        d = daisy(training_data, metric = "gower")
        hc.list = list(ward.hc = hclust(d, method = "ward.D"),
                       ward.d2.hc = hclust(d, method = "ward.D2"),
                       single.hc = hclust(d, method = "single"),
                       complete.hc = hclust(d, method = "complete"),
                       ave.hc = hclust(d, method = "average"),
                       med.hc = hclust(d, method = "median"),
                       centroid.hc = hclust(d, method = "centroid"),
                       mcquitty.hc = hclust(d, method = "mcquitty"),
                       agnes.hc = agnes(d, diss = T),
                       diana.hc = diana(d, diss = T))
        tops = lapply(hc.list, convert_to_phylo)
        print("a")
        results[k, ] = sapply(tops, FUN = L_score, original_data = test_data)
      }
    }
  }else{
    if(cl.method == "factors"){
      results = matrix(nrow = p, ncol = 8)
      for(k in (1:p)){
        test_data = as.data.frame(original_data[, k])
        colnames(test_data) = colnames(original_data)[k]
        rownames(test_data) = rownames(original_data)
        training_data = original_data[, -c(k)]
        types = sapply(training_data, class)
        bool =  (types == "integer" | types == "numeric")
        if (length(bool[bool != T]) == 0){ 
          d = dist(scale(training_data))
          hc.list = list(ward.hc = hclust(d, method = "ward.D"),
                         ward.d2.hc = hclust(d, method = "ward.D2"),
                         single.hc = hclust(d, method = "single"),
                         complete.hc = hclust(d, method = "complete"),
                         ave.hc = hclust(d, method = "average"),
                         mcquitty.hc = hclust(d, method = "mcquitty"),
                         agnes.hc = agnes(d, diss = TRUE),
                         diana.hc = diana(d, diss = TRUE))
          tops = lapply(hc.list, convert_to_phylo)
          results[k, ] = sapply(tops, FUN = L_score, original_data = test_data)
          print(k)
        }else{
          d = daisy(training_data, metric = "gower")
          hc.list = list(ward.hc = hclust(d, method = "ward.D"),
                         ward.d2.hc = hclust(d, method = "ward.D2"),
                         single.hc = hclust(d, method = "single"),
                         complete.hc = hclust(d, method = "complete"),
                         ave.hc = hclust(d, method = "average"),
                         mcquitty.hc = hclust(d, method = "mcquitty"),
                         agnes.hc = agnes(d, diss = T),
                         diana.hc = diana(d, diss = T))
          tops = lapply(hc.list, convert_to_phylo)
          results[k, ] = sapply(tops, FUN = L_score, original_data = test_data)
          print(k)
        }
      }
        }else{
    results = matrix(nrow = p, ncol = 8)
    for(k in (1:p)){
      test_data = as.data.frame(original_data[, k])
      colnames(test_data) = colnames(original_data)[k]
      rownames(test_data) = rownames(original_data)
      training_data = original_data[, -c(k)]
      types = sapply(training_data, class)
      bool =  (types == "integer" | types == "numeric")
      if (length(bool[bool != T]) == 0){ 
        d = dist(scale(training_data))
        hc.list = list(ward.hc = hclust(d, method = "ward.D"),
                       ward.d2.hc = hclust(d, method = "ward.D2"),
                       single.hc = hclust(d, method = "single"),
                       complete.hc = hclust(d, method = "complete"),
                       ave.hc = hclust(d, method = "average"),
                       med.hc = hclust(d, method = "median"),
                       centroid.hc = hclust(d, method = "centroid"),
                       mcquitty.hc = hclust(d, method = "mcquitty"))
        tops = lapply(hc.list, convert_to_phylo)
        results[k, ] = sapply(tops, FUN = L_score, original_data = test_data)
      }else{
        d = daisy(training_data, metric = "gower")
        hc.list = list(ward.hc = hclust(d, method = "ward.D"),
                       ward.d2.hc = hclust(d, method = "ward.D2"),
                       single.hc = hclust(d, method = "single"),
                       complete.hc = hclust(d, method = "complete"),
                       ave.hc = hclust(d, method = "average"),
                       med.hc = hclust(d, method = "median"),
                       centroid.hc = hclust(d, method = "centroid"),
                       mcquitty.hc = hclust(d, method = "mcquitty"))
        tops = lapply(hc.list, convert_to_phylo)
        results[k, ] = sapply(tops, FUN = L_score, original_data = test_data)
      }
      }
    }
  }
  return(results)
}

