# simulating data from bayesian network model
set.seed(2075)
d <- 6
k <- 8
cov_stab <- seq(from = 0.5, to = 1.5, length.out = 6)
evolve <- function(obs, time) 
  obs + rnorm(d, 0, cov_stab^time)

param <- list()
param[[1]] <- matrix(0, nrow = 1, ncol = d)
for(time in 1:(k-1))
{
  param[[time + 1]] <- matrix(NA, nrow = 2^time, ncol = d)
  n_param_ant = 2^(time-1)
  for(jj in 1:n_param_ant)
  {
    mut_1 = evolve(param[[time]][jj,], time)
    mut_2 = evolve(param[[time]][jj,], time)
    param[[time + 1]][jj,] = mut_1
    param[[time + 1]][jj + n_param_ant,] = mut_2
  }
}


# importing packages
library(ape)
library(dendextend)
library(cluster)
library(tidyverse)
library(phytools)
library(mltools)
library(data.table)
library(factoextra)
setwd("~/estatistica_UFSCAR/cv_cluster/modules")
source("convert_to_parenthesis.R")
source("cv_score.R")
library(tictoc)
library(mvMORPH)
library(RColorBrewer)
library(MASS)
library(plyr)
library(clValid)

sim_data = param[[k]] %>%
  as.data.frame()

test.list = list(hclust = c("ward.D2"))
dists = c("euclidean")
import_sim = L_cross_val_per_var_alt(sim_data, test.list, 
                                   dists, scale = T)

melted_sim = reshape2::melt(import_sim)
melted_sim$variable = as.factor(colnames(import_sim))
p1 = melted_sim %>%
  mutate(variable = fct_reorder(variable, value, .desc = T)) %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Features",
       y = "Importance score") +
  theme(text = element_text(size = 14, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(min(melted_sim$value) - 0.05, max(melted_sim$value) + 0.05)

p1

# RF
dend_ward = sim_data %>%
  scale() %>%
  dist() %>%
  hclust(method = "ward.D2")


sim_data %<>%
  mutate(clust3 = as.factor(cutree(dend_ward, k = 3)),
         clust6 = as.factor(cutree(dend_ward, k = 6)))

rf_clust_3 = ranger::ranger(formula = clust3 ~ . - clust6,
                            data = sim_data,
                            num.trees = 500,
                            importance = "impurity",
                            write.forest = TRUE,
                            verbose = FALSE,
                            probability = T)

rf_clust_6 = ranger::ranger(formula = clust6 ~ . - clust3,
                            data = sim_data,
                            num.trees = 500,
                            importance = "impurity",
                            write.forest = TRUE,
                            verbose = FALSE,
                            probability = T)

import = tibble(variable = c(names(ranger::importance(rf_clust_3)), 
                             names(ranger::importance(rf_clust_6))),
                importance = c(ranger::importance(rf_clust_3),
                               ranger::importance(rf_clust_6))) %>%
  mutate(cluster = as.factor(rep(c("k = 3", "k = 6"), each = dim(sim_data)[2] - 2))) %>%
  arrange(desc(importance))

import %>%
  ggplot(aes(x = variable,
             y = importance))+
  geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3") + 
  coord_flip()+
  labs(y = "Importance score",
       x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 15,
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient(low = "firebrick2", high = "dodgerblue3")+
  facet_wrap(~cluster)

dend_ward %<>% convert_to_phylo()
par(mfrow= c(2, 3))
x = setNames(round(sim_data[, 1], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="X1", ftype = "off")


x = setNames(round(sim_data[, 2], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="X2", ftype = "off")


x = setNames(round(sim_data[, 3], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="X3", ftype = "off")


x = setNames(round(sim_data[, 4], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="X4", ftype = "off")

x = setNames(round(sim_data[, 5], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="X5", ftype = "off")


x = setNames(round(sim_data[, 6], 2), gsub(" ", "", rownames(sim_data)))
reordered_x = x[dend_ward$tip.label]

obj = contMap(dend_ward, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1), outline=FALSE, lwd = c(2,5), leg.txt="X6", ftype = "off")

nb = NbClust::NbClust(sim_data[, -c(7, 8)], distance = "euclidean",
                 min.nc = 2, max.nc = 8, method = "ward.D2")

factoextra::fviz_nbclust(nb)
