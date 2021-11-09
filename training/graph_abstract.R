# code to make graphical abstract and boxplot/evolutionary dendrogram comparisson
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

# making simulated dataset for graphical abstract
# simulating dataset from multivariate normal
# according to a latent variable with 3 levels
set.seed(125)
n = 15
y = sample(c(1, 2, 3), size = n, replace = TRUE, prob = c(1/3, 1/3, 1/3))
mu_1 = c(1, 1, 1, 3)
mu_2 = c(4, 4, 1, 3)
mu_3 = c(-2, -1, 2, 2)
S = rbind(c(1, 0.5, 0.15, 0.15), c(0.5, 1, 0.15, 0.15), 
          c(0.15, 0.15, 1, 0.5), c(0.15, 0.15, 0.5, 1))
X = ((y == 1)*(mvrnorm(15, mu_1, S)) + (y == 2)*(mvrnorm(15, mu_2, S)) +
  (y == 3)*(mvrnorm(15, mu_2, S))) %>% 
  as.data.frame()
colnames(X) = c("X1", "X2", "X3", "X4")

# generating 3 different dendrograms with euclidean distance
ward_tree = hclust(dist(scale(X)), method = "ward.D2") %>%
  convert_to_phylo()
comp_tree = hclust(dist(scale(X)), method = "complete") %>%
  convert_to_phylo()
mcquitty_tree = hclust(dist(scale(X)), method = "mcquitty") %>%
  convert_to_phylo()


par(mfrow = c(1, 3))
plot(ward_tree,cex = 1)

plot(comp_tree,cex = 1)

plot(mcquitty_tree, cex = 1)

test.list = list(hclust = c("ward.D2", "complete","mcquitty"))
dists = c("euclidean")

tol = 1e-20

scores_arr = L_cross_val(X, cl.list = test.list, dists = dists,
                         tol = tol, seed = 99, scale =  T)

# importance score
test.list = list(hclust = c("mcquitty"))
dists = c("euclidean")
import_x = L_cross_val_per_var_alt(X, test.list, 
                        dists, scale = T)

ggplot_data = reshape2::melt(import_x)
ggplot_data$variable = as.factor(colnames(X))
p1 = ggplot_data %>%
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
  ylim(min(ggplot_data$value) - 0.05, max(ggplot_data$value) + 0.05)

p1


# evolutionary dendrograms for these two cases:
par(mfrow = c(1, 2))

# ploting X2 first
x = setNames(X$X2,gsub(" ", "", rownames(X)))
reordered_x = x[mcquitty_tree$tip.label]

obj = contMap(mcquitty_tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize= 0.7,outline=FALSE, lwd = c(2,5), leg.txt = "X2")

# X1 now
x = setNames(X$X1,gsub(" ", "", row.names(X)))
reordered_x = x[mcquitty_tree$tip.label]

obj = contMap(mcquitty_tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize= 0.7,outline=FALSE, lwd = c(2,5), leg.txt = "X1")


# generating boxplots for USArrest according to dendrogram partitions
arr_dend = hclust(dist(scale(USArrests), method = "manhattan"), method = "mcquitty")
# scaling and secting several partitions
arr_data = USArrests
arr_data %<>%
  scale() %<>%
  as.data.frame() %<>%
  mutate("k = 2" = as.factor(cutree(arr_dend, k = 2)),
         "k = 3" = as.factor(cutree(arr_dend, k = 3)),
         "k = 4" = as.factor(cutree(arr_dend, k = 4)))


melted_data = arr_data %>%
  pivot_longer(Murder:Rape, names_to = "variable", values_to = "values") %>%
  pivot_longer("k = 2":"k = 4", names_to = "k", values_to = "cluster")


melted_data %>%
  ggplot(aes(y = variable, x = values, fill = cluster)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Standardized values",
       y = "Variables",
       fill = "Cluster") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~k)

p1 = ggplot_data %>%
  ggplot(aes(y = variable, x = value, fill = clust_2)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Values",
       y = "Variables",
       fill = "Cluster",
       title = "k = 2") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette = "Set1")

p1
p2 = ggplot_data %>%
  ggplot(aes(y = variable, x = value, fill = clust_3)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Values",
       y = "Variables",
       fill = "Cluster",
       title = "k = 3") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette = "Set1")
p2

p3 = ggplot_data %>%
  ggplot(aes(y = variable, x = value, fill = clust_4)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Values",
       y = "Variables",
       fill = "Cluster",
       title = "k = 4") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_brewer(palette = "Set1")

p3

library(ggpubr)
ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = T)

