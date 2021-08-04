# code exploring USArrests using our score
library(ape)
library(dendextend)
library(cluster)
library(tibble)
library(tidyverse)
library(phytools)
library(mltools)
library(data.table)
library(factoextra)
source("C:/Users/lucru/Estatística_UFSCar/cv_cluster/modules/convert_to_parenthesis.R")
source("C:/Users/lucru/Estatística_UFSCar/cv_cluster/modules/cv_score.R")
library(tictoc)
library(mvMORPH)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(MASS)
library(kableExtra)
library(plyr)
library(clValid)

# importing USArrests
data("USArrests")
arr.data = USArrests
row.names(arr.data) = 1:nrow(arr.data)

# list of possible combinations:
# removing median and centroid linkage
test.list = list(hclust = c("ward.D", "single","ward.D2", "complete", "mcquitty", "average"), diana = NA)
dists = c("canberra", "euclidean", "manhattan")
tol = 1e-20

scores_arr = L_cross_val(arr.data, cl.list = test.list, dists = dists,
                         tol = tol, seed = 99, scale =  T)

id = which.min(scores_arr$V1)
min_id = scores_arr[id,]

# ordenando os scores
scores_arr %<>% mutate(names = row.names(scores_arr))
scores_arr %>% arrange(V1) %>% head(8)

arr_dend = hclust(dist(USArrests, method = "manhattan"), method = "mcquitty")
plot(arr_dend)


# selecting the best 
mcquitty.tree = convert_to_phylo(arr_dend)

# dendrogram for murder
par(mfrow = c(2,2))
x = setNames(arr.data[,1],gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]

obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1),outline=FALSE,lwd=c(3,7),leg.txt="Murder")

# dendrogram for assault
x = setNames(arr.data[,2],gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]

obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1),outline=FALSE,lwd=c(3,7),leg.txt="Assault")

# dendrogram for Urbanpop
x = setNames(arr.data[,3],gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]

obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1),outline=FALSE,lwd=c(3,7),leg.txt="Urbanpop")

# dendrogram for rape
x = setNames(arr.data[,4],gsub(" ", "", rownames(USArrests)))
reordered_x = x[mcquitty.tree$tip.label]

obj = contMap(mcquitty.tree, reordered_x, plot=FALSE)
obj = setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.7,1),outline=FALSE,lwd=c(3,7),leg.txt="Rape")

# importance by variable
test.list = list(hclust = c("mcquitty"))
dists = c("canberra", "euclidean", "manhattan")
tol = 1e-20
L_per_var = L_cross_val_per_var(arr.data, test.list, dists, scale = T)

ggplot_data = reshape2::melt(L_per_var[3, -5])
ggplot_data$variable = as.factor(colnames(arr.data))
p1 = ggplot_data %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Cross validated score values",
       title = "Scatterplot of cross validated score values versus 
       variables from USArrests dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0, max(ggplot_data$value) + 0.05)
p1


# correlation matrix
library(ggcorrplot)
ggcorrplot(cor(arr.data), hc.order = TRUE, type = "lower",
           outline.col = "white", lab = T) +
  ggsave(filename  = "corr_usa_arrests.pdf",
path = "C:/Users/lucru/Estatística_UFSCar/cv_cluster/figures",
width = 20.75, height = 12.5, units = "cm")


# using an alternative version of importance by variable
test.list = list(hclust = c("mcquitty"))
dists = c("canberra", "euclidean", "manhattan")
tol = 1e-20
L_per_var = L_cross_val_per_var_alt(arr.data, test.list, dists, scale = T)

ggplot_data = reshape2::melt(L_per_var[3, -5])
ggplot_data$variable = as.factor(colnames(arr.data))
p1 = ggplot_data %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Importance value",
       title = "Scatterplot of variable importance for USArrests dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0, max(ggplot_data$value) + 0.05)
p1


# adding and testing wheat seeds dataset
# wheat seeds dataset
wheat_data = read.delim("C:/Users/lucru/Estatística_UFSCar/cv_cluster/data/seeds_dataset.txt",
                        header = F)
wheat_data$V8 = as.factor(wheat_data$V8)
head(wheat_data)

# using an alternative version of importance by variable
test.list = list(hclust = c("mcquitty"))
dists = c("canberra", "euclidean", "manhattan")
tol = 1e-20
L_per_var_wheat = L_cross_val_per_var_alt(wheat_data[, -8], test.list, dists, scale = T)

ggplot_data = reshape2::melt(L_per_var_wheat[3, -8])
ggplot_data$variable = as.factor(colnames(wheat_data[,-8]))
p1 = ggplot_data %>%
  ggplot(aes(x = variable, y  = value, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Cross validated score values",
       title = "Scatterplot of cross validated score values versus 
       variables from wheat seeds dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5)) +
  ylim(0, max(ggplot_data$value) + 0.05)
p1


# variable alternative importance for simulated dataset and wheat seeds
# new test.list with more methods
test.list = list(hclust = c("ward.D", "single","ward.D2", "average", "complete", "mcquitty", "median"), 
                 agnes = c("weighted", "average", "ward"), diana = NA)

# wheat seeds dataset
wheat.l_cross_per_var = L_cross_val_per_var_alt(wheat_data[, -8], test.list, dists, scale = T)

# simulated dataset
data.generator = function(n, mu_0, mu_1, p = 2, S = NULL, seed = 500){
  set.seed(seed)
  Y = sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  if(is.null(S) == T) S = diag(nrow = p, ncol = p)
  X = matrix(nrow = n, ncol = p)
  X[Y == 0, ] = round(mvrnorm(sum(Y == 0), mu_0, S), 4)
  X[Y == 1, ] = round(mvrnorm(sum(Y == 1), mu_1, S), 4)
  
  sim.data = as.data.frame(cbind(X , Y))
  sim.data$Y = as.factor(Y)
  if(p == 2){
    colnames(sim.data) = c("X1", "X2", "Y")}
  row.names(sim.data) = c(1:nrow(sim.data))
  return(sim.data)
}


n = 170
p = 2
mu_0 = c(0, 0)
mu_1 = c(2, 2)
S = cbind(c(2, 1), c(1, 3))
sim.data = data.generator(n, mu_0, mu_1, S = S)

sim.data.l_cross_per_var_scaled = L_cross_val_per_var_alt(sim.data[, -3], test.list, dists, scale = T)


nb.cols = 33
mycolors = colorRampPalette(brewer.pal(33, "Paired"))(nb.cols)

ggplot_data = reshape2::melt(wheat.l_cross_per_var[, -8])
ggplot_data$dists = rep(wheat.l_cross_per_var$V8, 7)
ggplot_data$obs = rep(as.factor(1:33), 7)

p1 = ggplot_data %>%
  mutate(variable = as.factor(variable),
         dists = as.factor(dists)) %>%
  ggplot(aes(x = variable, y  = value, group = obs)) +
  geom_line(aes(color = obs)) +
  geom_point(aes(color = obs)) +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Importance score values",
       colour = "Combinations",
       title = "Scatterplot of variable importance for 
       wheat seeds dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
       plot.title = element_text(hjust = 0.5))

p1
ggplot_data = reshape2::melt(sim.data.l_cross_per_var_scaled[, -3])
ggplot_data$dists = rep(sim.data.l_cross_per_var_scaled$V3, 2)
ggplot_data$obs = rep(as.factor(1:33), 2)

p2 = ggplot_data %>%
  mutate(variable = as.factor(variable),
         dists = as.factor(dists)) %>%
  ggplot(aes(x = variable, y  = value, group = obs)) +
  geom_line(aes(color = obs)) +
  geom_point(aes(color = obs)) +
  scale_colour_manual(values = mycolors) +
  theme_minimal() +
  labs(x = "Variables names",
       y = "Importance score values",
       colour = "Combinations",
       title = "Scatterplot of variable importance for
       simulated dataset") +
  theme(text = element_text(size = 11, 
                            family ="serif"),
        plot.title = element_text(hjust = 0.5))
p2


library(ggpubr)
ggarrange(p1, p2, common.legend = T, nrow = 2, legend = "right")+
  ggsave(filename  = "wheat_seeds_simulated_importances.pdf",
         path = "C:/Users/lucru/Estatística_UFSCar/cv_cluster/figures",
         width = 22.75, height = 14.5, units = "cm")


# adding the phylogenetic simulated example
# simulating tree
set.seed(1234)
sim_tree = rtree(n = 24, tip.label = letters[1:24])
plotTree(sim_tree)

# generating discrete states:
disc_states = rTraitDisc(sim_tree, k = 2, freq = 0.45, rate = 0.65)

# plotting:
fitER<-ace(disc_states, sim_tree, model="ARD", type="discrete")
fitER
cols<-setNames(c("red","blue"),levels(disc_states))


plotTree(sim_tree, type="fan",fsize=0.9,ftype="i",lwd=1)
nodelabels(node=1:sim_tree$Nnode+Ntip(sim_tree),
           pie=fitER$lik.anc,
           piecol=cols,
           cex=0.4)
tiplabels(pie=to.matrix(disc_states[sim_tree$tip.label],
                        levels(disc_states)),piecol=cols,cex=0.25)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)

# generating continuous
cont_states = rTraitCont(sim_tree, sigma = 2)


obj<-contMap(sim_tree, cont_states, plot=FALSE)
plot(obj,lwd=7)




