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
       variables from wheat seeds dataset") +
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



