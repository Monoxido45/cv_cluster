# code used to convert a tree generated by clustering methods into parenthetic format

library(ape)
library(dendextend)
library(cluster)
library(tibble)
library(magrittr)
library(dplyr)
library(phytools)
library(stats)

# function converting generical clustering object into dendrogram object
to.dend = function(cl.obj){
  # converting generical clustering object into hclust object first
  cluster.obj = cl.obj %>% as.hclust()
  dend = cluster.obj %>% as.dendrogram(hang = -1, check = TRUE)
  return(dend)
}

# formating the obtained dendrogram into parentical format
convert_to_par <- function(dend, first_it =  TRUE)
{
  if (first_it == TRUE){
    dist = as.double(attr(dend, "height"))
    dend.object1 = dend[[1]]
    dist1 = as.double(attr(dend.object1, "height"))
    dend.object2 = dend[[2]]
    dist2 = as.double(attr(dend.object2, "height"))
    first_it = FALSE
    if (is.list(dend.object1) == TRUE & is.list(dend.object2) == TRUE){
    return(paste0("((",
                  convert_to_par(dend.object1, first_it),
                  "):",
                  dist - dist1,
                  ",",
                  "(",
                  convert_to_par(dend.object2, first_it),
                  "):",
                  dist - dist2,
                  ");"))}else{
                    if(is.list(dend.object1) == TRUE & is.list(dend.object2) == FALSE){
                      label = attr(dend.object2, "label")
                      return(paste0("((",
                                    convert_to_par(dend.object1, first_it),
                                    "):",
                                    dist - dist1,
                                    ",",
                                    label,
                                    ":",
                                    dist - dist2,
                                    ");"))
                    }else{
                      label = attr(dend.object1, "label")
                      return(paste0("(",
                                    label,
                                    ":",
                                    dist - dist1,
                                    ",",
                                    "(",
                                    convert_to_par(dend.object2, first_it),
                                    "):",
                                    dist - dist2,
                                    ");"))
                    }
                  }
  }else{
      dist = as.double(attr(dend, "height"))
      dend.object1 =  dend[[1]]
      dist1 = as.double(attr(dend.object1, "height"))
      dend.object2 = dend[[2]]
      dist2 = as.double(attr(dend.object2, "height"))
      if(is.list(dend.object1) == TRUE & is.list(dend.object2) == TRUE){
      return(paste0("(",
                    convert_to_par(dend.object1, first_it),
                    "):", dist - dist1,
                    ",",
                    "(",
                    convert_to_par(dend.object2, first_it),
                    "):", dist - dist2))
      }else{
        if(is.list(dend.object1) == FALSE & is.list(dend.object2) == TRUE){
          label1 = attr(dend.object1,"label")
          return(paste0(label1, ":",
                        dist - dist1,
                        ",",
                        "(",
                        convert_to_par(dend.object2, first_it),
                        "):", dist - dist2))
        }else{
          if(is.list(dend.object1) == TRUE & is.list(dend.object2) == FALSE){
            label2 = attr(dend.object2, "label")
            return(paste0("(",
                          convert_to_par(dend.object1, first_it),
                          "):", dist - dist1, ",",
                          label2, ":",
                          dist - dist2))
          }else{
            label1 = attr(dend.object1, "label")
            label2 = attr(dend.object2, "label")
            return(paste0(label1, ":", dist - dist1,
                          ",",
                          label2, ":", dist - dist2))
          }
        }
      }
    }
}


convert_to_phylo = function(cl.obj){
  dend = to.dend(cl.obj)
  dend.str = convert_to_par(dend)
  tree = read.tree(text = dend.str)
  return(tree)
}



