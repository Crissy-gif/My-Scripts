################################################################################################
# Script for Distance Matrix for visualization mutiple datasets in tableau                      #
################################################################################################

rm(list = ls(all=TRUE))

setwd("~/Documents/Pelago_Project/Multiple Datasets")
source("~/Documents/Functions:Templates/Functions.R")

library(tidyr)
library(dplyr)
library(stringr)
library(tidyverse)
library(magrittr)
#-----------------load data/input file--------------------------------------------
file_set <- read.csv("2DDR_AbbVie_1803_14L.csv", sep = ";")
n <- read.csv("", sep = ";")
n <- read.csv("", sep = ";")

#--------------------------variable selection-------------------------------------

names(file_set)
Var_Sel <- file_set %>%
  select("Gene","pValue","Amplitude", "Protein.Name")

#----------variable manipulation---------------------------
#use regular expression if need to change the data format/structure e.g if the dataframe contains commas replace them with dot
#########################################################################################################################

#----------variable manipulation---------------------------
#use regular expression if need to change the data format/structure e.g if the dataframe contains commas replace them with dot

library(stringr)
Var_Sel$pValue2 %<>%
  #str_replace_all(",", ".") %>%
  #str_replace_all("×10\\^", "E") %>%
  as.numeric

Var_Sel$Amplitude2 <- Var_Sel$Amplitude %>%
  str_split(", ") %>%
  lapply(function(xlist){as.numeric(xlist)})

bad_words <- c("", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "of", "protein", "and", "of")
Var_Sel$words <- Var_Sel$Protein.Name %>%
  str_split(" |-|,|\\[|\\]|/|\\(|\\)") %>%
  lapply(FUN=function(l) {l[!(str_to_lower(l) %in% bad_words)]})

all_words <- do.call("c", Var_Sel$words)
all_words_table <- table(all_words)
tail(sort(all_words_table), 50)

#   as.numeric(stringr::str_replace_all(stringr::str_replace_all(Var_Sel$pValue, ",", "."), "×10\\^", "E"))

#--------calulate the square roots of the proteins-------------------
#rmse(Root mean square error, also called root mean square deviation)

protein_dist <- function(protein1, protein2){
  amp_length <- length(protein1$Amplitude)
  
  amp1 <- protein1$Amplitude2[[1]]
  amp2 <- protein2$Amplitude2[[1]]
  amp1[which(is.nan(amp1))] <- 0
  amp2[which(is.nan(amp2))] <- 0
  rmse <- sqrt(mean((amp1 - amp2)^2))
  
  #rmse_ctl_diff <- sqrt(mean(((protein1$amp - protein1$amp_ctl) - (protein2$amp - protein2$amp_ctl))^2))
  #total <- rmse + rmse_ctl_diff
  
  total <- rmse
  total
}

pr_list <- list(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, pr9, pr10)
protein_dist(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, pr9, pr10) # <- a bit tedious
protein_dist(pr_list[1], pr_list[2], ...) # <- even more tedious
do.call("protein_dist", pr_list) # <- beautiful, sexy, and cute!


#-------------create the distance matrix------------------
install.packages("doParallel")
library(doParallel)
library(foreach)
n_threads <- 10
cl <- makeCluster(n_threads) # 10 threads
doParallel::registerDoParallel(cl)
foreach::getDoParRegistered()


par_chunks <- function(n, n_batches){
  assertthat::assert_that(n >= n_batches)
  chunks <- seq(1,n,as.integer(n/n_batches))
  if(chunks[length(chunks)]!=n){
    chunks <- c(chunks, n)
  }
  chunks[length(chunks)] <- chunks[length(chunks)] + 1
  chunks
}
dist_matrix <- function(proteins, protein_dist_function=protein_dist, n_batches=20){
  chunks <- par_chunks(nrow(proteins), n_batches)
  i_rows <- foreach(chunk_idx = 1:(length(chunks)-1)) %dopar% {
    i_rows_chunk <- lapply(chunks[chunk_idx]:(chunks[chunk_idx+1]-1), function(i){
      i_row = vector()
      for(j in i:nrow(proteins)) {
        if(i==j){
          i_row[j] <- 0
          next
        }
        i_row[j] <- protein_dist_function(proteins[i,], proteins[j,])
      }
      i_row
    })
    do.call("rbind", i_rows_chunk)
  }
  uppertri <- do.call("rbind", i_rows)
  mat <- uppertri
  for(i in 1:nrow(proteins)){
    for(j in i:nrow(proteins)){
      mat[j,i] <- mat[i,j]
    }
  }
  mat
}

#-----------------####################-------------------------
my_dist_mat <- dist_matrix(head(Var_Sel, 10), n_batches = 2)
my_dist_mat

xy <- my_mds(my_dist_mat)

Var_Sel[1:200, 'mds_x'] <- xy$x
Var_Sel[1:200, 'mds_y'] <- xy$y
plot(xy$x,xy$y)

# Warning: takes about 15 minutes to run (before parallelization: 40):

ptm <- proc.time()
my_dist_mat_all <- dist_matrix(Var_Sel, n_batches = 20)
timediff <- proc.time() - ptm
print(timediff)

xy <- my_mds(my_dist_mat_all)
Var_Sel[, 'mds_x'] <- xy$x
Var_Sel[, 'mds_y'] <- xy$y

# Amplitudes and words - long format

Var_Sel_long <- Var_Sel %>%
  select(c("Gene","pValue", "Amplitude2","mds_x","mds_y", "words")) %>%
  mutate(ampl_ind=list(1:10)) %>%
  unnest(Amplitude2,ampl_ind,.preserve = words) %>%
  unnest(words) %>%
  #melt(id=c("Gene","Stab.Destab","pValue2","mds_x","mds_y")) %>%
  rename(amplitude_index=ampl_ind, amplitude=Amplitude2)

# A small sample to check that everything looks OK:
Var_Sel_long[sample.int(nrow(Var_Sel_long), size = 10),]

write.csv(Var_Sel_long[,c("Gene", "pValue" , "mds_x", "mds_y", "amplitude_index", "amplitude", "words")],"tableau_extract_long2.txt")

write.csv(Var_Sel[1:200,c("Gene", "Stab.Destab", "pValue2", "mds_x", "mds_y")], "tableau_extract_small.txt")

write.csv(Var_Sel_long, "AbbVie_14L_tableau.txt")



