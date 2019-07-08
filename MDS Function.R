
# Functions for similarity distance

#--------calulate the square roots of the proteins-------------------
#rmse(Root mean square error, also called root mean square deviation)

protein_dist <- function(protein1, protein2){
  amp1 <- protein1$Amplitude2[[1]]
  amp2 <- protein2$Amplitude2[[1]]
  amp1[which(is.nan(amp1))] <- 0
  amp2[which(is.nan(amp2))] <- 0
  rmse <- sqrt(mean((amp1 - amp2)^2))
  rmse
}
#----------calutate the square roots of proteins in 3dim, 5 concs in 10 temps---------
protein_3dims_dist <- function(protein1, protein2){
  
  amp1 <- protein1$Amplitude2[[1]]
  amp2 <- protein2$Amplitude2[[1]]
  amp1[which(is.nan(amp1))] <- 0
  amp2[which(is.nan(amp2))] <- 0
  rmse <- sqrt(mean((amp1 - amp2)^2))
  rmse
}



#------------------Protein dist in 3d but with temp that has NAN in them--------------

protein_dist_3d_na2 <- function(protein1, protein2){
  amp_length <- length(protein1$T_combo)
  amp1 <- protein1$T_combo[[1]]
  amp2 <- protein2$T_combo[[1]]
  nans <- is.nan(amp1) | is.nan(amp2) # All indices with nans in either protein1 or protein 2
  rmse <- sqrt(mean((amp1[!nans] - amp2[!nans])^2))
  rmse
}

#-------------create the distance matrix------------------
dist_matrix <- function(proteins){
  mat <- matrix(nrow = nrow(proteins), ncol = nrow(proteins))
  for(i in 1:(nrow(proteins))){
    for(j in i:nrow(proteins)){
      if(i==j){
        mat[i,j] <- 0
      }
      mat[i,j] <- protein_dist(proteins[i,], proteins[j,])
      mat[j,i] <- protein_dist(proteins[i,], proteins[j,])
    }
  }
  mat
}


#-------------------Creating chunks----------------------------------------
par_chunks <- function(n, n_batches){
  assertthat::assert_that(n >= n_batches)
  chunks <- seq(1,n,as.integer(n/n_batches))
  if(chunks[length(chunks)]!=n){
    chunks <- c(chunks, n)
  }
  chunks[length(chunks)] <- chunks[length(chunks)] + 1
  chunks
}

#--------------------Creating distance matrix---------------------------

dist_matrix <- function(proteins, protein_dist_function=protein_dist, n_batches=20){
  chunks <- par_chunks(nrow(proteins), n_batches)
  i_rows <- foreach(chunk_idx = 1:(length(chunks)-1), .packages = "tcltk") %dopar% {
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


#---------------------

my_mds <- function(mat){
  fit <- cmdscale(mat, eig = TRUE, k = 2)
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  list(x=x,y=y)
}


