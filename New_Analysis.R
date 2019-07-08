###Analysis for the new updated genes frome zebrafish genome

rm(list=ls(all=TRUE))
setwd("~/Desktop/Master*s project/code")
#data_file contains all the varibles separated by tabs
data_file <- read.csv("for_canonical_analysis.txt", sep="\t", row.names = "barcode")

####
###variable selection of phenotype data from data_file
library(dplyr)
pheno_variables <- data_file %>%
  select(intldl_Res_r:MacNeu_Area_r)

###Variable selection of the genotype data from data_file

geno_variables <-data_file %>%
  select(D_apobb1_133353:D_ldlra_046297)
  
########################Canonical functions modified#####################################
#####Cancor2 functions combination of the yacca package and CCA package

cancor2 <- function(X, Y, xcols=NA, ycols=NA, subset_selection=TRUE){
  #####
  # Returns NULL if X or Y is rank deficient
  ####
  
  library(yacca)
  if(identical(ycols, NA)){
    ycols=names(Y)
  }
  if(identical(xcols, NA)){
    xcols=names(X)
  }
  library(tidyr) 
  library(dplyr)
  XY <- X %>%
    merge(Y, by="row.names") %>%
    select(ycols, xcols) %>%
    drop_na
  qrx <- qr(XY[xcols]) # QR decomposition with column pivoting
  qry <- qr(XY[ycols])
  if(subset_selection){
    selected_xcols <- colnames(qrx$qr)[1:qrx$rank]
    selected_ycols <- colnames(qry$qr)[1:qry$rank]
  } else {
    if(qrx$rank != length(xcols) || qry$rank != length(ycols)) {
      return(NULL)  
    } else  {
      selected_xcols <- xcols
      selected_ycols <- ycols
    }
  }
  #cca1 <- cancor(XY[xcols], XY[ycols])
  #selected_xcols <- row.names(cca1$xcoef)
  #selected_ycols <- row.names(cca1$ycoef)
  #cca2 <- yacca::cca(XY[selected_xcols], XY[selected_ycols])
  cca2 <- yaccacca2(XY[selected_xcols], XY[selected_ycols])
  fcca2 <- yacca::F.test.cca(cca2)
  list(cca2, fcca2)
  #sink(cca1,"CCA1.txt")
}

#### Modified yacca function
uncomplex <- function(m){
  oldrownames <- row.names(m)
  newstuff <- as.matrix(apply(m, 2, as.numeric))
  row.names(newstuff) <- oldrownames
  newstuff
}

yaccacca2 <- function (x, y, xlab = colnames(x), ylab = colnames(y), xcenter = TRUE, 
                       ycenter = TRUE, xscale = FALSE, yscale = FALSE, standardize.scores = TRUE, 
                       use = "complete.obs", na.rm = TRUE) 
{
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (is.data.frame(y)) 
    y <- as.matrix(y)
  if (xcenter || xscale) 
    x <- scale(x, scale = xscale, center = xcenter)
  if (ycenter || yscale) 
    y <- scale(y, scale = yscale, center = ycenter)
  if (is.null(dim(x))) 
    x <- matrix(x, ncol = 1)
  if (is.null(dim(y))) 
    y <- matrix(y, ncol = 1)
  nx <- dim(x)[2]
  ny <- dim(y)[2]
  ncv <- min(nx, ny)
  cvlab <- paste("CV", 1:ncv)
  o <- list()
  cxx <- cov(x, use = use)
  cyy <- cov(y, use = use)
  cxy <- cov(x, y, use = use)
  cyx <- t(cxy)
  ey <- eigen(qr.solve(cyy, cyx) %*% qr.solve(cxx, cxy))
  ex <- list(values = ey$values, vectors = qr.solve(cxx, cxy) %*% 
               (ey$vec))
  o$corr <- (ex$val[1:ncv])^0.5
  # Custom code starts here
  o$corr <- as.numeric(o$corr)
  # Custom code ends here
  names(o$corr) <- cvlab
  o$corrsq <- o$corr^2
  names(o$corrsq) <- cvlab
  o$xcoef <- ex$vec[, 1:ncv, drop = FALSE]
  rownames(o$xcoef) <- xlab
  colnames(o$xcoef) <- cvlab
  # Custom code starts here
  o$xcoef <- uncomplex(o$xcoef)
  # Custom code ends here
  o$ycoef <- ey$vec[, 1:ncv, drop = FALSE]
  rownames(o$ycoef) <- ylab
  colnames(o$ycoef) <- cvlab
  # Custom code starts here
  o$ycoef <- uncomplex(o$ycoef)
  # Custom code ends here
  o$canvarx <- x %*% o$xcoef
  rownames(o$canvarx) <- rownames(x)
  colnames(o$canvarx) <- cvlab
  o$canvary <- y %*% o$ycoef
  rownames(o$canvary) <- rownames(y)
  colnames(o$canvary) <- cvlab
  if (standardize.scores) {
    sdx <- apply(o$canvarx, 2, sd)
    sdy <- apply(o$canvary, 2, sd)
    o$canvarx <- sweep(o$canvarx, 2, sdx, "/")
    o$canvary <- sweep(o$canvary, 2, sdy, "/")
    o$xcoef <- sweep(o$xcoef, 2, sdx, "/")
    o$ycoef <- sweep(o$ycoef, 2, sdy, "/")
  }
  # Custom code starts here
  o$canvarx <- uncomplex(o$canvarx)
  o$canvary <- uncomplex(o$canvary)
  # Custom code ends here
  o$xstructcorr <- cor(x, o$canvarx, use = use)
  rownames(o$xstructcorr) <- xlab
  colnames(o$xstructcorr) <- cvlab
  o$ystructcorr <- cor(y, o$canvary, use = use)
  rownames(o$ystructcorr) <- ylab
  colnames(o$ystructcorr) <- cvlab
  o$xstructcorrsq <- o$xstructcorr^2
  rownames(o$xstructcorrsq) <- xlab
  colnames(o$xstructcorrsq) <- cvlab
  o$ystructcorrsq <- o$ystructcorr^2
  rownames(o$ystructcorrsq) <- ylab
  colnames(o$ystructcorrsq) <- cvlab
  o$xcrosscorr <- cor(x, o$canvary, use = use)
  rownames(o$xcrosscorr) <- xlab
  colnames(o$xcrosscorr) <- cvlab
  o$ycrosscorr <- cor(y, o$canvarx, use = use)
  rownames(o$ycrosscorr) <- ylab
  colnames(o$ycrosscorr) <- cvlab
  o$xcrosscorrsq <- o$xcrosscorr^2
  rownames(o$xcrosscorrsq) <- xlab
  colnames(o$xcrosscorrsq) <- cvlab
  o$ycrosscorrsq <- o$ycrosscorr^2
  rownames(o$ycrosscorrsq) <- ylab
  colnames(o$ycrosscorrsq) <- cvlab
  o$xcancom <- apply(o$xstructcorrsq, 1, sum, na.rm = na.rm)
  names(o$xcancom) <- xlab
  o$ycancom <- apply(o$ystructcorrsq, 1, sum, na.rm = na.rm)
  names(o$ycancom) <- ylab
  o$xcanvad <- apply(o$xstructcorrsq, 2, mean, na.rm = na.rm)
  names(o$xcanvad) <- cvlab
  o$ycanvad <- apply(o$ystructcorrsq, 2, mean, na.rm = na.rm)
  names(o$ycanvad) <- cvlab
  o$xvrd <- o$xcanvad * o$corrsq
  names(o$xvrd) <- cvlab
  o$yvrd <- o$ycanvad * o$corrsq
  names(o$yvrd) <- cvlab
  o$xrd <- sum(o$xvrd, na.rm = na.rm)
  o$yrd <- sum(o$yvrd, na.rm = na.rm)
  bartvbase <- -(NROW(x) - 1 - (nx + ny + 1)/2)
  o$chisq <- bartvbase * (sum(log(1 - o$corr^2)) - c(0, cumsum(log(1 - 
                                                                     o$corr^2))[-ncv]))
  o$df <- (nx + 1 - (1:ncv)) * (ny + 1 - (1:ncv))
  names(o$chisq) <- cvlab
  names(o$df) <- cvlab
  o$xlab <- xlab
  o$ylab <- ylab
  class(o) <- "cca"
  o
}


####Multiple geno vs multiple pheno CCA
CCA_M_geno_vs_M_pheno <- cancor2(X=geno_variables, Y=pheno_variables)
library(CCA)
######################################################################
#Evaluation of the 2 methods!
######################################################################
# Are there any significant combinations?
#install.packages("doParallel")

find_significant_associations <- function(phenos, genos, alpha, geno_data, pheno_data, all_data){
  alpha_adjusted <- alpha/(length(genos)*length(phenos))# Set to 0.05, remove the B.C
  print(paste0("Adjusted alpha: ", alpha_adjusted))
  library(ggm)
  pset_geno <- powerset(set = genos)
  pset_geno2 <- pset_geno[sapply(pset_geno, function(x) {length(x) > 1})]
  print(paste0("Number of genos to test: ", length(pset_geno2)))
  pset_pheno <- powerset(set = phenos)
  pset_pheno2 <- pset_pheno[sapply(pset_pheno, function(x) {length(x) > 1})]
  print(paste0("Number of phenos to test: ", length(pset_pheno2)))
  total_iters <- length(pset_pheno2)*length(pset_geno2)
  print(paste0("Total number of combinations to test: ", total_iters))
  degfr <- length(genos)*length(phenos)
  
  result <- data.frame(phenos=numeric(total_iters), genos=numeric(total_iters),
      p_adjusted_CCA=numeric(total_iters), p_adjusted_LM=numeric(total_iters))
  result$phenos <- list(rep(list(),total_iters))
  result$genos <- list(rep(list(),total_iters))
  result$phenos <- NA
  result$genos <- NA
  
  sign_LM <- vector()
  sign_CCA <- vector()
  n_iter <- 0
  n_non_deficient <- 0
  min_LM_p <- Inf
  min_CCA_p <- Inf
  #library(doParallel)
  #registerDoParallel(cores=2)
  #foreach(geno_set=pset_geno2) %dopar% {
  for(geno_set in pset_geno2){
    for(pheno_set in pset_pheno2){
      n_iter <- n_iter + 1
      result$phenos[[n_iter]] <- pheno_set
      result$genos[[n_iter]] <- geno_set
      if(n_iter %% 1000 == 0){
        print(paste0("Iteration: ", n_iter, ". % complete: ", format(100*n_iter/total_iters, digits=4), 
                     ". Min CCA p-value: ", format(min_CCA_p, digits=6), ", adjusted: ", format(min_CCA_p*degfr, digits=6),
                     ". Min LM p-value: ", format(min_LM_p, digits=6), ", adjusted: ", format(min_LM_p*degfr, digits=6)))
      }
      # CCA
      res <- cancor2(geno_data, pheno_data, xcols = geno_set, ycols=pheno_set, subset_selection = FALSE)
      if(is.null(res)){
        next
      }
      n_non_deficient <- n_non_deficient + 1
      f_res <- res[[2]]
      p_value <- f_res$p.value[1]
      if(is.nan(p_value)){
        p_value <- Inf # Bug fix
      }
      if(p_value < min_CCA_p){
        min_CCA_p <- p_value
        min_p_CCA_model <- res
      }
      is_significant <- p_value < alpha_adjusted
      sign_CCA <- c(sign_CCA, is_significant)
      if(is_significant){
        print("Genos: ");print(geno_set);print("Phenos: ");
        print(pheno_set);print("-------------")
      }
      result[n_iter, "p_adjusted_CCA"] <- p_value*degfr
      
      # LM
      any_significant <- FALSE
      for(geno in geno_set){
        b <- paste(pheno_set, collapse = "+")
        lm_formula <- as.formula(paste0(geno, " ~ ", b))# Rethink
        tg_geno <- lm(formula = lm_formula, data = all_data)
        sm1 <- summary(tg_geno)
        p_value <- pf(sm1$fstatistic[1],sm1$fstatistic[2],sm1$fstatistic[3],lower.tail=FALSE)
        if(p_value < alpha_adjusted) {any_significant <- TRUE}
        if(p_value < min_LM_p){
          min_LM_p <- p_value
          min_p_LM_model <- tg_geno
        }
      }
      sign_LM <- c(sign_LM, any_significant)
      result[n_iter, "p_adjusted_LM"] <- p_value*degfr
    }
  }
  print(paste0("Number of iterations: ", n_iter))
  print(paste0("Number of non-deficient (not skipped) iterations: ", n_non_deficient))
  print(paste0("Minimum p-value CCA: ", min_CCA_p))
  print(paste0("Minimum p-value LM: ", min_LM_p))
  print(paste0("Minimum adjusted p-value CCA: ", min_CCA_p*degfr))
  print(paste0("Minimum adjusted p-value LM: ", min_LM_p*degfr))
  total_LM <- sum(sign_LM)
  total_CCA <- sum(sign_CCA)
  print(paste0("Total significant sets for LM: ", total_LM))
  print(paste0("Total significant sets for CCA: ", total_CCA))
  print(paste0("N where both significant:", sum(sign_LM[sign_CCA]))) # Both are significant
  print(paste0("N where only one is significant:", sum(sign_CCA + sign_LM == 1))) # Only one is significant
  print(paste0("N where none is significant:", sum(sign_CCA + sign_LM == 0))) # None is significant
  print(paste0("N where only CCA significant:", sum(sign_CCA & !sign_LM))) # CCA is significant but LM is not
  print(paste0("N where only LM significant:", sum((!sign_CCA) & sign_LM))) # LM is significant but CCA is not
  result
  #list(min_p_CCA_model=min_p_CCA_model, min_p_LM_model=min_p_LM_model, n_iters=n_iter)
}

phenos <- c("intldl_Res_r", "inthdl_Res_r", "inttc_Res_r", "inttg_Res_r", "intglucose_Res_r", 
            "intMac_Area_r", "intNeu_Area_r", "cld_r", "MacLip_Area_r", "NeuLip_Area_r", "MacNeu_Area_r")
genos <- c("D_apobb1_133353", "D_apoea_172219" , "D_apoeb_058965", "D_ldlra_046297")
#alpha <- 0.05
alpha <- 0.1
geno_data <- data_file[genos]
pheno_data <- data_file[phenos]
all_data <- data_file

significants <- find_significant_associations(phenos, genos, alpha, geno_data, pheno_data, all_data)

# It was significant only on the 0.1 level. Let's investigate the significant association.

significants <- significants %>%
  mutate(p_ratio = p_adjusted_LM/p_adjusted_CCA)
hist(significants$p_ratio)
max_p_ratio <- significants[which.max(significants$p_ratio),c("genos", "phenos")]
max_p_ratio_cca <- cancor2(geno_data, pheno_data, 
                           xcols = max_p_ratio$genos[[1]], ycols=max_p_ratio$phenos[[1]],
                           subset_selection = FALSE)
plot(max_p_ratio_cca[[1]])

s2 <- significants[significants$p_adjusted_CCA < 0.15,]
head(s2[sort.list(s2$p_ratio, decreasing = T),])
 
nrow(significants)

#plot(significants$min_p_CCA_model[[1]])
#plot(significants$min_p_LM_model)
#summary(significants$min_p_LM_model)
#significants$min_p_LM_model$terms

####################### Analysis with homozygous larva data####################################
library(dplyr)
homo_data <- filter(data_file, homozygous == 1)
geno_homo_var <- homo_data %>%
  select(D_apobb1_133353:D_ldlra_046297)
pheno_homo_var <-homo_data %>%
  select(intldl_Res_r:MacNeu_Area_r)

######## Combinations for the homo data#########################################################
geno_data <- geno_homo_var[genos]
pheno_data <- pheno_homo_var[phenos]
all_data <- homo_data

significants_homo <- find_significant_associations(phenos, genos, alpha, geno_data, pheno_data, all_data)


###########Testing with the old data#####################################################################################
data1 <- read.csv("LDLR_APOE_APOBb1_for_analysis.txt", sep="\t", row.names = "barcode")
geno <-data1 %>%
  select(apoba:apoe_ldlr, -ends_with("_homo"), -miss, -apoe_ldlr, - apob_ldlr, -apoe_apob)
pheno_orig <- data1 %>%
  select(cld:intvolume_Res, ends_with("_Res"), -line, -sequencing, 
         -starts_with("Annotation"), -mutant, starts_with("int"), -starts_with("Nr"), 
         -batch)
pheno_int_transformed <- pheno_orig %>%
  select(dorsal_area_Res:intvolume_Res, -starts_with("tc_Res"), -starts_with("tg_Res"), 
         -starts_with("volume_Res"), -starts_with("lateral_area_Res"), -starts_with("ldl_Res"),
         -starts_with("dorsal_area_Res"), -starts_with("hdl_Res"), -starts_with("glucose_Res"), 
         starts_with("inttg_Res"))
phenos <- c("intMac_Area","intNeu_Area","intldl_Res", "inthdl_Res", "inttc_Res",          
          "inttg_Res", "intglucose_Res","intlength","intdorsal_area_Res","intlateral_area_Res","intvolume_Res")
genos <- c("apoba","apobb1","apobb2","apoea","ldlra","apob","apoe")

geno_data <- geno[genos]
pheno_data <- pheno_int_transformed[phenos]
all_data <- data1


significants_old_data <- find_significant_associations(phenos, genos, alpha, geno_data, pheno_data, all_data)



#plot(significants_old_data$min_p_CCA_model[[1]])
#plot(significants$min_p_LM_model)
#summary(significants$min_p_LM_model)
#significants$min_p_LM_model$terms
significants_old_data$p_adjusted_CCA[significants_old_data$p_adjusted_CCA == 0] <- NA
significants_old_data$p_adjusted_LM[significants_old_data$p_adjusted_LM == 0] <- NA

significants_old_data <- significants_old_data %>%
  mutate(p_ratio = p_adjusted_LM/p_adjusted_CCA)
hist(log(significants_old_data$p_ratio), breaks=1000)
max_p_ratio <- significants_old_data[which.max(significants_old_data$p_ratio),]
max_p_ratio_cca <- cancor2(geno_data, pheno_data, 
                           xcols = max_p_ratio$genos[[1]], ycols=max_p_ratio$phenos[[1]],
                           subset_selection = FALSE)
plot(max_p_ratio_cca[[1]])
min(significants_old_data$p_adjusted_CCA, na.rm = T)
max_p_ratio
min_p <- significants_old_data[which.min(significants_old_data$p_adjusted_CCA),]
min_p

smpl <- sample.int(nrow(significants_old_data), size=10000)
sod2 <- significants_old_data[smpl,]
sod2$is_sign <- (sod2$p_adjusted_CCA < alpha) + 1
plot(sod2$p_adjusted_CCA, log(sod2$p_ratio), col=sod2$is_sign)
plot(sod2$p_adjusted_LM, log(sod2$p_ratio), col=sod2$is_sign)

#################
X <- geno_data
Y <- pheno_data
xcols <- genos
ycols <- phenos
library(tidyr) 
library(dplyr)
XY <- X %>%
  merge(Y, by="row.names") %>%
  select(ycols, xcols) %>%
  drop_na
cca1 <- cancor(XY[xcols], XY[ycols])
qrx <- qr(XY[xcols])
qry <- qr(XY[ycols])
xcols_new <- colnames(qrx$qr)[1:qrx$rank]
ycols_new <- colnames(qry$qr)[1:qry$rank]
cca2 <- cancor(XY[xcols_new], XY[ycols_new])


nr <- nrow(XY)
inp <- qr.qty(qrx, qr.qy(qry, diag(1, nr, qry$rank)))[1L:qrx$rank, , drop = FALSE]
z <- svd(inp, qrx$rank, qry$rank)

