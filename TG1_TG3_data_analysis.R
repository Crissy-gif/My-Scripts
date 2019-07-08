rm(list=ls(all=TRUE))
setwd("~/Desktop/Master*s project/code")
#data_file contains all the varibles separated by tabs
data_TG1 <- read.csv("TG1_for_canonical_analysis.txt", sep="\t", row.names = "barcode")
###variable selection of phenotype data from data_file
library(dplyr)
pheno_vars <- data_TG1 %>%
  select(intldl_Res_r:MacNeu_Area_cat)
pheno_vars_cat <- data_TG1 %>%
  select(intldl_Res_r:MacNeu_Area_cat, -starts_with("cld_r"),-starts_with( "MacLip_Area_r"),-starts_with("NeuLip_Area_r"),
         -starts_with("MacNeu_Area_r"))
pheno_vars_r <- pheno_vars <- data_TG1 %>%
  select(intldl_Res_r:MacNeu_Area_cat, -ends_with("cat"))

###Variable selection of the genotype data from data_file

geno_vars <-data_TG1 %>%
  select(Dosagemap3k1_160345:Dosagevegfabmm_132431)

#######CCA TG1 cat###################
CCA_TG1_cat <- cancor2(X=geno_vars, Y=pheno_vars_cat)
print(CCA_TG1_cat)
###single to multple_cat
CCAs_TG1_cat <- lapply(names(geno_vars), function(xcol){
  cancor2(X=geno_vars, Y=pheno_vars_cat, xcols=c(xcol))
})
print(CCAs_TG1_cat)
#######CCA TG1 r ###################
CCA_TG1_r <- cancor2(X=geno_vars, Y=pheno_vars_r)
print(CCA_TG1_r)
###single to multple_r
CCAs_TG1_r <- lapply(names(geno_vars), function(xcol){
  cancor2(X=geno_vars, Y=pheno_vars_r, xcols=c(xcol))
})
print(CCAs_TG1_r)

library(CCA)

##############################CCA FUNCTIONS############################################################################
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


###############################TG3 data analysis######################################
data_TG3 <- read.csv("TG3_for_canonical_analysis.txt", sep="\t", row.names = "barcode")
######variable separation
pheno_vars3 <- data_TG3 %>%
  select(intldl_Res_r:MacNeu_Area_cat)
pheno_vars_cat3 <- data_TG3 %>%
  select(intldl_Res_r:MacNeu_Area_cat, -starts_with("cld_r"),-starts_with( "MacLip_Area_r"),-starts_with("NeuLip_Area_r"),
         -starts_with("MacNeu_Area_r"))
pheno_vars_r3 <- data_TG3 %>%
  select(intldl_Res_r:MacNeu_Area_cat, -ends_with("cat"))

###Variable selection of the genotype data from data_file
geno_vars3 <-data_TG3 %>%
  select(Dosagegatad2ab_033103:Dosagelpar2b_062426)

#######CCA Multiple gene to multiple phenotype###################
CCA_TG3_cat <- cancor2(X=geno_vars3, Y=pheno_vars_cat3)
print(CCA_TG3_cat)
#Single gene to multiple pheno
CCAs_TG3_cat3 <- lapply(names(geno_vars3), function(xcol){
  cancor2(X=geno_vars3, Y=pheno_vars_cat3, xcols=c(xcol))
})
print(CCAs_TG_cat3)

#######
CCA_TG3_r3 <- cancor2(X=geno_vars3, Y=pheno_vars_r3)
print(CCA_TG3_r3)
###single to multple_r
CCAs_TG3_r <- lapply(names(geno_vars3), function(xcol){
  cancor2(X=geno_vars3, Y=pheno_vars_r3, xcols=c(xcol))
})
print(CCAs_TG3_r)

###############################Insulin Resistence data analysis######################################
data_IR <- read.csv("IR_POP1_for_canonical_analysis.txt", sep="\t", row.names = "barcode")
######variable separation
pheno_IR<- data_IR %>%
  select(intldl_Res_r:bodyfat_r)

###Variable selection of the genotype data from data_file
geno_IR <-data_IR %>%
  select(D_ins_051222:D_pparg_172834)

#######CCA Multiple gene to multiple phenotype###################
CCA_IR <- cancor2(X=geno_IR, Y=pheno_IR)
print(CCA_IR)
#Single gene to multiple pheno
CCAs_IR <- lapply(names(geno_IR), function(xcol){
  cancor2(X=geno_IR, Y=pheno_IR, xcols=c(xcol))
})
print(CCAs_IR)


##################### Apoe dataset######################################
data_Apoe_file <- read.csv("for_canonical_analysis_clean.txt", sep="\t", row.names = "barcode")

######variable separation
pheno_vars_alls <- data_Apoe_file %>%
  select(intldl_Res_r:MacNeu_Area_cat)
pheno_vars1_cat <- data_Apoe_file %>%
  select(intldl_Res_r:MacNeu_Area_cat, -starts_with("cld_r"),-starts_with( "MacLip_Area_r"),-starts_with("NeuLip_Area_r"),
         -starts_with("MacNeu_Area_r"))
pheno_vars1_r <- data_Apoe_file %>%
  select(ldl:MacNeu_Area, -ends_with("cat"))

###Variable selection of the genotype data from data_file
geno_vars1 <- data_Apoe_file %>%
  select(apobb1:ldlra)

#######CCA Multiple gene to multiple phenotype r###################
CCA_Apoe_r <- cancor2(X=geno_vars1, Y=pheno_vars1_r)
print(CCA_Apoe_r)
#Single gene to multiple pheno
CCAs_Apoe_r <- lapply(names(geno_vars1), function(xcol){
  cancor2(X=geno_vars1, Y=pheno_vars1_r, xcols=c(xcol))
})
print(CCAs_Apoe_r)
