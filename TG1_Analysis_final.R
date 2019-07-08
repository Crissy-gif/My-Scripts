rm(list=ls(all=TRUE))
setwd("~/Desktop/Master*s project/code")
#data_file contains all the varibles separated by tabs
data_TG1 <- read.csv("TG1_for_canonical_analysis_clean.txt", sep="\t", row.names = "barcode")


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

yacca_plot2 <- function (x, ...) 
{
  ncv <- length(x$corr)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(ceiling(sqrt(ncv)-1), ceiling(sqrt(ncv))))
  for (i in 1:ncv) {
    plot(x$canvarx[, i], x$canvary[, i], xlab = "X", ylab = "Y", 
         main = paste("Canonical Variate Plot - Variate", 
                      i, sep = " "), ...)
    abline(mean(x$canvary[, i], na.rm = TRUE) - x$corr[i] * 
             mean(x$canvarx[, i], na.rm = TRUE), x$corr[i])
    text(mean(x$canvarx[, i], na.rm = TRUE), mean(x$canvary[, 
                                                            i], na.rm = TRUE), label = paste("r=", round(x$corr[i], 
                                                                                                         digits = 2), sep = ""), pos = 1, srt = 180/pi * atan(x$corr[i]))
  }
  par(mfrow = c(1, 1), ask = TRUE)
  h <- cbind(x$xvrd, x$yvrd)
  h <- rbind(c(x$xrd, x$yrd), h)
  colnames(h) <- c("X Given Y", "Y Given X")
  rownames(h) <- c("Total", paste("CV", 1:ncv))
  barplot(h, beside = TRUE, main = "Canonical Variate Redundancy Plot", 
          ylim = c(0, max(h)*1.25), ylab = "Fraction of Variance Explained", 
          legend.text = rownames(h), col = rainbow(ncv + 1))
  par(mfrow = c(ceiling(sqrt(ncv)-1), ceiling(sqrt(ncv))))
  for (i in 1:ncv) {
    helio.plot(x, cv = i, main = paste("Structural Correlations for CV", 
                                       i))
  }
  par(mfrow = c(ceiling(sqrt(ncv)-1), ceiling(sqrt(ncv))))
  for (i in 1:ncv) {
    helio.plot(x, cv = i, main = paste("Explained Variance for CV", 
                                       i), type = "variance", axis.circ = c(0.5, 1), range.rad = 25)
  }
}

###################################TG1 Analysis############################################
###variable selection of phenotype data from data_file
library(dplyr)
pheno_vars <- pheno_vars <- data_TG1 %>%
  select(ldl:MacNeu_Area, -ends_with("cat"))

###Variable selection of the genotype data from data_file
geno_vars <-data_TG1 %>%
  select(map3k1:vegfb)
#######CCA TG1 cat###################
CCA_TG1 <- cancor2(X=geno_vars, Y=pheno_vars)
print(CCA_TG1)

yacca::plot.cca(CCA_TG1[[1]])
                

corr_TG1 <- CCA_TG1[[1]]$xcrosscorr
print(corr_TG1)


###single to multple_cat
CCAs_TG1 <- lapply(names(geno_vars), function(xcol){
  cancor2(X=geno_vars, Y=pheno_vars, xcols=c(xcol))
})
print(CCAs_TG1)



###
####################################################
###r^2 =  (SS_tot - SS_yhat)/SS_tot
#########################################################################################
#Subsetting the parameters from the CCA results
#Taking all the xcoef
Y_coef <- data.frame(lapply(CCAs_TG1, function(x) x[[1]]$ycoef))
names(Y_coef) <- names(geno_vars)
#Taking all corrs
all_corr <- data.frame(lapply(CCAs_TG1, function(x) x[[1]]$corr))
names(all_corr) <- names(geno_vars)
#Taking all the xcrosscorrs
Y_ycrosscorrs <- data.frame(lapply(CCAs_TG1, function(x) x[[1]]$ystructcorr))
names(Y_ycrosscorrs) <- names(geno_vars)
#Taking all the var
Y_varex <- data.frame(lapply(CCAs_TG1, function(x) x[[1]]$ystructcorrsq))
names(Y_varex) <- names(geno_vars)
Y_crosscorr <- data.frame(lapply(CCAs_TG1, function(x) x[[1]]$ycrosscorr))
names(Y_crosscorr) <- names(geno_vars)
write.csv(Y_coef, "Y_crosscorr.csv")
P_vals <- data.frame(lapply(CCAs_TG1, function(x) x[[2]]$p.value))
names(P_vals) <- names(geno_vars)
#######################################################################################
library(grid)     ## Need to attach (and not just load) grid package
library(pheatmap)


## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

#install.packages("pheatmap")
library(pheatmap)
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

myheat <- function(data){
  # Center the color scale around 0
  test <- data
  paletteLength <- 100
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
  
  ## Try it out
  library(RColorBrewer)
  pheatmap(data, color = myColor, breaks = myBreaks, display_numbers = T, fontsize = 20)
}
myheat(Y_coef)
myheat(Y_crosscorr)
myheat(Y_varex)
myheat(Y_ycrosscorrs)
myheat(all_corr)
myheat(corr_mvsm)
myheat(Y_varex)
library(statmod)

