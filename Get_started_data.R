#Author ***Crisencia Sizyoongo***

rm(list=ls(all=TRUE))
setwd("~/Desktop/Master*s project/code")
#data1 contains all the varibles separated by tabs
data1 <- read.csv("LDLR_APOE_APOBb1_for_analysis.txt", sep="\t", row.names = "barcode")


str(data1)

library(dplyr)

pheno <- data1 %>%
  select(matches("apo|^ldl$"))

data3 <- pheno[,sapply(X = pheno, FUN = function(x) sum(is.na(x))) < 10]

lm2 <- lm(data=data3, formula = length ~ .)
summary(lm2)

lm3 <- lm(data=data3, formula = length ~ apoea)
summary(lm3)
plot(data3$apoea, data3$length)
abline(lm3)


data3_batch5 <- data3[data1$batch==5,]
lm4 <- lm(data=data3_batch5, formula = length ~ apoea)
summary(lm4)
plot(lm4)
plot(data3_batch5$apoea, data3_batch5$length)
abline(lm4)

lm5 <- lm(data=data3, formula = ldl ~ .)
summary(lm5)
# R squared adjusted: 0.0154
plot(data3$apoea, data3$ldl)
abline(lm5)

data_rf <- pheno[!is.na(pheno$ldl),]
data_train_idx = sample(1:nrow(data_rf), size=nrow(data_rf)*0.8, replace = FALSE)
data_test_idx = (1:nrow(data_rf))[-data_train_idx]
library(party)
cf1 <- cforest(formula = ldl ~ ., data = data_rf[data_train_idx,])
cf1_prediction <- predict(cf1, data_rf[data_test_idx,], OOB = TRUE)
cor(data_rf[data_test_idx, 'ldl'], cf1_prediction)
# [1,] 0.1477008   <- weak positive correlation
# squared: 0.02181553
cor(data_rf[data_train_idx, 'ldl'], predict(cf1))
# [1,] 0.4113123   <- overfitting on the training dataset

pt <- prettytree(cf1@ensemble[[1]], names(cf1@data@get("input"))) 
nt <- new("BinaryTree") 
nt@tree <- pt 
nt@data <- cf1@data 
nt@responses <- cf1@responses 

plot(nt, type="simple")


data_lm <- data_rf[,sapply(X = data_rf, FUN = function(x) sum(is.na(x))) < 10]
lm6 <- lm(data=data_lm[data_train_idx,], formula = ldl ~ .)
summary(lm6)
# R squared adjusted: 0.02043
lm6_prediction <- predict(lm6, data_lm[data_test_idx,])
cor(data_lm[data_train_idx,'ldl'], predict(lm6))
# [1] 0.1910108
cor(data_lm[data_test_idx,'ldl'], lm6_prediction)
# 1] 0.02870602
# squared: [1] 0.0008240357



ct1 <- ctree(data= data_rf, formula = ldl ~ ., 
             controls = ctree_control(mincriterion = 0.1, minbucket = 0, maxdepth = 0))
plot(ct1)



######## CCA #########
#%>% match finder
#X <- data1 %>%
  select(starts_with("apo"))
#X2 <- X[,sapply(X = X, FUN = function(x) sum(is.na(x))) < 1]

#Y <- data1 %>%
  select(matches("hdl|ldl|mac"))
#Y2 <- Y[,sapply(X = Y, FUN = function(x) sum(is.na(x))) < 1]

#cca1 <- cancor(X2, Y2)
#####################################################################
#CCA Part 1                                                         #
#####################################################################

#Step 1. Varible sorting and distribution check#

# colRange <- colnames(data1) %>% {match("cld", .):match("intvolume_Res", .)}
#Phenotype variable selesction
#pheno <- data1 %>%
#  select(cld:intvolume_Res, -ends_with("_Res"), -line,
#         -sequencing, -starts_with("Annotation"), -mutant)
#head(pheno)
#unique(data1$sequencing)
##library(magrittr)
#pheno %$% unique(batch)
#hist(pheno$cld, breaks = 30)
#cld not bell shaped
#hist(log(pheno$cld), breaks = 30)
#Transform cld before use
#install.packages("e1071")
#install.packages("moments")
#library(moments)
#log_cdl <- log(pheno$cld) %>%
#  .[is.finite(.)]
#skewness(log_cdl)
#kurtosis(log_cdl)
#qqnorm(log_cdl)
#qqline(log_cdl, col = 2)
#Variable transformation by log and squared
#par(mfrow=c(3,3))
#for(i in 1:ncol(pheno)){
#  hist((pheno[,i]), main = paste0(names(pheno)[i]), breaks = 30)
#  hist(log(pheno[,i]), main = paste0("log_", names(pheno)[i]), breaks = 30)
#  hist(sqrt(pheno[,i]), main = paste0("sqrt_", names(pheno)[i]), breaks = 30)
#}
  
#Select transformed varaibles
#sapply(pheno, hist)
#sapply(pheno, class)
#pheno$mutant
#Genotype variables
#geno <-data1 %>%
#  select(apobb1:apoe_ldlr,-miss)
#names(geno)
#Check the distribution of the variables
#for(i in 1:ncol(geno)) {
 # hist((geno[,i]), main = paste0(names(geno)[i]), breaks = 30)
#  hist(log(geno[,i]), main = paste0("log_", names(geno)[i]), breaks = 30)
#}
#Select the variables that has been checked for distribution
#All the variables that starts with int..follow a normal distribution with raw data.
#The int.variables have been tranformed using the inverse-normal-transformation
#Batch varible is just numbers put not from the observations

#transformations <- rep("", ncol(pheno))
#names(transformations) <- names(pheno)
#transformations[c("bld", "cld", "Mac_Neu_Area", "Mac_Lip_Area", "NeuLip_N", 
#                  "NeuLip_Area", "Mac_NeuLip_Area",
#                  "Volume", "tc")] <- "log"
#transformations[c("Mac_Area", "Neu_N", "Neu_Area", "bld", "glucose", "ldl")] <- "sqrt"
#pheno_trans <- data.frame(mapply(FUN = function(trans, pheno_col){
 # if(trans == ""){
 #   pheno_col
 # } else if(trans == "sqrt"){
 #   sqrt(pheno_col)
 # } else {
#    log(pheno_col)
#  }
#}, transformations, pheno))

#names(pheno_trans) <- paste0(transformations, "_", names(pheno))
#par(mfrow=c(3,3))
#mapply(function(x, name) hist(x, main=name), pheno_trans, names(pheno_trans))

#Step 2. Pairwise testing Gene vs dependant varibles

#########################################################################################################
#                                    Steps in the Analysis                                              #                               
#########################################################################################################
#Step 1#Variable sorting #
##########################

#pheno_orig contains untouched raw phenotype variables
pheno_orig <- data1 %>%
  select(cld:intvolume_Res, ends_with("_Res"), -line, -sequencing, 
         -starts_with("Annotation"), -mutant, starts_with("int"), -starts_with("Nr"), 
         -batch)
#pheno_derived contains inverse-normal-tansformed phenotype variables
pheno_derived <- data1 %>%
  select(tc_Res:intvolume_Res, -line, -intvolume_Res, -intlateral_area_Res, -lateral_area_Res
         , -intvolume_Res, - volume_Res, -intdorsal_area_Res, -intlength, -intlateral_area_Res, -intvolume_Res, -lateral_area_Res
       -sequencing, -starts_with("Annotation"), -mutant, -starts_with("Nr"))
#pheno_covariants contains body size (length, dorsal, lateral, volume & batch)
pheno_covariants <- data1 %>%
  select(starts_with("lateral_area") ,starts_with("length"), starts_with("dorsal_area"),
         starts_with("volume"), -ends_with("_Res"))
#pheno_not_transformed contains the raw data of the variables 
pheno_not_tranformed <- pheno_orig %>%
   select(Mac_Area:hdl, -ends_with("_N"), -starts_with("tod"), -starts_with("protein"), 
          starts_with("glucose"), -ends_with("_Res"))
pheno_int_transformed <- pheno_orig %>%
  select(dorsal_area_Res:intvolume_Res, -starts_with("tc_Res"), -starts_with("tg_Res"), 
         -starts_with("volume_Res"), -starts_with("lateral_area_Res"), -starts_with("ldl_Res"),
         -starts_with("dorsal_area_Res"), -starts_with("hdl_Res"), -starts_with("glucose_Res"), 
         starts_with("inttg_Res"))
pheno_proteins <- pheno_not_tranformed %>%
  select()

#*************************************************************************************************
#genotype variables
#************************************************************************************************
#pheno_final contains varibles selected to be used in the final analysis
#geno contains the genes selected excluding the genes that end with "homo"
geno <-data1 %>%
 select(apoba:apoe_ldlr, -ends_with("_homo"), -miss, -apoe_ldlr, - apob_ldlr, -apoe_apob)
geno2 <- geno%>%
   select(ldlr:apoe)
geno3 <- geno%>%
  select(apoba:ldlrb)

#Step 2#
# Merge the independent and dependent variables. We do this in order to drop the same
# observations from both data sets
#phenogeno <- merge(geno,pheno_derived, by = "row.names")
#row.names(phenogeno) <- phenogeno$Row.names
#phenogeno <- phenogeno %>%
#  select(-matches("Row.names"))

#**********************************************************************************************************
#Step 2# Cancor 2 function
#***********************************************************************************************************
#CCA can only be performed with 2 independent datasets
#Making the barcode to be the row.names
#install.packages("yacca")
library(yacca)
cancor2 <- function(X, Y, xcols=NA, ycols=NA){
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
  #print(summary(XY))
  # TODO: how are variables selected when the matrix is singular? (i.e. when variables are too correlated)
  cca1 <- cancor(XY[xcols], XY[ycols])
  selected_xcols <- row.names(cca1$xcoef)
  selected_ycols <- row.names(cca1$ycoef)
  #print((cca1))
  #cca2 <- yacca::cca(XY[selected_xcols], XY[selected_ycols])
  cca2 <- yaccacca2(XY[selected_xcols], XY[selected_ycols])
  #print(cca2)
  fcca2 <- yacca::F.test.cca(cca2)
  list(cca2, fcca2)
  #sink(cca1,"CCA1.txt")
}
library(car)
#****************************************************************************************************************************


##############################################################################################
#Step 3 CCA single pheno vs multiple geno#
library(CCA)
#Single trait vs multiple genes loop:
ccas <- lapply(names(pheno_int_transformed), function(ycol){
  cancor2(X=geno3, Y=pheno_int_transformed, ycols=c(ycol))
})

#sink("RESULTz.txt")
ccas
sink()
#Subsetting the parameters from the CCA results
#Taking all the xcoef
X_coef <- data.frame(lapply(ccas, function(x) x[[1]]$xcoef))
names(X_coef) <- names(pheno_int_transformed)
#Taking all corrs
all_corr <- data.frame(lapply(ccas, function(x) x[[1]]$corr))
names(all_corr) <- names(pheno_int_transformed)
#Taking all the xcrosscorrs
X_xcrosscorrs <- data.frame(lapply(ccas, function(x) x[[1]]$xstructcorr))
names(X_xcrosscorrs) <- names(pheno_int_transformed)
#Taking all the var
X_varex <- data.frame(lapply(ccas, function(x) x[[1]]$xstructcorrsq))
names(X_varex) <- names(pheno_int_transformed)
X_crosscorr <- data.frame(lapply(ccas, function(x) x[[1]]$xcrosscorr))
names(X_crosscorr) <- names(pheno_int_transformed)
write.csv(X_coef, "X_crosscorr.csv")
P_vals <- data.frame(lapply(ccas, function(x) x[[2]]$p.value))
names(P_vals) <- names(pheno_int_transformed)
write.csv(P_vals, "P.vals.csv")

p_vals <- sapply(ccas, function(x) x[[2]]$p.value)
corrs <-sapply(ccas, function(x) x[[2]]$corr)
library(xlsx)
write.xlsx(corrs,"corrs.xlsx")


sign_level <- 0.05
#results <- names(pheno_derived)[p_vals <= sign_level]
#results


results <- data.frame(pheno = names(pheno_derived)[p_vals <= sign_level],
           p_value = p_vals[p_vals <= sign_level], 
           corr = corrs[p_vals <= sign_level])
write.table(results, file = "RESULTS1.txt", row.names = F, quote = F, sep = "\t")
#cca1 <- cancor(phenogeno2$cld, phenogeno2[c("apoba", "apobb1", "apobb2", "apoea", "ldlr")])

#yacca::plot.cca(ccas[[c(1,2)]])
yacca::plot.cca(ccas[[c(2,1)]])

par(mfrow=c(2,3))
for(i in 1:length(ccas)){
  yacca::helio.plot(ccas[[c(i,1)]], cv = 1, 
          main = paste("Structural Correlations for CV", 1))
}

par(mfrow=c(2,3))
for(i in 1:length(ccas)){
  yacca::helio.plot(ccas[[c(i,1)]], cv = 1, 
                    main = paste("Canonical Variate Plot-Variate 1", 1))
}

par(mfrow=c(2,3))
for(i in 1:length(ccas)){
  yacca::helio.plot(ccas[[c(i,1)]], cv = 1, 
                    main = paste("Explained Variance for CV 1", 1))
}


install.packages("gplots")
library(gplots)

heatmap.2(as.matrix(X_crosscorr), col=bluered(10), srtCol=25)


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

install.packages("pheatmap")
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
  pheatmap(data, color = myColor, breaks = myBreaks, display_numbers = T)
}
myheat(X_coef)
myheat(X_crosscorr)
myheat(X_varex)
myheat(X_xcrosscorrs)
myheat(all_corr)
myheat(corr_mvsm)
library(HiveR)

###########################################################################################
#Meeting with Ci
#Extract the xcoef and in one file: Put then 
#Extract all the p_value  and put them in a single file
#

#####################################################################################
#Step 4 CCA multiple pheno vs multiple geno#
cca1 <- cancor2(X=geno, Y=pheno_derived, ycols = c("tc_Res", "inttc_Res", "intldl_Res", "ldl_Res","tg_Res", "glucose_Res" ,
                                                   "inttg_Res",  "inttg_Res", "intNeu_Area", "intglucose_Res", "hdl_Res","tg_Res" , "intMac_Area"  ))
#######################################################################################################################################################
#/ genes vs non tansformed variables
pheno_nt_t <- sapply(pheno_not_tranformed,log)
cca6 <- cancor2(X=geno3, Y=pheno_not_tranformed, ycols = c("tc", "ldl", "glucose" ,"tg", "hdl" , "Mac_Area", "MacNeuLip_Area","Neu_Area","MacNeu_Area",
                                                           "NeuLip_Area","MacLip_Area" ))
cca6_2 <- cancor2(X=geno3, Y=pheno_nt_t, ycols = c("tc", "ldl", "glucose" ,"tg", "hdl" , "Mac_Area", "MacNeuLip_Area","Neu_Area","MacNeu_Area",
                                                           "NeuLip_Area","MacLip_Area" )) # the zero values gives inf..so cant compute

cca_not_bodysize <- cancor2(X=geno3, Y=pheno_not_tranformed, ycols = c("length", "volume","lateral_area","dorsal_area")) #Not signifcant
cca_not_cell_area <- cancor2(X=geno3, Y=pheno_not_tranformed, ycols = c("Mac_Area", "MacNeuLip_Area","Neu_Area","MacNeu_Area",
                                                                       "NeuLip_Area","MacLip_Area")) #significant with p = 0.05024

yacca::plot.cca(cca_not_cell_area[[1]])                                                                    
yacca::plot.cca(cca6[[1]])
yacca::plot.cca(ccas[[c(1,2)]])

#cca6 for multple vs multiple gives the f test of p values of 0.004978 in CV1
cca_nt_transformed <- cancor2(geno3, pheno_not_tranformed, ycols = c("tc", "tg","ldl","hdl","glucose"))## cca_nt_transformed gives the p = 0.004073


#cca1 <- cancor2(X=geno, Y=pheno_derived, ycols = c("tc_Res", "inttc_Res", "intldl_Res", "ldl_Res","tg_Res", "glucose_Res" ,
#                                                   "inttg_Res",  "inttg_Res", "intNeu_Area", "intglucose_Res", "hdl_Res","tg_Res" , "intMac_Area"  ))

cca_mvsm <- cancor2(X=geno3, Y=pheno_int_transformed, ycols = c( "intldl_Res" , "inthdl_Res",  "inttg_Res","intglucose_Res","intNeu_Area",
                                                                "intMac_Area", "intlateral_area_Res", "intvolume_Res",
                                                                "intlength", "intdorsal_area_Res")) #cca_int give p = 0.01258 *
corr_mvsm <- cca_mvsm[[1]]$xcrosscorr
myheat()

#cca_int 7 genes vs 11 traits show a p value of 0.01096 *
yacca::plot.cca(cca_mvsm[[1]])                       



p_val <- cca1[[2]]$p.value[1]
p_val # Unexpectedly high. Higher the more variables we include as phenos.
yacca::plot.cca(cca1[[1]])


yacca::plot.cca(cca1[[2]])


###########################################################################################
#New steps/2018.03.8
#Pheno derived should contain the int and the res variables without the body size factors e.g length, volume, dorsal, lateral etc
#Steps
#Select a single trait from the pheno derived and merge it to the geno by barcode. Both the dataset to be merged must contain barcode.
#drop the na to the dataset created
#Do the same for the multiple traits (pheno_derived) and the mulitple geno as above
# Drop observations with missing values
#########################################################################################
library(tidyr)
#phenogeno2 <- phenogeno %>%
#  drop_na

# Run CANCOR
#Single trait vs multple genes
#cca1 <- cancor(phenogeno2$cld, phenogeno2[c("apoba", "apobb1", "apobb2", "apoea", "ldlr")])
#cca1

# Run multiple linear regression
#lm1 <- lm(data = phenogeno2, formula = reformulate(c("apoba", "apobb1", "apobb2", "apoea", "ldlr"), response="cld"), )
#summary(lm1)

# The following shows how CCA with one dependent variable is essentially 
# achieving the same thing as multiple linear regression (in a different way)
#cca11s <- cca1$ycoef[,1]/cca1$xcoef
#cca11sn <- scale(cca11s)

#cca11sn
#scale(lm1$coefficients[2:length(lm1$coefficients)])
###########################################################################################
#############################################################
library(xlsx)
trial <- matcor(X=geno, Y=pheno_derived)
#lapply(trial, function(trial) write.csv( data.frame(trial), 'test4.csv'  , append= T, sep= " " ))
trial2 <- matcor(geno3, pheno_not_tranformed)
trail3 <-as.data.frame(trial2, row.names = NULL)
img.matcor(trial2,type=2)





############################################################################
library(gmodel)
library(swirl)
########################################################################
cca4 <- cancor2(geno2, pheno_not_tranformed, ycols = c("tc", "tg","ldl","hdl","glucose"))
cca4_1 <- cancor2(geno3, pheno_not_tranformed, ycols = c("tc", "tg","ldl","hdl","glucose"))

cca5 <- cancor2(geno3, pheno_not_tranformed, ycols = c("dorsal_area", "lateral_area", "volume"))
yacca::plot.cca(cca4_1[[1]])


yacca::plot.cca(cca4_1[[2]])


cca5_1 <- cancor2(geno2, pheno_not_tranformed, ycols = c("dorsal_area", "lateral_area", "volume")) 

#############################################################################################################################
#Multiple linear regression
tg_geno <-lm(formula =tg ~ apoba + apobb1 + apobb2 + apoea + apoeb + ldlra + ldlrb , data = data1) #P = 0.49
tc_geno <- lm(formula =tc ~ apoba + apobb1 + apobb2 + apoea + apoeb + ldlra + ldlrb , data = data1) #p = 0.06
hdl_geno <- lm(formula =hdl ~ apoba + apobb1 + apobb2 + apoea + apoeb + ldlra + ldlrb , data = data1) #p = 0.09
glucose_geno <- lm(formula =glucose ~ apoba + apobb1 + apobb2 + apoea + apoeb + ldlra + ldlrb , data = data1) #P =0.3212
ldl_geno <- lm(formula =ldl ~ apoba + apobb1 + apobb2 + apoea + apoeb + ldlra + ldlrb , data = data1) #  p = 0.064
######################################################################################
# Custom yacca functions
######################################################################################

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
########################################################################
#o keep it consistent with the analyses performed earlier, 
#it would be great if you could focus on [cld, MacLip_Area, NeuLip_Area] 
#and/or [tg, ldl, hdl, glucose] and [length, dorsal_area]

#Non transformed
final_var1_nt <- pheno_orig%>% 
   select (starts_with("tg"), starts_with("ldl"),starts_with ("hdl"),starts_with("glucose")) 
#"tg, ldl, hdl, glucose"
final_var2_nt <-pheno_orig%>%  #"cld, MacLip_Area, NeuLip_Area"
 select(starts_with ("cld"), starts_with ("MacLip_Area"), starts_with ("NeuLip_Area"))
final_var3_nt <- pheno%>%  #"length, dorsal_area",)
 select(starts_with ("length"),starts_with ("dorsal_area"))
#Inversly transformed
 final_var1_int <- pheno_orig%>% 
   select (starts_with("inttg_Res"), starts_with("intldl_Res"),starts_with ("inthdl_Res"),starts_with("intglucose_Res")) 
 #"tg, ldl, hdl, glucose"
 final_var2_int <-pheno_orig%>%  #"cld, MacLip_Area, NeuLip_Area"
   select(starts_with ("cld"), starts_with ("intMac_Area"), starts_with ("intNeu_Area"))
 final_var3_int <- pheno_orig%>%  #"length, dorsal_area",)
 select(starts_with ("intlength"),starts_with ("intdorsal_area_Res"))
 
#Run the cancor2 with the inversely tranformed varibles multiple geno vs multiple pheno!
CCA_final_int_lipids <- cancor2(X=geno3, Y=final_var1_int, ycols = c("inttg_Res","intldl_Res", "inthdl_Res","intglucose_Res" ))
#CCA-final_int_lipids has a p.value of 0.006622 **
ccas_lipds_svm <- lapply(names(final_var1_nt), function(ycol){
  cancor2(X=geno3, Y=final_var1_int, ycols=c("intldl_Res"),row.names("intldl_Res")) #,"intldl_Res", "intglucose_Res","inthdl_Res"))
})

CCA_final_int_lipids_3genes <- cancor2(X=geno2, Y=final_var1_int, ycols = c("inttg_Res","intldl_Res", "inthdl_Res","intglucose_Res" ))
#3 genes + the lipids!
CCA_final_int_cellareas <- cancor2(X=geno3, Y=final_var2_int, ycols = c("cld","intMac_Area", "intNeu_Area")) #p.value is 0.006514 **

CCA_final_int_bodysize <- cancor2(X=geno3, Y=final_var3_int, ycols = c("intlength","intdorsal_area_Res" )) #No association

######################################################################
#Evaluation of the 2 methods!
######################################################################
phenos <- c("hdl", "ldl", "glucose","tg")
genos <- c("apoba","apobb1","apobb2","apoea","apoeb","ldlra","ldlrb")
alpha <- 0.025
#alpha <- 0.0005


library(ggm)
pset_geno <- powerset(set = phenos)
pset_geno2 <- pset_geno[sapply(pset_geno, function(x) {length(x) > 1})]
pset_pheno <- powerset(set = genos)
pset_pheno2 <- pset_pheno[sapply(pset_pheno, function(x) {length(x) > 1})]
sign_LM <- vector()
sign_CCA <- vector()
n_iter <- 0
for(geno_set in pset_geno2){
  for(pheno_set in pset_pheno2){
    n_iter <- n_iter + 1
    # CCA
    res <- cancor2(geno3, pheno_not_tranformed, xcols = geno_set, ycols=pheno_set)
    f_res <- res[[2]]
    p_value <- f_res$p.value[1]
    is_significant <- p_value < alpha
    sign_CCA <- c(sign_CCA, is_significant)
    #if(is_significant){print("Genos: ");print(geno_set);print("Phenos: ");print(pheno_set);print("-------------")}

    # LM
    any_significant <- FALSE
    for(pheno in pheno_set){
      b <- paste(geno_set, sep ="+")
      lm_formula <- as.formula(paste0(pheno, " ~ ", b))
      tg_geno <- lm(formula = lm_formula, data = data1)
      sm1 <- summary(tg_geno)
      p_value <- pf(sm1$fstatistic[1],sm1$fstatistic[2],sm1$fstatistic[3],lower.tail=FALSE)
      if(p_value < alpha) {any_significant <- TRUE}
    }
    sign_LM <- c(sign_LM, any_significant)
  }
}
print(paste0("Number of iterations: ", n_iter))
total_LM <- sum(sign_LM)
total_CCA <- sum(sign_CCA)
print(paste0("Total significant sets for LM: ", total_LM))
print(paste0("Total significant sets for CCA: ", total_CCA))
sum(sign_LM[sign_CCA]) # Both are significant
sum(sign_CCA + sign_LM == 1) # Only one is significant
sum(sign_CCA + sign_LM == 0) # None is significant
sum(sign_CCA & !sign_LM) # CCA is significant but LM is not
sum((!sign_CCA) & sign_LM) # LM is significant but CCA is not

