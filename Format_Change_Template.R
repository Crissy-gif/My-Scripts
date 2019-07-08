
library(tidyverse)
library(dplyr)
library(tidyr)


data <-read.csv("MC_Servier_1809_X01_peptideExport.csv", sep = ";")
data2 <- read.csv("MC_Servier_1809_X01_org.csv", sep = ";")


#Select variables to transform 
sel <- data %>%
select(Master,TM,Sample)
# converting commas to dots
#sel$TM2 <- sel$TM %>%
#str_replace_all(",", ".") 

sel_nest <- sel %>%
 group_by(Master,Sample) %>%
  summarise_all(funs(trimws(paste(., collapse = ","))))
   #summarise(TM = paste(TM, collapse = ","))
########################

wide_format <- spread(sel_nest, Sample, TM)

wide_format$peptide.TMs.Control <- paste(wide_format$ctrl1, wide_format$ctrl2)
wide_format$peptide.TMs.Drug <- paste(wide_format$drug1, wide_format$drug2)

wide_format <- wide_format %>%
 select(-starts_with("V1" ), -starts_with( "ctrl1"),-starts_with("ctrl2"),
         -starts_with( "drug1"),-starts_with( "drug2"))

wide_merged <- wide_format %>%
 merge(data2, by.x = "Master", by.y="Acc.Number", all.y = TRUE)

#wide_merged <- wide_merged%>%
 # select(-starts_with( "Ctrl1.A"),-starts_with("Ctrl1.B"),-starts_with("Drug1.A"),-starts_with("Ctrl2.B"), -starts_with("Ctrl2.A"),
   #      -starts_with("Drug2.B"), -starts_with("Ctrl1.scal"), -starts_with("Ctrl1.scal"),-starts_with("Drug1.scal"))
wide_merged$peptide.TMs..ctrl <- NULL
wide_merged$peptide.TMs..drug <- NULL
wide_merged2 <- wide_merged[, c(1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,2,3)]
write.csv(wide_merged2, "wide_merged.csv", row.names = FALSE)
#---------------------------------------------------------------------------#

  

 






