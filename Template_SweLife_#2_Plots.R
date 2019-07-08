########################################################################################
# R scripts for generating volcanol and correlation plots in the meantime              #
########################################################################################

rm(list = ls())
setwd("~/Documents/Data Analysis/SweLife MCR Analysis/Swelife_PB0074/PBMC#2_Plots")
#++++Load the libraries
library(ggplot2)
library(plotly)
library(gridExtra)
library(magrittr)
library(dplyr)
library(tidyverse)
library(GGally)
library(shiny)

#++++load data
file_ <- read.csv(file = "MCR_SweLife_1810_PBMC_3_sheet.csv")

#++++Coverting varibles with class factors to numeric
file_ %<>%
  mutate_at(vars(ends_with("_MTX")), funs(as.numeric(as.character(.)))) %>%
  mutate_at(vars(ends_with("_MP6")), funs(as.numeric(as.character(.)))) %>%
  mutate_at(vars(ends_with("_CytA")), funs(as.numeric(as.character(.)))) 
# 
# file_2 <- file_
# file_2["Gene" == "SERPINE1", ] <- NULL #Try to remove this gene

#++++Volcano plot for each compd. 
MTX_var <- file_ %>%
  select(Gene,Amp_MTX,pVals_MTX) #Replace the compd name 

MTX_df <- data.frame(Amp =MTX_var$Amp_MTX,
                     Pval =MTX_var$pVals_MTX) #Run each compd one at a time

MTX_f <- cbind(MTX_var,MTX_df)

MTX_f <-MTX_f %>%
  select(-pVals_MTX,-Amp_MTX)
#++++Add a grouping column
MTX_f["group"] <- "Not.Significant"

MTX_f[which(MTX_f[3] < 0.05 &MTX_f[2] <= -0.15), "group"] <- "Dest.significant"
MTX_f[which(MTX_f[3] < 0.05 &MTX_f[2] >= 0.15), "group"] <- "Stab.significant"

#++++make the Plot.ly plot
p <- plot_ly(data =MTX_f, x = ~Amp, y = ~-log10(Pval), text = ~Gene, mode = "markers", color = ~group) %>%  #Choose colour of your choose colors = c("red","green","blue")
layout(title ="PBMC#3_MTX Plot") #Change the PBMC #accordly
p
#++++Save as a html interactive plot
htmlwidgets::saveWidget(as.widget(p),"PBMC#3_MTX Plot.html")#Change the PBMC #accordly
#++++filter out the significant proteins with the amplitude of 0.15 and P value of 0.05
Prot_filtered <- filter(MTX_f, Pval < 0.05, Amp <= -0.15 | Amp >= 0.15)
write.csv(Prot_filtered, "PBMC#3_MTX_Protein_filtered.csv") #Change the PBMC# accordly
#++++Pairwise correlations
corr_data <- file_ %>%
  select(Gene,Amp_VCR, Amp_MTX,Amp_MP6,Amp_CytA)
corr_plot <- corr_data %>%
  rename(VCR = Amp_VCR, MTX = Amp_MTX, MP6 = Amp_MP6, CytA = Amp_CytA)

corr <- ggpairs(corr_plot[,2:5], title = "PBMC#3 Chart Correlation ", mapping=ggplot2::aes(text = corr_plot[,1])) #Works well :)
corr_It <- ggplotly(corr)
corr_It
htmlwidgets::saveWidget(as.widget(corr_It), "PBMC#3 Chart Correlation.html")

#***********************************************************************************************************
#************************For tableau file preparartion for Volcano plot*************************************
#setwd("~/Documents/Data Analysis/SweLife MCR Analysis/PB0075/PB0075 SweLife")
library (ggplot2)

file_ <- read.csv(file = "MCR_SweLife_1811_PBMC_4_sheet.csv")

#MTX_df <- data.frame(log10Amp = file_$Amp_MTX,
#                   logPval = -log10(file_$pVals_MTX)) # For comp MTX

MTX_df <- data.frame(log10Amp = file_$Amp_MTX,
                     logPval = -log10(file_$pVals_MTX)) # For comp MTX 

#MTX_df <- data.frame(log10Amp = file_$Amp_MTX,
#                    logPval = -log10(file_$pVals_MTX)) # For comp MTX

#MTX_df <- data.frame(log10Amp = file_$Amp_MTX,
#                   logPval = -log10(file_$pVals_MTX)) # For comp MTX 
#################################################################################################################

ggplot(MTX_df, aes(x=log10Amp, y=logPval)) +
  geom_point() # Test the volcano plot

MTX_dfs <- cbind(file_ , MTX_df)

write.csv(MTX_dfs, "PB0072_MTX.csv")
#------------------filter the proteins------------------------------------------------------------------------
#Filter all the proteins above logPV 2 and log10FC below -0.5 and above  +0.5

Val_Sel <-MTX_dfs %>%
  select(Acc.Number, Gene, Peptides, log10Amp, logPval)

#prot_filtered2 <- prot_filtered %>%
Prot_filtered <- filter(Val_Sel, logPval > -log10(0.05), log10Amp <= -0.15 | log10Amp >= 0.15)

#Prot_filtered <- filter(Val_Sel, logPval > -log10(0.01), log10Amp < -0.15 | log10Amp > 0.15)
write.csv(Prot_filtered, "PB0072_MTX_Protein_filtered.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Correlation
corr_data <- file_ %>%
  select(Gene,Amp_MP6, Amp_MTX,Amp_MP6,Amp_CytA)

corr_plot <- corr_data %>%
  rename(MP6 = Amp_MP6, MTX = Amp_MTX, MP6 = Amp_MP6, CytA =  Amp_CytA)

corr <- ggpairs(corr_plot[,2:5], title = "PBMC#2 Chart Correlation ", mapping=ggplot2::aes(text = corr_plot[,1])) #Works per
corr <- ggpairs(corr_plot[,2:5], title = "PB0072- Chart correlation ", mapping=ggplot2::aes(text = corr_plot[,1])) #Works per

#corr <- plot_ly::ggpairss(corr_plot, title = "PB0072- Chart correlation ")

corr_It <- ggplotly(corr)
corr_It
htmlwidgets::saveWidget(as.widget(corr_It), "PBMC#2 Chart Correlation.html")

########################################################################
#correlation plot with performance Analytics
# library("PerformanceAnalytics")
# chart.correlation (corr_data, histogram = TRUE, pch=19) works well 
#chart.correlation (corr_plot, histogram = TRUE, pch=25, main = "PB0072")
##########################################################################

# Correlation with Shiny!!!Working in progress! R server is need to share the results
corr_data <- file_ %>%
  select(Gene,Amp_MP6, Amp_MTX,Amp_MP6,Amp_CytA)

corr_plot <- corr_data %>%
  rename(MP6 = Amp_MP6, MTX = Amp_MTX, MP6 = Amp_MP6, CytA =  Amp_CytA)

cor_mat1 <- as.matrix(corr_plot)
cor_mat <- cor(corr_plot[, 2:5], use = "complete.obs") 
cor_hmisc <- rcorr(cor_mat1)

#++++compute a correlation  matrix
corr_plot$Gene <- NULL
correlation <- round(cor(corr_plot, use = "complete.obs"),3)
#correlation <- round(cor(mtcars), 3)

nms <- names(corr_plot)
#nms <- names(mtcars)

ui <- fluidPage(
  mainPanel(
    plotlyOutput("heat"),
    plotlyOutput("scatterplot")
  ),
  verbatimTextOutput("selection")
)

server <- function(input, output, session) {
  output$heat <- renderPlotly({
    plot_ly(x = nms, y = nms, z = correlation , 
            key = correlation , type = "heatmap", source = "heatplot") %>%
      layout(xaxis = list(title = ""), 
             yaxis = list(title = ""))
  })
  
  output$selection <- renderPrint({
    s <- event_data("plotly_click")
    if (length(s) == 0) {
      "Click on a cell in the heatmap to display a scatterplot"
    } else {
      cat("You selected: \n\n")
      as.list(s)
    }
  })
  
  output$scatterplot <- renderPlotly({
    s <- event_data("plotly_click", source = "heatplot")
    if (length(s)) {
      vars <- c(s[["x"]], s[["y"]])
      print(vars)
  #    d <- setNames(mtcars[vars], c("x", "y"))
      d <- setNames(corr_plot[vars], c("x", "y"))
      d <- d[complete.cases(d),]
      yhat <- fitted(lm(y ~ x, data = d))
      plot_ly(d, x = ~x) %>%
        add_markers(y = ~y) %>%
        add_lines(y = ~yhat) %>%
        layout(xaxis = list(title = s[["x"]]), 
               yaxis = list(title = s[["y"]]), 
               showlegend = FALSE)
    } else {
      plotly_empty()
    }
  })
  
}

runApp(host="0.0.0.0",port=5050)
shinyApp(ui, server)

#***************************************************************************************
#Boxplots for the compds
library(stringr)
file_format <- proteinData %>%
  as.data.frame %>%
  select(1:5) %>%
  gather(key=compound, value=amps, -Acc.Number) %>%
  separate(amps, into=as.character(seq(8)), sep=",") %>%
  gather(key=concentration_index, value=amplitude, -c(Acc.Number, compound)) %>%
  mutate(amplitude=as.numeric(str_trim(amplitude))) %>%
  write_csv("./PB0072_PD_MCR_PBMC#2_concentration_boxplots_tableau.csv")
#***************************************************************************************
###################################################################################
#Correlation across datasets
d.set2 <- read.csv("~/Documents/Data Analysis/SweLife MCR Analysis/PB0072/MCR_SweLife_1810_PBMC_2_sheet.csv")
d.set3 <- read.csv("~/Documents/Data Analysis/SweLife MCR Analysis/Swelife_PB0074/MCR_SweLife_1810_PBMC_3_sheet.csv")
d.set4 <- read.csv("~/Documents/Data Analysis/SweLife MCR Analysis/PB0075/PB0075 SweLife/MCR_SweLife_1811_PBMC_4_sheet.csv")
# d.set4 <- read.csv(~/Documents/Data Analysis/SweLife MCR Analysis/PB0072/MCR_SweLife_1810_PBMC_2_sheet.csv")
library(plyr)
dset2_var <-d.set2 %>%
  select(Gene,Amp_VCR,Amp_MTX,Amp_MP6,Amp_CytA) %>%
  rename(c("Amp_VCR" ="VCR_2","Amp_MTX" = "MTX_2", "Amp_MP6" = "MP6_2","Amp_CytA" = "CytA_2"))
 
dset3_var <-d.set3 %>%
  select(Gene,Amp_VCR,Amp_MTX,Amp_MP6,Amp_CytA) %>%
  rename(c("Amp_VCR"="VCR_3","Amp_MTX" = "MTX_3", "Amp_MP6" = "MP6_3","Amp_CytA" = "CytA_3"))

dset4_var <-d.set4 %>%
  select(Gene,Amp_VCR,Amp_MTX,Amp_MP6,Amp_CytA) %>%
  rename(c("Amp_VCR" ="VCR_4","Amp_MTX" = "MTX_4", "Amp_MP6" = "MP6_4","Amp_CytA" = "CytA_4"))

dset_common <- dset2_var %>%
  inner_join(dset3_var, by='Gene') %>%
  inner_join(dset4_var, by='Gene')

cormat <- cor(dset_common %>% mutate(Gene=NULL), use="pairwise.complete.obs")
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()  +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#  filter(Gene %in% dset3_var$Gene & Gene %in% dset4_var$Gene)

#++++clean the data.
#The interest is the amplitude for the compds and their corresponding genes 
#Make correlation across datasets per compd for all datasets
#Find the common proteins in all the datasets
gene_fil <- file_ %>%
  filter(Gene =="SERPINE1", del)

file_["Gene" == "SERPINE1", ] <- NULL
 



sessionInfo() 


