
###  Multigroup CFA on vascular model

# To explore age group invaraince across cross-sectional sample. 
# For supplementary section of the paper. 

#------------------------------------------------------------------------------#
### Load packages

library(lavaan) # for SEM
library(semPlot) # for plotting SEM diagram. 

# install.packages("rgl" ) # needed for 'qpcr'. 
library(rgl)

# install.packages("qpcR") # install.packages('qpcR)
library(qpcR) # for akaike weights # 

library(plyr)

library(ggplot2); library(GGally)

#------------------------------------------------------------------------------#
### Load data

setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")
mydata <- readRDS("3_mydata_cfa_vascular_output_n668.rds")

#------------------------------------------------------------------------------#
### Replicate preprocessing steps from script 4 of main analysis, on vascular EFA

## Normalise vascular vars with log:
mydata$HRV_lf <- log(mydata$ecgwin1_lf)
mydata$HRV_hf <- log(mydata$ecgwin1_hf)
mydata$BMI <- log(mydata$bmiW1)

## Standardise (aka Z-score)
mydata$BP<- scale(mydata$BP, center=TRUE, scale=TRUE)
mydata$PP<- scale(mydata$PP, center=TRUE, scale=TRUE)
mydata$HR<- scale(mydata$HR, center=TRUE, scale=TRUE)
mydata$HRV_lf<- scale(mydata$HRV_lf, center=TRUE, scale=TRUE)
mydata$HRV_hf<- scale(mydata$HRV_hf, center=TRUE, scale=TRUE)
mydata$BMI<- scale(mydata$BMI, center=TRUE, scale=TRUE)

#------------------------------------------------------------------------------#
### Split age into two groups 

# Rank order by Age 
attach(mydata)
mydata <- mydata[order(Age),]
print(mydata$Age)

# Create a variable for age groups labelled "Age_1" for young, "Age_2" for old 
# Split participants equally per group, n=334 in each.

mydata$Age_Group_2 <- "Age_1"
mydata[335:668,]$Age_Group_2  <- "Age_2"

# Check
print(mydata$Age_Group_2 )
print(mydata[330:340,]$Age_Group_2 )

count(mydata$Age_Group_2=="Age_1")
count(mydata$Age_Group_2=="Age_2")

# Report age range
range(subset(mydata, Age_Group_2=="Age_1", select=("Age")))
range(subset(mydata, Age_Group_2=="Age_2", select=("Age")))

#------------------------------------------------------------------------------#
## Specify the CFA, based on previous 3-factor EFA output

model_cfa <- 'LV1 =~ BP + BMI
              LV2 =~ PP + HR
              LV3 =~ HRV_hf + HRV_lf
            '
#------------------------------------------------------------------------------#
### Fit multigroup CFA, over two age groups

fit_constrained <- cfa(model_cfa,
             data=mydata,
             estimator='MLR', missing='fiml',
             group="Age_Group_2" ,
             group.equal = c("lv.variances","intercepts","lv.covariances",
                             "residual.covariances", "loadings", "residuals"))

fit_loadings_free <- cfa(model_cfa,
                       data=mydata,
                       estimator='MLR', missing='fiml',
                       group="Age_Group_2" ,
                       group.equal = c("lv.variances","intercepts","lv.covariances",
                                       "residual.covariances",  "residuals"))

#------------------------------------------------------------------------------#
### Check MG-CFA model fit:

# Models compared on change in robust CFI, which is argued to be superior
# to X^2, likelihood ratio (i.e. anova), and other goodness of fit indices, 
# because it is independent of sample size. More details in 
# Cheung and Rensvold, 2002 paper on evaluating goodness of fit. 

summary(fit_constrained, fit.measures=TRUE, standardized=T)

summary(fit_loadings_free, fit.measures=TRUE, standardized=T)

#------------------------------------------------------------------------------#
### Plot free to vary model
  
par(mfrow=c(1,2))
semPaths(fit_loadings_free,  
         what='est',
         rotation = 2, # default rotation = 1 with four options
         curve = 1, # pull covariances' curves out a little
         nCharNodes = 10,
         nCharEdges = 10, # don't limit variable name lengths
         sizeMan = 10, # font size of manifest variable names
         style = "lisrel", # single-headed arrows vs. # "ram"'s 2-headed for variances
         edge.label.cex=2, # size of numbers on path arrows
         curvePivot = TRUE,
         fade=FALSE,
         residuals = FALSE, intercepts = FALSE,
         sizeInt = 10,sizeLat=10)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
### ENDS
