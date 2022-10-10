# Splitting vascular EFA per two age groups 

# 2022-10-04

#------------------------------------------------------------------------------#
## Load 

library(GGally) # for correlation plot
library(plyr) # for count().
library(psych) # for the EFA function: fa().
# install.packages("GPArotation")
library(GPArotation) # for rotations.

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
### Fit and compare models on Young Age Group 

## One factor EFA 
fit_efa_1_y <- fa(data.frame(subset (mydata, Age_Group_2=="Age_1", 
                                   select=c (BP, PP, BMI, HR, HRV_lf, HRV_hf ))) ,
                nfactors= 1 ,  # would give a summary score, eg QRisk or Framingham Stroke
                scores=TRUE, missing=TRUE) # imputes missing

## Two factor EFA 
fit_efa_2_y <- fa(data.frame(subset (mydata, Age_Group_2=="Age_1",
                                    select=c (BP, PP, BMI, HR, HRV_lf, HRV_hf ))) ,
               nfactors= 2 , scores=TRUE, missing=TRUE)
# Error. Model not identified. The estimated weights for the factor scores are probably incorrect.

## Three factor EFA
fit_efa_3_y <- fa(data.frame(subset (mydata, Age_Group_2=="Age_1", 
                                   select=c (BP, PP, BMI, HR, HRV_lf, HRV_hf ))) ,
                nfactors= 3 ,
                scores=TRUE, missing=TRUE) # imputes missing
# Warning: The estimated weights for the factor scores are probably incorrect.


### Compare models on Age Group 1

anova(fit_efa_1_y, fit_efa_2_y) 

anova(fit_efa_2_y, fit_efa_3_y) 

# The 3-factor EFA wins in the young age group


#------------------------------------------------------------------------------#
### Fit and compare models on Old Age Group

## One factor EFA 
fit_efa_1_o <- fa(data.frame(subset (mydata, Age_Group_2=="Age_2", 
                                     select=c (BP, PP, BMI, HR, HRV_lf, HRV_hf ))) ,
                  nfactors= 1 ,  # would give a summary score, eg QRisk or Framingham Stroke
                  scores=TRUE, missing=TRUE) # imputes missing

## Two factor EFA 
fit_efa_2_o <- fa(data.frame(subset (mydata, Age_Group_2=="Age_2",
                                     select=c (BP, PP, BMI, HR, HRV_lf, HRV_hf ))) ,
                  nfactors= 2 , scores=TRUE, missing=TRUE)

# Error. Model not identified. The estimated weights for the factor scores are probably incorrect.

## Three factor EFA
fit_efa_3_o <- fa(data.frame(subset (mydata, Age_Group_2=="Age_2", 
                                     select=c (BP, PP, BMI, HR, HRV_lf, HRV_hf ))) ,
                  nfactors= 3 ,
                  scores=TRUE, missing=TRUE) # imputes missing


### Compare models on Age Group 2

anova(fit_efa_1_o, fit_efa_2_o) 

anova(fit_efa_2_o, fit_efa_3_o) 

# The 3-factor EFA wins in the old age group


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
### ENDS
