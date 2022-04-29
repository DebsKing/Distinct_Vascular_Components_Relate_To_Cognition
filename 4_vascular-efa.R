### 4. Vascular exploratory factor analysis (EFA)

#------------------------------------------------------------------------------#
## Load 

library(GGally) # for correlation plot
library(plyr) # for count().
library(psych) # for the EFA function: fa().
# install.packages("GPArotation")
library(GPArotation) # for rotations.

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis/")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

 mydata <- readRDS("3_mydata_cfa_vascular_output_n668.rds")

#------------------------------------------------------------------------------#
### Inspect latent variables 'BP', 'PP' and 'HR' output from CFA

## Distributions
par(mfrow=c(1,3))

qqnorm(mydata$BP) # normality a bit off
qqnorm(mydata$PP)
qqnorm(mydata$HR) 

hist(mydata$BP)
hist(mydata$PP)
hist(mydata$HR)

hist(mydata$BP, breaks=100) # outliers? 
hist(mydata$PP, breaks=100)
hist(mydata$HR, breaks=100)

# Possible non-normality & outliers in the latent variable BP. 
# Is this a problem that may bias factor formation in EFA?

## Skew and kurtosis
describe(mydata[,cbind( "BP", "PP",  "HR")])  
# acceptable 

#------------------------------------------------------------------------------#
### Inspect observed variables BMI, HRV-lf & HRV-hf

## Distributions
par(mfrow=c(1,3))

qqnorm(mydata$bmiW1)
qqnorm(mydata$ecgwin1_lf) 
qqnorm(mydata$ecgwin1_hf)

hist(mydata$bmiW1)
hist(mydata$ecgwin1_lf)
hist(mydata$ecgwin1_hf)

hist(mydata$bmiW1, breaks=100)
hist(mydata$ecgwin1_lf, breaks=100)
hist(mydata$ecgwin1_hf, breaks=100)

# All have non-normal, positively skewed distributions

## Check range
describe(mydata[,cbind( "bmiW1", "ecgwin1_lf", "ecgwin1_hf")]) 
# Ranges do not span zero; suitable for log transformation to normal distribution.

### Normalise vascular vars with log:
mydata$HRV_lf <- log(mydata$ecgwin1_lf)
mydata$HRV_hf <- log(mydata$ecgwin1_hf)
mydata$BMI <- log(mydata$bmiW1)

### Check distributions, skew and kurtosis after log 
qqnorm(mydata$BMI)
qqnorm(mydata$HRV_lf)
qqnorm(mydata$HRV_hf)

hist(mydata$BMI, breaks=100) 
hist(mydata$HRV_lf, breaks=100)
hist(mydata$HRV_hf, breaks=100)

describe(mydata[,cbind( "bmiW1", "ecgwin1_lf", "ecgwin1_hf")]) 

# Distributions are imperfect but much improved after log transformation. 
# Kurtosis remains high for HRV_lf.

#------------------------------------------------------------------------------#
### Standardise (aka Z-score)

# Unlike in Lavaan's cfa(), EFA with fa() does not z-score input variables. 
# Therefore specify manually here.

mydata$BP<- scale(mydata$BP, center=TRUE, scale=TRUE)
mydata$PP<- scale(mydata$PP, center=TRUE, scale=TRUE)
mydata$HR<- scale(mydata$HR, center=TRUE, scale=TRUE)

mydata$HRV_lf<- scale(mydata$HRV_lf, center=TRUE, scale=TRUE)
mydata$HRV_hf<- scale(mydata$HRV_hf, center=TRUE, scale=TRUE)
mydata$BMI<- scale(mydata$BMI, center=TRUE, scale=TRUE)

## Check it worked: mean=0, sd=1. 
describe(mydata[,cbind( "BP", "PP",  "HR")])  
describe(mydata[,cbind( "BMI", "HRV_lf", "HRV_hf")]) 

#------------------------------------------------------------------------------#
### Correlation plot of EFA input

### Correlation plot of LVFs, CDS & Age
# (vascular variables are Z-scored)
par(mfrow=c(1,1))
corr_pca <- data.frame( mydata$BP, mydata$PP, mydata$HR,
                        mydata$BMI, mydata$HRV_lf, mydata$HRV_hf,
                        mydata$Age)
# Rename cols
colnames(corr_pca)<-c( "BP", "PP", "HR", "BMI", "HRV_LF", "HRV_HF", "Age")

# Detailed corrplot:
plotcor<-ggpairs(corr_pca,digits=3)
plotcor

# More detailed plots:
my_fn <- function (data, mapping, method = "loess"){
  p <- ggplot( data = data, mapping = mapping) + 
    geom_jitter(height = 0.2, width = 0.2,  # jitter to visualise discrete data bit easier
                alpha=0.4) +  # greyscale to visualise data intensity
    geom_smooth (method=lm, color = "blue", fill="blue" , # add association line
                 se=FALSE) # removes confidence intervals, to keep plot uncluttered
  p
}
# Plot resulting graph:
ggpairs( corr_pca, lower = list (continuous =  my_fn  ))

## Save plot
# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Outputs/Plots")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Outputs/Plots")

ggsave(plotcor,filename = "Corrplot-VFs-Age.png")

rm(corr_pca, plotcor)

#------------------------------------------------------------------------------#
### Fit EFA models 

### One factor EFA
fit_efa_1 <- fa(data.frame(mydata$BP , mydata$PP, mydata$HR,
                          mydata$BMI , mydata$HRV_lf , mydata$HRV_hf) ,
                          nfactors= 1 ,  # would give a summary score, eg QRisk or Framingham Stroke
                          scores=TRUE, missing=TRUE) # imputes missing

### Two factor EFA 
# fit_fa_2 <- fa(data.frame(mydata$BP , mydata$PP, mydata$HR,
#                           mydata$BMI , mydata$HRV_lf , mydata$HRV_hf) ,
#                nfactors= 2 , scores=TRUE, missing=TRUE)
# Error. Model not identified.
# The estimated weights for the factor scores are probably incorrect. 

### Three factor EFA
fit_efa_3 <- fa(data.frame(mydata$BP , mydata$PP, mydata$HR,
                        mydata$BMI , mydata$HRV_lf , mydata$HRV_hf) ,
                        nfactors= 3 ,
                        scores=TRUE, missing=TRUE) # imputes missing
### Four factor EFA
# fit_efa_4 <- fa(data.frame(mydata$BP , mydata$PP, mydata$HR,  
#                            mydata$BMI , mydata$HRV_lf , mydata$HRV_hf) , 
#                 nfactors= 4 , 
#                 scores=TRUE, missing=TRUE)

#------------------------------------------------------------------------------#
### Compare models 

anova(fit_efa_1, fit_efa_3) 
# For factor analytic models, anova compares Chi Square

### Report values for 1-factor model
round(c(fit_efa_1$dof, fit_efa_1$STATISTIC, fit_efa_1$PVAL), digits=3)
# Chi sq (9) = 340.924 , p<0.001.

### Report values for 3-factor model
round(c(fit_efa_3$dof, fit_efa_3$STATISTIC, fit_efa_3$PVAL), digits=3)
# Chi sq (0) = 9.314, p=NA.

### Remove losing model
rm(fit_efa_1)

#------------------------------------------------------------------------------#
### Inspect winning model

### Diagram 
par(mfrow=c(1,1))
fa.diagram(fit_efa_3,   digits=3) 
# Shows largest factor loading per observed variable for simplicity of 
# presentation & understanding.


### Factor loadings & inter-factor correlations
# Including residual variance between observed and latent variables that are the 
# primary loadings.
print(fit_efa_3, digits=3)
# Report in sup mat table

#------------------------------------------------------------------------------#
### Subject scores

## Extract subject scores
scores_efa <- data.frame (fit_efa_3$scores)
# head(fit_efa_3$scores)
# summary(fit_efa_3$scores)

## Check scores for missing data
count (scores_efa$MR1 != "NaN" | scores_efa$MR2 != "NaN" | scores_efa$MR3 != "NaN") 
count(scores_efa$MR1 == "NaN" | scores_efa$MR2 == "NaN" | scores_efa$MR3 == "NaN") 
# EFA output is complete for n=668; imputation worked. 

## Input subject scores to main dataframe
# Re-order for interpretation in paper 
# LVF1 is total blood pressure
mydata$LVF1 <- scores_efa$MR2
# LVF2 is pulse pressure
mydata$LVF2 <- scores_efa$MR3
# LVF3 is heart rate variability
mydata$LVF3 <- scores_efa$MR1

### Correlate EFA factors with Age
cor.test(mydata$LVF1, mydata$Age) # r=0.3383757   , p<0.001
cor.test(mydata$LVF2, mydata$Age) #  r=0.4665545, p<0.001
cor.test(mydata$LVF3, mydata$Age) # r= -0.6122282 , p<0.001
# Could also do r and pair correlation tests - see Deliah's scripts. 

#-----------------------------------------------------------------------------#
### Save model

# Save winning model's fit structure of results and loadings

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

saveRDS(fit_efa_3, "fit_efa_three_factor.rds")

#-----------------------------------------------------------------------------#
### Save data 

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

# Save as RDS for use in future scripts. 
saveRDS(mydata, file="4_mydata_efa_vascular_output_n668.rds")

# Save as CSV for viewing externally.
write.csv(mydata, "4_mydata_efa_vascular_output_n668.csv")


# Check if CSV looses dp precision, which happens sometimes.
mydata <- read.csv("4_mydata_efa_vascular_output_n668.csv")
# Check for decimal places in some vars:
head(mydata$bmiW1)
head(mydata$ecgwin1_hf)



#-----------------------------------------------------------------------------#
rm(mydata)
rm(mydata, fit_efa_3, scores_efa)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
### ENDS




