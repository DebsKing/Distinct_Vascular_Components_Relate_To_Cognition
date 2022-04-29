### 5. CFA on cognitive factors

#------------------------------------------------------------------------------#
### Load 

# install.packages("lavaan")
library(lavaan)

# install.packages("GGally")
library(GGally) # correlation plot

# install.packages("tidyverse") # for ggplot
library(tidyverse)

# install.packages("ggplot2")
library(ggplot2)


library(plyr) # count() function
library(psych) 

# install.packages("semPlot")
library(semPlot) 

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis/")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

mydata <- readRDS("2_mydata_preprocessed_cognitive_n678.rds")

#------------------------------------------------------------------------------#
### Inspect data 

head(mydata$Proverbs_Score)
head(mydata$STW_total)
head(mydata$Cattell_SubScore1)
# Discrete data

# Check for normality
par(mfrow=c(1,2))
qqnorm(mydata$Proverbs_Score)
qqnorm(mydata$STW_total)

par(mfrow=c(2,2))
qqnorm(mydata$Cattell_SubScore1)
qqnorm(mydata$Cattell_SubScore2)
qqnorm(mydata$Cattell_SubScore3)
qqnorm(mydata$Cattell_SubScore4)

# Check for skew
par(mfrow=c(1,2))
hist(mydata$Proverbs_Score)
hist(mydata$STW_total)

par(mfrow=c(2,2))
hist(mydata$Cattell_SubScore1)
hist(mydata$Cattell_SubScore2)
hist(mydata$Cattell_SubScore3)
hist(mydata$Cattell_SubScore4)

describe(mydata[,cbind("Proverbs_Score", "STW_total",
                       "Cattell_SubScore1", "Cattell_SubScore2", 
                       "Cattell_SubScore3", "Cattell_SubScore4"
                       )])

# Negative skews and differing extents of kurtosis. Non-continuous data. Do not log. 

# Also no need to z-score, because observed variables are standardised 
# within the Lavaan function: cfa(..., std.ov=T).

#---------------------------------------------------------------------------------------------------#
### Correlation plot of observed cognitive vars & Age

par(mfrow=c(1,1))
mydata_cor <- data.frame( mydata$STW_total, mydata$Proverbs_Score,
                        mydata$Cattell_SubScore1, mydata$Cattell_SubScore2,
                        mydata$Cattell_SubScore3, mydata$Cattell_SubScore4,
                        mydata$Age)
# Rename cols
colnames(mydata_cor)<-c( "STW", "Proverbs",
                       "Cattell1", "Cattell2", "Cattell3", "Cattell4",
                       "Age")


# Detailed corrplot:
plotcor<-ggpairs(mydata_cor,digits=3)
plotcor

### Editing plot
# Useful resource: https://ggobi.github.io/ggally/reference/ggpairs.html. 
# Co-author suggested (Jan-2022), add jitter, grayscale and 'lm' line of best fit for associations.
# Add minimal jitter bcse it is important that data visualisation is accurate, showing data to be discrete. 

# ggpairs(mydata_cor, 
#         digits=3,# Decimal places to 3dp
#         lower= list(continuous= wrap ("smooth",     # "points", ## smooth adds association line; points is just the data points.
#                                      method = "lm",  # i.e. "lm" is linear, rather than method="loess" for non-linear. 
#                                       alpha=0.4, # visualise grayscale of data intensity
#                                       position = position_jitter (height = 0.2, width = 0.2) # jitter
#                                      ))) 

## Editing plot with a function funciton to plot with regression lines: 
my_fn <- function (data, mapping, method = "loess"){
  p <- ggplot( data = data, mapping = mapping) + 
    geom_jitter(height = 0.2, width = 0.2,  # jitter to visualise discrete data bit easier
                alpha=0.4) +  # greyscale to visualise data intensity
    geom_smooth (method=lm, color = "blue", fill="blue" , # add association line
                 se=FALSE) # removes confidence intervals, to keep plot uncluttered
  p
}
# Plot resulting graph:
ggpairs( mydata_cor, 
         lower = list (continuous =  my_fn  ),
         upper = list(continuous = wrap( "cor", method = "spearman")))

## Save plot
# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Outputs/Plots")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Outputs/Plots")

ggsave(plotcor,filename = "Corrplot_Cognitive-Obs.png")

#------------------------------------------------------------------------------#
### Run 2-factor CFA on cognitive data

### Develop CFA model
# 2-factor driven by theory and previous Cam-CAN papers. 

model_CF <- 'LCF1 =~ Proverbs_Score + STW_total
              LCF2 =~ Cattell_SubScore1 + Cattell_SubScore2 + Cattell_SubScore3 + Cattell_SubScore4
             '

### Fit CFA model
fit_CF <- cfa(model_CF, data=mydata, 
              std.ov=T,  # std.ov - observed variables are standardized
              std.lv=T,  # std.lv - variances of LVs fixed to 1
              missing='fiml', estimator='mlr') # estimate missing data.

### Model fit
summary(fit_CF) # Chi-sq (8) = 19.237, p=0.014.
# The summary shows the estimates which are the path loadings.

### Plot model
# par(mfrow=c(1,1))
# semPaths(fit_CF, "par",
#          sizeMan = 15, sizeInt = 30, sizeLat=10, fade=FALSE,
#          edge.label.cex=1.5,# Make edge labels readable
#          residuals = FALSE,  # Remove the residual variance, within each variable
#          intercepts = FALSE,# Remove intercepts
#          rotation = 3, nCharNodes=9)

### Extract subject scores per factor:
CF <- as.data.frame( predict(fit_CF) )

### Check imputation worked 
count(CF$LCF1 == "NaN" |   CF$LCF2 == "NaN")

### Put latent variable subject scores into main dataframe:
mydata <- data.frame (mydata[,], CF$LCF1, CF$LCF2)

### Rename cols
names(mydata)[names(mydata) == "CF.LCF1"] <- "LCF1"
names(mydata)[names(mydata) == "CF.LCF2"] <- "LCF2"

rm(CF, fit_CF, model_CF)

#------------------------------------------------------------------------------#
### Sanity checks 

## Check that the factors show the anticipated relationships with age

# LCF1 representing crystallised intelligence should increase slightly with age
ggplot(data = mydata, aes(x = Age, y =LCF1)) +
  geom_point(alpha=0.5) +  geom_smooth(method="lm",se=TRUE) 

plot(mydata$Age, mydata$LCF1)

# LCF2 representing fluid intelligence should decline with age
ggplot(data = mydata, aes(x = Age, y =LCF2)) +
  geom_point(alpha=0.5) +  geom_smooth(method="lm",se=TRUE) 

plot(mydata$Age, mydata$LCF2)

# #------------------------------------------------------------------------------#
## Notes

# do not repeat z-scoreing inbetween CFA and calculation of CDS; because
# McDonough et al., 2016 did not do this in original paper.
# And I z-score at a later stage, ie before linear regression.

#------------------------------------------------------------------------------#
### Calculate the ability discrepancy score per participant

# DS = CI (proverbs and STW) - FI (cattell)
mydata$DS <- mydata$LCF1 - mydata$LCF2

describe(mydata$DS)
# mean=0, SD=1.01.  
# range: -3.66 - 3.2.
# Range in McDonough et al., 2016: 2.97-2.32. 

# count(mydata$DS>0)
# McDonugh et al., 2016 included only participants with DS>0. I will include all.

#------------------------------------------------------------------------------#
### Sanity checks on DS

# correlates positively with age

cor.test(mydata$DS,mydata$Age) # r= 0.6844341 , p<0.001.

ggplot(data = mydata, aes(x = Age, y =DS)) +
  geom_point(alpha=0.5) +  geom_smooth(method="lm",se=TRUE) 

plot(mydata$Age, mydata$DS)

#------------------------------------------------------------------------------#
### Save data

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis/")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

# Save as RDS for use in future scripts. 
saveRDS(mydata, file="5_mydata_discrepancy-score_output_n678.rds")

# Save as CSV for viewing externally.
write.csv(mydata, "5_mydata_discrepancy-score_output_n678.csv")

rm(mydata, corr_pca, plotcor)
