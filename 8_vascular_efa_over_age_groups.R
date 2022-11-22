### EFA over mutiple age sub-groups

# Run in November 2022
# Revisions on manuscript
# Deborah King


### Description of approach 

# The vascular EFA model was originally run over the whole sample, which has a broad age range.
# It is therefore possible that the mean age trends within the vascular observations
# input to EFA may influence the EFA model structure. 
# To investigate the model's robustness with age, we re-ran the EFA over sub-groups split
# by age. We extracted the resulting subject factor scores. 
# We stacked these factor scores, putting them back into the original dataframe. 
# We then input the factor scores to the regression models. We asked, are the regression
# results consistent with those from the EFA, as run originally over the whole sample?
# We showed that they are. This showed that the EFA mdoel structure is 
# fairly robust to age.


### Notes on approach

# Before EFA, we need to standardise the vascular observations. This refers to
# centering the mean = 0, also known as de-meaning, and scaleing the SD = 1 . 

# In the whole sample EFA, we z-scored across the whole sample (n=668), immediately
# before running EFA. 

# In the multi-group EFA, we z-scored across each individual sub-group, as split by age. 

# However, by z-scoring across sub-groups we removed the mean age trends
# within the vascular observations. These were important for regression.
# Therefore, after z-scoring and modelling EFA on sub-groups, we needed to re-inject the 
# age-group means. 
# To do this, we extracted from the whole sample (n=668), the mean of a given
# latent vascular factor, where the mean was calculated using only those participants
# corresponding to a given age specific sub-group. The resulting mean was then added
# to the factor scores produced in the EFA run on that sub-group. 
# This was repeated for all 3 latent factors, and all sub-groups.


#------------------------------------------------------------------------------#
## Load 

rm(list=ls())

## packages
library(GGally) # for correlation plot
library(ggplot2)
library(plyr) # for count()
library(psych) # for describe()
library(sjPlot); library(sjmisc); library(sjlabelled) # lm results in table
library(stats) # for p.adjust(), to correct for multiple comparisons.
library(MASS) # for rlm(), to run robust linear models.
library(sfsmisc)

#------------------------------------------------------------------------------#
### Load data

rm(list=ls())

# setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")
# mydata <- readRDS("3_mydata_cfa_vascular_output_n668.rds")

# Load data output from script "3_vascular_cfa.R". 

## Normalise vascular vars which were not input to CFA models with log:
mydata$HRV_lf <- log(mydata$ecgwin1_lf)
mydata$HRV_hf <- log(mydata$ecgwin1_hf)
mydata$BMI <- log(mydata$bmiW1)


#------------------------------------------------------------------------------#
### Age groups

# Here, mimicking the age groups in the demographics table and figure 4 of initial manuscript. 

# # count(mydata$Age_Group)
# 1 Middle  311
# 2    Old  199
# 3  Young  158

#------------------------------------------------------------------------------#
### Fit EFA models on whole sample, and on young, iddle and old sub-groups

# note - scale vars within the EFA function.

# Variables to run EFA over: 
varname = c("BP","PP","HR","BMI","HRV_lf","HRV_hf")

Age_Group=mydata$Age_Group

# Whole sample n=668
fit_efa3_all <- fa(data.frame(scale(mydata[,varname], center=TRUE, scale=TRUE)) ,
                nfactors= 3 ,
                scores=TRUE, missing=TRUE)

## Young
fit_efa3_y <- fa(data.frame(scale(mydata[Age_Group=="Young",varname], center=TRUE, scale=TRUE))  ,
                  nfactors= 3 ,
                  scores=TRUE, missing=TRUE) # imputes missing
# Warning: The estimated weights for the factor scores are probably incorrect.

## Middel
fit_efa3_m <- fa(data.frame(scale(mydata[Age_Group=="Middle",varname], center=TRUE, scale=TRUE))  ,
                 nfactors= 3 ,
                 scores=TRUE, missing=TRUE) # imputes missing
# Warning: The estimated weights for the factor scores are probably incorrect.

## Old
fit_efa3_o <- fa(data.frame(scale(mydata[Age_Group=="Old",varname], center=TRUE, scale=TRUE))  ,
                  nfactors= 3 ,
                  scores=TRUE, missing=TRUE) # imputes missing

### Diagram 
par(mfrow=c(1,4))
fa.diagram(fit_efa3_all, digits=3, main="Whole Sample EFA")
fa.diagram(fit_efa3_y,   digits=3, main="Young Age Sub-Group") 
fa.diagram(fit_efa3_m,   digits=3, main="Middle Age Sub-Group") 
fa.diagram(fit_efa3_o,   digits=3, main="Old Age Sub-Group") 


#------------------------------------------------------------------------------#
### Subject scores

## Extract subject scores
scores_all <- data.frame (fit_efa3_all$scores)
scores_y <- data.frame (fit_efa3_y$scores)
scores_m <- data.frame (fit_efa3_m$scores)
scores_o <- data.frame (fit_efa3_o$scores)

## Input subject scores to main dataframe

# Note that factors MR1-3 are numbered differently across the whole sample and sub-group EFAs
# I input to dataframe based on the factors reported in the paper, ie 
# LVF1 expresses predominantly BP; LVF2 expresses predominantly pulse pressure; 
# and LVF3 expresses HRV.

# Factor scores on latent vascular factors (LVF) across all participants:
mydata$LVF1_all <- scores_all$MR2 # blood
mydata$LVF2_all <- scores_all$MR3 # pulse
mydata$LVF3_all <- scores_all$MR1 # HRV

# Create empty multi-group factors
mydata$LVF1_mg <- NaN
mydata$LVF2_mg <- NaN
mydata$LVF3_mg <- NaN

# Factor scores on latent vascular factors (LVF) for the multigroup EFAs, stacked:


# Young
mydata[Age_Group=="Young",]$LVF1_mg <- scores_y$MR3 # blood
mydata[Age_Group=="Young",]$LVF2_mg <- scores_y$MR1 # pressure
mydata[Age_Group=="Young",]$LVF3_mg <- scores_y$MR2 # hrv

# Middle
mydata[Age_Group=="Middle",]$LVF1_mg <- scores_m$MR3 # blood
mydata[Age_Group=="Middle",]$LVF2_mg <- scores_m$MR2 # pulse 
mydata[Age_Group=="Middle",]$LVF3_mg <- scores_m$MR1 # hrv

# Old
mydata[Age_Group=="Old",]$LVF1_mg <- scores_o$MR3 # ... mostly heart rate
mydata[Age_Group=="Old",]$LVF2_mg <- scores_o$MR2 # pulse ... also blood
mydata[Age_Group=="Old",]$LVF3_mg <- scores_o$MR1 # hrv

#---------------------------------------------------------------#
### Investigate whether the factor scores output from EFA were automatically z-scored

describe(mydata[,cbind( "LVF1_all", "LVF2_all",  "LVF3_all","LVF1_mg", "LVF2_mg",  "LVF3_mg")])  

#---------------------------------------------------------------#
### Add the group mean of each extract factor score, to the factor scores

## Add means of the the whole sample EFA scores, for corresponding participants,
# in main dataframe 'mydata'

# first, create empty variables
mydata$LVF1_mg_mean  <- NaN
mydata$LVF2_mg_mean <- NaN
mydata$LVF3_mg_mean <- NaN

# For young group:
mydata[Age_Group=="Young",]$LVF1_mg_mean  <- (mydata[Age_Group=="Young",]$LVF1_mg + (mean(mydata[Age_Group=="Young",]$LVF1_all)))
mydata[Age_Group=="Young",]$LVF2_mg_mean <- (mydata[Age_Group=="Young",]$LVF2_mg + (mean(mydata[Age_Group=="Young",]$LVF2_all)))
mydata[Age_Group=="Young",]$LVF3_mg_mean <- (mydata[Age_Group=="Young",]$LVF3_mg + (mean(mydata[Age_Group=="Young",]$LVF3_all)))

# For middel group:
mydata[Age_Group=="Middle",]$LVF1_mg_mean  <- (mydata[Age_Group=="Middle",]$LVF1_mg + (mean(mydata[Age_Group=="Middle",]$LVF1_all)))
mydata[Age_Group=="Middle",]$LVF2_mg_mean <- (mydata[Age_Group=="Middle",]$LVF2_mg + (mean(mydata[Age_Group=="Middle",]$LVF2_all)))
mydata[Age_Group=="Middle",]$LVF3_mg_mean <- (mydata[Age_Group=="Middle",]$LVF3_mg + (mean(mydata[Age_Group=="Middle",]$LVF3_all)))

# For Old group:
mydata[Age_Group=="Old",]$LVF1_mg_mean  <- (mydata[Age_Group=="Old",]$LVF1_mg + (mean(mydata[Age_Group=="Old",]$LVF1_all)))
mydata[Age_Group=="Old",]$LVF2_mg_mean <- (mydata[Age_Group=="Old",]$LVF2_mg + (mean(mydata[Age_Group=="Old",]$LVF2_all)))
mydata[Age_Group=="Old",]$LVF3_mg_mean <- (mydata[Age_Group=="Old",]$LVF3_mg + (mean(mydata[Age_Group=="Old",]$LVF3_all)))


## Check it worked
count(mydata$LVF2_mg==mydata$LVF2_mg_mean) # all differ
describe(mydata[,cbind( "LVF1_mg", "LVF1_mg_mean","LVF2_mg",  "LVF2_mg_mean","LVF3_mg",   "LVF3_mg_mean")])  

#-----------------------------------------------------------------------------#
### Look at how much the MG-EFA+means factors differ from those in 
# the whole group EFA

par(mfrow=c(1,1))
corr_pca <- data.frame( mydata$LVF1_all, mydata$LVF2_all,mydata$LVF3_all,
                        mydata$LVF1_mg_mean, mydata$LVF2_mg_mean, mydata$LVF3_mg_mean,
                        mydata$Age) 

# Rename cols
colnames(corr_pca)<-c( "LVF1_all", "LVF2_all","LVF3_all", 
                       "LVF1_mg_mean", "LVF2_mg_mean", "LVF3_mg_mean",
                       "Age")

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
ggpairs( corr_pca, lower = list (continuous =  my_fn  ),
         title="Factor Scores From All vs Multi-Group EFA ") # uses Pearson.

#---------------------------------------------------------------------------#
### Stack data and Plot for supplementary material of paper:

## Stack LVF scores from EFA on whole sample
mydata_all<-(data.frame(mydata$LVF1_all, mydata$LVF2_all, mydata$LVF3_all))
colnames(mydata_all)
colnames(mydata_all) <- c("Vacsular Factor 1", 
                          "Vacsular Factor 2", 
                          "Vacsular Factor 3")
mydata_all<-stack(mydata_all)
head(mydata_all)
colnames(mydata_all) <- c("Scores_all", "LVF")


## Stack LVF scores from EFA on 2 groups
mydata_mg<-(data.frame(mydata$LVF1_mg_mean, mydata$LVF2_mg_mean, mydata$LVF3_mg_mean))
colnames(mydata_mg)
colnames(mydata_mg) <- c("LVF1_mg_mean", "LVF2_mg_mean", "LVF3_mg_mean")
mydata_mg<-stack(mydata_mg)
head(mydata_mg)
colnames(mydata_mg) <- c("Scores_mg", "Factors_mg")


## Stack the 2 stacked dataframes
mydata_stck<- data.frame(mydata_all, mydata_mg)
head(mydata_stck)


## Plot

ggplot(mydata_stck, aes(Scores_all, Scores_mg)) + 
  geom_point(alpha=0.5)+
  facet_grid(cols=vars(LVF))+
  stat_smooth(method="lm",se=FALSE) +
  theme(text=element_text(size=21),
        strip.text.x=element_text(size=22)) +
  # scale_x_continuous(breaks=c(-4,-2,0,2,4,6))+
  labs(
    x="Whole Sample Factor Scores",
    y="Multigroup Factor Scores") 

## R values
cor.test(mydata$LVF1_all, mydata$LVF1_mg_mean )
cor.test(mydata$LVF2_all, mydata$LVF2_mg_mean )
cor.test(mydata$LVF3_all, mydata$LVF3_mg_mean )
# 

### Factor scores can be saved and input to regression models, in place of factors
# produced on EFA run over the whole sample, and as in script 7. 