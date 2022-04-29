#### 7. Multiple linear regressions

#------------------------------------------------------------------------------#
### Load

## packages
library(GGally) # for correlation plot
library(ggplot2)
library(plyr) # for count()
library(psych) # for describe()
library(sjPlot); library(sjmisc); library(sjlabelled) # lm results in table
library(stats) # for p.adjust(), to correct for multiple comparisons.
library(MASS) # for rlm(), to run robust linear models.
library(sfsmisc)

## Data
# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis/")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")
mydata <- readRDS("6_mydata_for_regressions_n655.rds")

#------------------------------------------------------------------------------#
### Preliminary checks

# 'mydata' was created by merging 2 dataframes, in script 6. 
# The output dataframe contains the IDs from both input dataframes. 
# Check that they match up exactly.

sum(mydata$SubCCIDc == mydata$SubCCIDc_c) # 665
sum(mydata$SubCCIDc != mydata$SubCCIDc_c) # 0

# Check that the data contains vars of interest as expected
#tail(colnames(mydata))

# Check that data is complete for variables of interest
count(mydata$DS != "NaN")  
count(mydata$LVF1 != "NaN" | mydata$LVF2 != "NaN" | mydata$LVF3 != "NaN")

#------------------------------------------------------------------------------#
### Prepare data

## Z-score numerical variables
# Name with 'z' so can easily check linear models are specified with z-scored variables.

mydata$LVF1z <- scale(mydata$LVF1, center=TRUE, scale=TRUE)
mydata$LVF2z <- scale(mydata$LVF2, center=TRUE, scale=TRUE)
mydata$LVF3z <- scale(mydata$LVF3, center=TRUE, scale=TRUE)
mydata$DSz <-scale(mydata$DS, scale=TRUE, center=TRUE)
mydata$Agez <-scale(mydata$Age, scale=TRUE, center=TRUE)

# sex, education and medications are categorical, therefore not z-scored. 

### check z-score. It should have SD=1. 
#describe(mydata[,cbind("LVF1z", "LVF2z", "LVF3z")])
#describe(mydata[,cbind("Agez", "DSz")])

## Education
# Use categorical variable.
mydata$Education <- mydata$qual_Cognition
#head(mydata$Education)


#------------------------------------------------------------------------------#
### Correlation plot of LVFs, DS & Age
# (LVFs and DS are Z-scored, but not age which has useful units for this plot)

par(mfrow=c(1,1))
corr_pca <- data.frame( mydata$LVF1, mydata$LVF2, mydata$LVF3,
                        mydata$DS, mydata$Age)
# Rename cols
colnames(corr_pca)<-c( "LVF1", "LVF2", "LVF3",
                       "Discrepancy",
                       "Age" )
# Detailed corrplot:
plotcor<-ggpairs(corr_pca,digits=3)
plotcor

-----------------------------------------------------------------------------#
  ### Correlation plot of EFA output
  
par(mfrow=c(1,1))
corr_pca <- data.frame( mydata$LVF1, mydata$LVF2, mydata$LVF3,
                        mydata$DS,
                        mydata$Age)
# Rename cols
colnames(corr_pca)<-c( "LVF1", "LVF2", "LVF3", 
                       "Discrepancy", 
                       "Age")

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
setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Outputs/Plots/")
ggsave(plotcor,filename = "Corrplot_LVFs_Discrepancy_Age.png")

rm(plotcor, corr_pca)

# #---------------------------------------------------------------------------------------------------#
# ### Model 1 with ordinary least squares regression
# 
# lm.1.ols <-lm(DSz ~ LVF1z +  LVF2z +  LVF3z + Sex + Education,data=mydata)
# 
# ### Check assumptions
# 
# ## Residuals v fitted
# # Checks the linearity of the data. Aim: no pattern & approx horizontal red line.
# plot(lm.1.ols,1) # acceptable
# 
# ## QQ for normal distribution
# plot(lm.1.ols,2) # outliers 
# 
# ## Spread-location plot
# # Checks homogeneity of variance. Aim: horizontal line & equally spread points, 
# # showing that the variance of residuals does not change in line with the predictors.
# plot(lm.1.ols,3) # acceptable
# 
# ### Outliers and influence.
# # Outliers are okay if they are the tail end of a normal distribution, but not if
# # they have undue influence and differ from the pattern of the data.
# 
# ## Residuals vs leverage plot.
# # Aim: standardised residuals <3 absolute, or considered outliers. 
# # Aim: leverage <0.02 & no points in upper/lower R corners, or considered influential. 
# plot(lm.1.ols,5) # Indicates influential outliers. 
# 
# # Residuals and leverage are combined in the Cook's Distance metric of outlier influence.
# # An observation has high influence if Cook's Distance exceeds 
# # 4/(n - p - 1)(P. Bruce and Bruce 2017), where n is the number of observations 
# # and p the number of predictor variables. 
# cooksd <- cooks.distance(lm.1.ols)
# sample_size<- 655
# plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  
# # abline(h = 4/sample_size, col="red")  
# abline(h = 4/(sample_size -3 -1), col="blue")  
# 
# rm(cooksd, sample_size)
# 
# ## Residuals v Dependent Variable
# # Also checks for outliers.
# # Standardized residuals with absolute values >3 are outliers.
# par(mfrow=c(1,1))
# plot(mydata$DS,
#      rstandard(lm.1.ols),
#      ylab='Standardized Residuals', xlab='Cognitive Discrepancy Score')
# 
# ### Summary
# # All plots show influential outliers, violating the assumptions of OLS
# # regression and necessitating robust regression in all further models. 
# # Robust regression uses an M-estimator (maximum-likelihood) and down weights outliers.
# 
# # There are multiple functions for robust regression: 
# # lmrob(), robustbase package; lmRob(), robust package; rlm(), MASS package.
# # lmrob() was chosen for this analysis because it appears to have the most
# # recent references and updates (eg Koller and Stahel, 2011 & 2012). 
# # lmrob() is described in Yohai (1987) and Koller and Stahel (2011).
# 
# #---------------------------------------------------------------------------------------------------#
# ### Model 1 with robust regression
# 
# ## Fit robust regression
# lm.1.rob <- rlm(DSz ~ LVF1z +  LVF2z +  LVF3z + Sex + Education,
#                   data=mydata)
#
# ## Compare fit of OLS vs robust regression
# # Residual standard error (RSE) measures the SD of a model's residuals. Lower 
# # RSE indicates better model fit. 
# summary(lm.1.ols)$sigma
# summary(lm.1.rob)$sigma 
# # lmrob provides a marginally lower RSE indicating better fit.
# 
# # Do not check assumptions of lmrob(); it runs despite violations such as outliers.
# 
# ## Remove models 
# rm(lm.1.ols, lm.1.rob)


#------------------------------------------------------------------------------#
### Model 1

lm.1 <- rlm(DSz ~ LVF1z +  LVF2z +  LVF3z +
                Sex + Education,
              data=mydata)

#------------------------------------------------------------------------------#
### Model 2


lm.2 <- rlm(DSz ~ LVF1z* scale(poly(Agez,2)) + 
                 LVF2z* scale(poly(Agez,2)) +  
                 LVF3z* scale(poly(Agez,2))  + 
                 Sex + Education,
               data=mydata )

#---------------------------------------------------------------------------------------------------#
### Model 3

lm.3 <- rlm(DSz ~ 
                 # LVF1 * each medication * age squared
                 LVF1z*AntiH_Binary*scale(poly(Agez,2)) + 
                 LVF1z*drugclass2BB*  scale(poly(Agez,2)) + 
                 LVF1z*drugclass2Diu* scale(poly(Agez,2)) + 
                 LVF1z*drugclass2Stat*  scale(poly(Agez,2)) +
                 # LVF2 * each medication * age squared
                 LVF2z*AntiH_Binary* scale(poly(Agez,2)) + 
                 LVF2z*drugclass2BB*  scale(poly(Agez,2)) + 
                 LVF2z*drugclass2Diu*  scale(poly(Agez,2)) + 
                 LVF2z*drugclass2Stat*  scale(poly(Agez,2)) + 
                 # LVF * each medication * age squared
                 LVF3z*AntiH_Binary*  scale(poly(Agez,2)) + 
                 LVF3z*drugclass2BB*  scale(poly(Agez,2)) + 
                 LVF3z*drugclass2Diu*  scale(poly(Agez,2)) + 
                 LVF3z*drugclass2Stat*  scale(poly(Agez,2)) + 
                 Sex + Education,
               data=mydata) 

#---------------------------------------------------------------------------------------------------#
### Model 4

lm.4 <- rlm (DSz ~ LVF1z*   scale(poly(Agez,2)) * Sex  + 
                LVF2z*   scale(poly(Agez,2)) *Sex  +   
                LVF3z*   scale(poly(Agez,2)) *Sex + 
                Sex + Education, 
              data=mydata)

#---------------------------------------------------------------------------------------------------#
### Model 5

# In model 5, al terms interact. 
# Do not include sex interactions, if not significant in previous step.

lm.5 <- rlm (DSz ~ LVF1z * scale(poly(Agez,2))  + 
                LVF2z * scale(poly(Agez,2))  + 
                LVF3z * scale(poly(Agez,2)) +
                LVF1z * LVF2z * LVF3z +
                LVF1z * LVF2z * scale(poly(Agez,2))  +
                LVF1z * LVF3z * scale(poly(Agez,2))  +
                LVF2z * LVF3z * scale(poly(Agez,2))  +
                LVF1z * LVF2z * LVF3z * scale(poly(Agez,2))  + 
                Sex + Education,
              data=mydata)

#-----------------------------------------------------------------------------#
### Model Results

# Report uncorrected p-values in sup mat tables.

#When using scale(poly(Age,2)), the 'tab_model' function will not specify
# which age term is linear and which is quadratic. 
# However, this can be clarified using 'tab_model' in combination with
# 'summary'.

## Results of Model 1
tab_model(lm.1, show.se=TRUE, digits=2)
tab_model(lm.1, show.se=TRUE, digits=2, file="Model1_results.doc") 

tab_model(lmr.3, show.se=TRUE, digits=3)


## Results Model 2 
# summary(lm.2)
tab_model(lm.2, show.se=TRUE, digits=2)
tab_model(lm.2, show.se=TRUE, digits=2, file="Model2_results.doc") 

## Results Model 3
# summary(lm.3)
tab_model(lm.3, show.se=TRUE, digits=2)
tab_model(lm.3, show.se=TRUE, digits=2, file="Model3_results.doc") 

## Results Model 4 
# summary(lm.4)
tab_model(lm.4, show.se=TRUE, digits=2)
tab_model(lm.4, show.se=TRUE, digits=2, file="Model4_results.doc") 

## Results Model 5 
# summary(lm.5)
tab_model(lm.5, show.se=TRUE, digits=2, file="Model5_results.doc") 

#-----------------------------------------------------------------------------#
### Comparing Models 

## Extract AIC fit metrics

aic <- AIC(lm.1,lm.2,lm.3, lm.4, lm.5)
aic

## Save AIC
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")
saveRDS(aic, file="7_output_AIC.rds")

## Extract BIC fit metrics 

bic <- BIC(lm.1,lm.2,lm.3, lm.4, lm.5)
bic 


## Save BIC
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")
saveRDS(bic, file="7_output_BIC.rds")


### Model comparisons with SoS
# Running anova() on rlm() regression models, gives sum of squares. 
anova(lm.1, lm.2)
anova(lm.2, lm.3)
anova(lm.2, lm.4)
anova(lm.2, lm.5)

### Sup Mat Table 4. Model Comparisons. 

## Create table to store results in
comparisons <- as.data.frame(matrix (nrow=4, ncol=3))
rownames(comparisons)<- c("Model.1v2",
                          "Model.2v3",
                          "Model.2v4",
                          "Model.2v5")
colnames(comparisons) <- c("Difference in AIC",
                           "Difference in BIC",
                           "Difference in Sum of Squares")

### Calculations to report:

# Model 1 vs 2
comparisons$`Difference in AIC`[1] <- round((aic$AIC[1]-aic$AIC[2]),2)
comparisons$`Difference in BIC`[1] <- round((bic$BIC[1]-bic$BIC[2]),2)
temp <- round((anova(lm.1, lm.2)),2)
comparisons$`Difference in Sum of Squares`[1] <- temp$`Sum of Sq`[2]

# Model 2 vs 3
comparisons$`Difference in AIC`[2] <-round((aic$AIC[2]-aic$AIC[3]),2)
comparisons$`Difference in BIC`[2] <- round((bic$BIC[2]-bic$BIC[3]),2)
temp <- round((anova(lm.2, lm.3)),2)
comparisons$`Difference in Sum of Squares`[2] <- temp$`Sum of Sq`[2]

# Model 2 vs 4
comparisons$`Difference in AIC`[3] <-round((aic$AIC[2]-aic$AIC[4]),2)
comparisons$`Difference in BIC`[3] <- round((bic$BIC[2]-bic$BIC[4]),2)
temp <-round((anova(lm.2, lm.4)),2)
comparisons$`Difference in Sum of Squares`[3] <- temp$`Sum of Sq`[2]


# Model 2 vs 5
comparisons$`Difference in AIC`[4] <-round((aic$AIC[2]-aic$AIC[5]),2)
comparisons$`Difference in BIC`[4] <-round((bic$BIC[2]-bic$BIC[5]),2)
temp <-round((anova(lm.2, lm.5)),2)
comparisons$`Difference in Sum of Squares`[4] <- temp$`Sum of Sq`[2]

rm(temp)

# Save model comparisons
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")
saveRDS(comparisons, file="7_output_comparisons.rds")

#-----------------------------------------------------------------------------#
### Inspect winning model 2

tab_model(lm.2)

# DS is significantly predicted by LVF2*Age^2

## Visualise DS~LVF2*Age^2 

plotLVF2 <- ggplot(mydata,aes(y=DSz ,x=LVF2z, color=factor(Age_Group))) +
  geom_point(alpha=0.5)+
  stat_smooth(method="lm",se=FALSE) +
  xlab("LVF2 Factor Scores") +  ylab("Ability Discrepancy Score") +
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
  scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
  guides(col=guide_legend("Age Group")) + 
  scale_colour_manual(values = c("darkorange4", "darkorange3", "darkorange"),
                      name="Age Group",
                      breaks=c("Old", "Middle", "Young"),
                      labels=c("Old", "Middle", "Young"))
plotLVF2

### Save plot
setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Outputs/Plots/")
ggsave(plotLVF2,filename = "LVF2-DS.png")

rm(plotLVF2)

### Add age-group specific correlations to the plot.
# Subset data into age groups, then correlate VFS-DS
# Report r-values to show effect size, 
# don't report p-values because categorisation into ages is arbitrary.

# Young:
youngdata <- subset(mydata, Age_Group == 'Young',select=c(LVF2z, DSz, Age_Group))
cor.test(youngdata$DSz,youngdata$LVF2z)  # r=0.006820874 
count(mydata$Age_Group=="Young") # 154

# Middle:
middata <- subset(mydata, Age_Group == 'Middle',select=c(LVF2z, DSz, Age_Group))
cor.test(middata$DSz,middata$LVF2z)  # r=0.0172425 
count(mydata$Age_Group=="Middle") # 307

# Old:
olddata <- subset(mydata, Age_Group == 'Old',select=c(LVF2z, DSz, Age_Group))
cor.test(olddata$DSz,olddata$LVF2z)  # r=0.1335735 
count(mydata$Age_Group=="Old") # 194

rm(youngdata, middata, olddata)
rm(lm.2)

# This plot does not account for covariates or quadratic interaction, which we
# could implement using for eg sliding window visualisations, however the plot
# as it stand is adequate for reporting and ease of audience understanding.

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

rm(mydata)

### ENDS

