### 3. Vascular Confirmatory Factor Analysis (CFA)

#-----------------------------------------------------------------------------#
### Load

library(GGally) # correlation plot
library(lavaan)
library(psych)
library(semPlot) 

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis/")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

mydata <- readRDS("2_mydata_preprocessed_vascular_n668.rds")

#------------------------------------------------------------------------------#
### Correlation plot of Diastolic, Systolic and Age

# Variables are in raw units, of mmHg

par(mfrow=c(1,1))
corr_pca <- data.frame( mydata$bp_dia1W1, mydata$bp_dia2W1, mydata$bp_dia3W1,
                        mydata$bp_sys1W1, mydata$bp_sys2W1, mydata$bp_sys3W1,
                        mydata$Age)

# Rename cols
colnames(corr_pca)<-c( "Diastolic1", "Diastolic2", "Diastolic3",
                       "Systolic1", "Systolic2", "Systolic3",
                       "Age")

# Detailed corrplot:
plotcor<-ggpairs(corr_pca)
plotcor

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
ggpairs( corr_pca, lower = list (continuous =  my_fn  ))


## Save plot
# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Outputs/Plots/R")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Outputs/Plots")
ggsave(plotcor,filename = "Corrplot-Diastolic-Systolic-Age.png")

## Count how many per variable, for figure legend
sum(mydata$bp_dia1W1 != "NaN") # 578
sum(mydata$bp_dia2W1 != "NaN") # 579
sum(mydata$bp_dia3W1 != "NaN") # 577

sum(mydata$bp_sys1W1 != "NaN") # 578
sum(mydata$bp_sys2W1 != "NaN") # 579
sum(mydata$bp_sys3W1 != "NaN") # 577

# #-----------------------------------------------------------------------------#
### Prepare blood pressure data for CFA

# Calculate total blood pressure
mydata$BP1 <- as.numeric(unlist(data.frame( mydata$bp_sys1W1 + mydata$bp_dia1W1 )))
mydata$BP2 <- as.numeric(unlist(data.frame( mydata$bp_sys2W1 + mydata$bp_dia2W1 )))
mydata$BP3 <- as.numeric(unlist(data.frame( mydata$bp_sys3W1 + mydata$bp_dia3W1 )))

# Calculate pulse pressure
mydata$PP1 <- as.numeric(unlist( data.frame( mydata$bp_sys1W1 - mydata$bp_dia1W1 )))
mydata$PP2 <- as.numeric(unlist( data.frame( mydata$bp_sys2W1 - mydata$bp_dia2W1 )))
mydata$PP3 <- as.numeric(unlist( data.frame( mydata$bp_sys3W1 - mydata$bp_dia3W1 )))

# Visualise the new variables
par(mfrow=c(1,1))
hist(mydata$BP1) 
qqnorm(mydata$BP1)

# Check for skew
describe(mydata[,cbind("BP1", "BP2", "BP3", "PP1", "PP2", "PP3")])
# All positively skewed so needs log transformation
# Range doesn't cover zero; suitable for input to log.

## Normalise with log:
mydata$BP1 <- log(mydata$BP1)
mydata$BP2 <- log(mydata$BP2)
mydata$BP3 <- log(mydata$BP3)

mydata$PP1 <- log(mydata$PP1)
mydata$PP2 <- log(mydata$PP2)  
mydata$PP3 <- log(mydata$PP3)

## Check log normalisation
hist(mydata$PP1)
hist(mydata$PP1, breaks=100) 
describe(mydata[,cbind("BP1", "BP2", "BP3", "PP1", "PP2", "PP3")])
# Skew and kurtosis has been reduced.

#-----------------------------------------------------------------------------#
### Prepare heart rate data for CFA

# Visualise the 3x observations
par(mfrow=c(1,3))
# Check for normality
qqnorm(mydata$pulse1W1)
qqnorm(mydata$pulse2W1)
qqnorm(mydata$pulse3W1)

# Check for skew
hist(mydata$pulse1W1) 
hist(mydata$pulse2W1) 
hist(mydata$pulse3W1) 

describe(mydata[,cbind("pulse1W1", "pulse2W1", "pulse3W1")])
# Positive skew & kurtosis
# Range does not span zero; suitable for log.

# Log transformation to normalise skew
mydata$HR1 <- log(mydata$pulse1W1)
mydata$HR2 <- log(mydata$pulse2W1)
mydata$HR3 <- log(mydata$pulse3W1)

## Check normalisation:
qqnorm(mydata$HR1)
hist(mydata$HR1, breaks=100) 
describe(mydata[,cbind( "HR1", "HR2", "HR3")])
# skew and kurtosis much reduced after log.

#-----------------------------------------------------------------------------#
### Confirmatory Factor Analysis 

# CFA on the 3x values of blood pressure.
model_BP <- 'BP =~ BP1 + BP2 + BP3'
fit_BP <- cfa(model_BP, data=mydata, 
              std.ov=T,  # std.ov - observed variables are standardized
              std.lv=T,  # std.lv - variances of LVs fixed to 1
              missing='fiml', 
              estimator='mlr')  # 'ml' is default in Lavaan for continuous data
                                # Borgeest & Fuhrmann scripts use both 'mlr' and 'fiml'. 
                                # estimator='MLR' is for both complete and incomplete data. 

# CFA on the 3x values of pule pressure.
model_PP <- 'PP =~ PP1 + PP2 + PP3'
fit_PP <- cfa(model_PP, data=mydata, std.lv=T,std.ov=T,missing='fiml', estimator='mlr') 

# Run CFA on the 3x values of HR
model_hr <- 'HR =~ HR1 + HR2 + HR3'
fit_HR <- cfa(model_hr, data=mydata,
              std.lv=T, std.ov=T, missing='fiml', estimator='mlr') 

## Draw out the models:
# par(mfrow=c(1,3))
# semPaths(fit_HR, # change to each fit
#          "par",
#          sizeMan = 20, sizeInt = 20, sizeLat=20, fade=FALSE,
#          edge.label.cex=3,# Make edge labels readable
#          residuals = FALSE,  # Remove the residual variance, within each variable
#          intercepts = FALSE,# Remove intercepts
#          rotation = 3, nCharNodes=9)

## Inspect models
# summary(fit_HR)
# summary(fit_BP, rsquare=TRUE) 

# Note that the model fit indices are unhelpful here because there are no constraints
# imposed on the models that are exactly identified (degrees of freedom=0). Such
# constraints would otherwise be used to assess model fit. Here, we do not need 
# to investigate fit because we were using CFA to essentially take a robust average. 

# Extract values for latent vars created in CFA:
BP <- predict(fit_BP) 
PP <- predict(fit_PP) 
HR <- predict(fit_HR) 

# Put latent variables into main dataframe:
mydata <- data.frame (mydata[,],BP, PP ,HR)

### Recode CFA output from 'NA' to 'NaN' for consistency
mydata$BP[is.na(mydata$BP)]<-NaN
mydata$PP[is.na(mydata$PP)]<-NaN
mydata$HR[is.na(mydata$HR)]<-NaN

# Tidy workspace
rm(model_BP, model_PP, model_hr, 
   fit_BP, fit_PP, fit_HR,
   BP, PP, HR )

#-----------------------------------------------------------------------------#
### Save data

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

# Save as RDS for use in future scripts. 
saveRDS(mydata, file="3_mydata_cfa_vascular_output_n668.rds")

# Save as CSV for viewing externally.
write.csv(mydata, "3_mydata_cfa_vascular_output_n668.csv")

# # # Check if CSV looses dp precision, which happens sometimes. 
# mydata <- read.csv("3_mydata_cfa_vascular_output_n668.csv")
# # # Check for decimal places in some vars:
# head(mydata$bmiW1)
# head(mydata$ecgwin1_hf)


rm(mydata) 

# rm(corr_pca, plotcor)
#-----------------------------------------------------------------------------#

