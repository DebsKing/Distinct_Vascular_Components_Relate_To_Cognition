###  Multigroup CFA on vascular model

# Emulating De Mooij et al., 2018. 

# Here, I explore if models fit better when the covariances of latent vascular factors 
# are free to vary with age, vs constrained. 

#------------------------------------------------------------------------------#
### Load packages

library(lavaan) # for SEM
library(semPlot) # for plotting SEM diagram. 
# install.packages("qpcR") # install.packages('qpcR)
library(rgl) # needed for 'qpcR'
library(qpcR) # for akaike weights  

#------------------------------------------------------------------------------#
### Load data 

setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

mydata <- readRDS("3_mydata_cfa_vascular_output_n668.rds")

#------------------------------------------------------------------------------#

### Normalise vascular vars with log:
mydata$HRV_lf <- log(mydata$ecgwin1_lf)
mydata$HRV_hf <- log(mydata$ecgwin1_hf)
mydata$BMI <- log(mydata$bmiW1)

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

#------------------------------------------------------------------------------#
## Fit the CFA

# do not stipulate covariances as in EFA

model_cfa <- 'LV1 =~ 0.98*BP + 0.30*BMI
              LV2 =~ 0.98*PP + -0.37*HR
              LV3 =~ 0.85*HRV_hf + 0.82*HRV_lf
            '

# ## Check that the model looks as expected: 
# fit_cfa <- cfa(model_cfa,
#                data=mydata,
#                estimator='MLR', missing='fiml')
# # Warning: error on ov negative variances; we want HR to load negatively onto the latent variable.
# # semPaths(fit_cfa)
# # summary(fit_cfa, fit.measures=TRUE)

#------------------------------------------------------------------------------#
### Multi-group CFA, fit models

# All paths equality constrained, apart from lv covariances, which are free to vary w age
fit_free <- cfa(model_cfa,
                data=mydata,
                estimator='MLR', missing='fiml',
                group="Age_Group",
                group.equal = c("intercepts", "residuals", "lv.variances",
                                "residual.covariances",  "loadings"))

# All paths equality constrained
fit_constrained <- cfa(model_cfa,
                 data=mydata,
                 estimator='MLR', missing='fiml',
                 group="Age_Group",
                 group.equal = c("intercepts", "residuals", "lv.variances",
                                 "residual.covariances", "lv.covariances", "loadings"))

#------------------------------------------------------------------------------#
### Multi-group CFA, compare models

anova(fit_free,fit_constrained)
# because the pval is >0.05, we say models are not significantly different.

#get 2dp AIC values 
AIC(fit_free)
AIC(fit_constrained)


# Calculate Akaike weights (ie adjusted AIC values)
aics <- akaike.weights(c(AIC(fit_free), 
                         AIC(fit_constrained)))$weights # weights akiake
aics # report in table

aics[[1]]/aics[[2]] # Estimate of total preference

# Inpust sample adjusted BIC, for each model: 
bics <- akaike.weights(c(9026.267, # saBIC for the free to vary model
                         9019.033) # saBIC for the constrained model
                       )$weights  

bics # Wi(BIC), report in table 

bics[[1]]/bics[[2]] # Estimate of total preference

#------------------------------------------------------------------------------#
# ## Plot winning model, put into supplementary  material
# 
# # As both models are similar, this is not particularly interesting.
# par(mfrow=c(1,3))
# semPaths(fit_1, 
#          what='est', 
#          rotation = 2, # default rotation = 1 with four options
#          curve = 2, # pull covariances' curves out a little
#          nCharNodes = 4,
#          nCharEdges = 5, # don't limit variable name lengths
#          sizeMan = 20, # font size of manifest variable names
#          style = "lisrel", # single-headed arrows vs. # "ram"'s 2-headed for variances
#          edge.label.cex=2, 
#          curvePivot = TRUE, 
#          fade=FALSE, 
#          residuals = FALSE, intercepts = FALSE,
#          sizeInt = 10,sizeLat=10)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


