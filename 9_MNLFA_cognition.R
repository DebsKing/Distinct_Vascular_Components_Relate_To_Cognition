## -------------------------------
# Moderated non-linear factor analysis

# Script investigates invariance with age, in the confirmatory factor analysis
# model of cognitive variables: fluid and crystallized intelligence. 

# Script follows tutorial: 

### Laura Kolbe
### Last updated: 22 February 2022
### MNLFA in OpenMx with mokken::DS14 data.
### Supplementary material



## -------------------------------
## Step 1: Install and load OpenMx
## -------------------------------

library(OpenMx)


## -----------------------------
## Step 2: Load and prepare data
## -----------------------------

library(mokken)

rm(list = ls())

## Load data

# We use the participants who have data on vascular and cognitive measures, 
# because this is most relevant for understanding relationship between
# latent vascular factors and the discrepancy score.

mydata <- readRDS("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis/3_mydata_cfa_vascular_output_n668.rds")


### Note
# I also repeated teh MNFLA on complete case data only: n=504.
# ## Repeat MNFLA, on only n=504 with complete vascular and cognitive data - same resslt, of invariance.
# 
# setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")
# mydata <- readRDS("1_mydata_preprocessed.rds")
# 
# mydata <- subset(mydata, Cattell_SubScore1!="NaN" & Cattell_SubScore2!="NaN" &
#                    Cattell_SubScore3!="NaN" & Cattell_SubScore4!="NaN" &
#                    STW_total !="NaN" & Proverbs_Score!="NaN" &
#                    bp_dia1W1!="NaN" & bp_dia2W1!="NaN" & bp_dia3W1 !="NaN"&
#                    bp_sys1W1 != "NaN" & bp_sys2W1 != "NaN" & bp_sys3W1 != "NaN" &
#                    pulse1W1!= "NaN" & pulse2W1!= "NaN" & pulse3W1!= "NaN" &  # HR
#                    bmiW1 != "NaN" &  ecgwin1_lf != "NaN"& ecgwin1_hf != "NaN")
# 
# # n=504.

# #------------------------------------------------------------------------------#
##
DC <- mydata[,c("Age","Cattell_SubScore1","Cattell_SubScore2","Cattell_SubScore3","Cattell_SubScore4",
                "STW_total","Proverbs_Score")]


## Save as data frame
DC <- data.frame(DC)

## Standardize age
DC$Age <- (DC$Age - mean(DC$Age))/sd(DC$Age)
# mean-centering age is another option, but caused convergence difficulties with 
# this dataset... in the tutorial.

DC <-scale(DC,center = TRUE, scale = TRUE)

## Check data
head(DC)

## Create data object
mxdata1 <- mxData(observed=DC, type="raw")

## Indicate names and number of indicators 
manVars <- colnames(DC[,-c(1)]) # all col names apart form col1=Age.
nv <- length(manVars)


## ------------------------------------
## Step 3: Test full scalar invariance
## ------------------------------------

## Specify matrices for configural model
matT0 <- mxMatrix(type="Full", nrow=1, ncol=nv,
                  free=TRUE,
                  values=1,
                  name="matT0")
matB1 <- mxMatrix(type="Full", nrow=1, ncol=nv,
                  free=TRUE,
                  values=0,
                  name="matB1")

# this is where we specify the Cattell on LF 1 and STW&Proverbs on LF2.
matL0 <- mxMatrix(type="Full", nrow=nv, ncol=2,
                  free=c(rep(c(TRUE,FALSE),4), # i.e. the first 4 vars load onto LF1
                         rep(c(FALSE,TRUE),2)), # the last 2 vars load onto LF2
                  values=c(rep(c(1,0),4),
                           rep(c(0,1),2)),
                  byrow=TRUE,
                  name="matL0")
matC1 <- mxMatrix(type="Full", nrow=nv, ncol=2,
                  free=c(rep(c(TRUE,FALSE),4),rep(c(FALSE,TRUE),2)),
                  byrow=TRUE,
                  values=0,
                  name="matC1")

matE0 <- mxMatrix(type="Diag", nrow=nv, ncol=nv,
                  free=TRUE,
                  values=1,
                  name="matE0")
matD1 <- mxMatrix(type="Diag", nrow=nv, ncol=nv,
                  free=TRUE,
                  values=0,
                  name="matD1")

matP0 <- mxMatrix(type="Symm", nrow=2, ncol=2,
                  free=c(FALSE,TRUE,TRUE,FALSE), 
                  values=c(1,0,0,1),
                  name="matP0")
matH1 <- mxMatrix(type="Symm", nrow=2, ncol=2,
                  free=c(FALSE,TRUE,TRUE,FALSE), # to identify the model
                  values=0,
                  name="matH1")

matA0 <- mxMatrix(type="Full", nrow=2, ncol=1,
                  free=FALSE,
                  values=0,
                  name="matA0")
matG1 <- mxMatrix(type="Full", nrow=2, ncol=1,
                  free=FALSE, # to identify the model
                  values=0,
                  name="matG1")

matV1 <- mxMatrix(type="Full", nrow=1, ncol=1, 
                  free=FALSE, 
                  labels="data.Age", 
                  name = "Age")

matIa <- mxMatrix(type="Diag", nrow=2, ncol=2, 
                  free=FALSE,
                  values=1, 
                  name="matIa")

matIb <- mxMatrix(type="Full", nrow=2, ncol=2, 
                  free=FALSE, 
                  values=c(0,1,1,0),
                  name="matIb")

## Specify algebra for the dependent parameters
matT <- mxAlgebra(expression=matT0+matB1*Age, 
                  name="matT")
matL <- mxAlgebra(expression=matL0+matC1*Age, 
                  name="matL")
matE <- mxAlgebra(expression=matE0*exp(matD1*Age), 
                  name="matE")
matA <- mxAlgebra(expression=matA0+matG1*Age, 
                  name="matA")

## Specify algebra for covariance matrix of factors (transformed to ensure positive definite matrices)
matVar <- mxAlgebra(expression=(matP0*exp(matH1*Age)), 
                    name="matVar")
matR <- mxAlgebra(expression=(exp(2*(matP0+matH1*Age))-1)/
                    (exp(2*(matP0+matH1*Age))+1), 
                  name="matR")
matCov <- mxAlgebra(expression=(matIa*sqrt(matVar))%*%matR%*%(matIa*sqrt(matVar)), 
                    name="matCov")
matP <- mxAlgebra(expression=matIa*matVar+matIb*matCov, 
                  name="matP")

## Specify model-implied matrices
matC <- mxAlgebra(expression=matL%*%matP%*%t(matL)+matE, 
                  name="matC") 
matM <- mxAlgebra(expression=matT+t(matL%*%matA), 
                  name="matM") 

## Specify expectation and fit function
expF <- mxExpectationNormal(covariance="matC", 
                            means="matM",
                            dimnames=manVars)
fitF <- mxFitFunctionML() 

## Make mxModel object and run the model
modConfig <- mxModel(model="Configural", 
                     matT, matT0, matB1, 
                     matL, matL0, matC1,  
                     matE, matE0, matD1, 
                     matP, matP0, matH1, 
                     matA, matA0, matG1,   
                     matIa,matIb, matV1,  
                     matVar, matR, matCov, matM, matC, 
                     expF, fitF, mxdata1)
fitConfig <- mxRun(modConfig)
fitConfigCognition <- fitConfig
summary(fitConfigCognition) 

## Specify matrices scalar model
matB1 <- mxMatrix(type="Full", nrow=1, ncol=nv,
                  free=FALSE,
                  values=0,
                  name="matB1")

matC1 <- mxMatrix(type="Full", nrow=nv, ncol=2,
                  free=FALSE,
                  values=0,
                  name="matC1")

matH1 <- mxMatrix(type="Symm", nrow=2, ncol=2,
                  free=TRUE, 
                  values=0,
                  name="matH1")

matG1 <- mxMatrix(type="Full", nrow=2, ncol=1,
                  free=TRUE, 
                  values=0,
                  name="matG1")


## Make mxModel object and run the model
modScalar <- mxModel(model="Scalar", 
                     matT, matT0, matB1,
                     matL, matL0, matC1,
                     matE, matE0, matD1,
                     matP, matP0, matH1,
                     matA, matA0, matG1,
                     matIa, matIb, matV1, 
                     matVar, matR, matCov, matM, matC, 
                     expF, fitF, mxdata1)
fitScalar <- mxRun(modScalar)
summary(fitScalar)

## Compare fit of unconstrained model with constrained model
miTest <- mxCompare(fitConfig, fitScalar)
miTest[2,c(7,8,9)] # numbers to report. 
miTest$p[2] < 0.001

# p>0.05 gives no evidence for variance with age. 
# Therefore we do not need to continue with the tutorial, exploring partial invariance.

### END

# 
# 
# ## ---------------------------------
# ## Step 4: Select anchor indicators
# ## ---------------------------------
# 
# ## Run unconstrained model for each indicator
# fitAbo <- list()
# 
# for (i in 1:6){
#   freeparT <- matrix(FALSE, nrow=1, ncol=6)
#   freeparT[i] <- TRUE
#   freeparL <- matrix(FALSE, nrow=6, ncol=2)
#   freeparL[i,ifelse(i<8,1,2)] <- TRUE
#   matB1 <- mxMatrix(type="Full", nrow=1, ncol=nv,
#                     free=freeparT,
#                     values=0,
#                     name="matB1")
# 
#   matC1 <- mxMatrix(type="Full", nrow=nv, ncol=2,
#                     free=freeparL,
#                     values=0,
#                     byrow=TRUE,
#                     name="matC1")
#   
#   modAbo <- mxModel(model=paste0("All_but_", i), 
#                    matT, matT0, matB1, 
#                    matL, matL0, matC1, 
#                    matE, matE0, matD1, 
#                    matP, matP0, matH1, 
#                    matA, matA0, matG1, 
#                    matIa, matIb, matV1,
#                    matVar, matR, matCov, matM, matC, 
#                    expF, fitF, mxdata1)
#   fitAbo[[i]] <- mxRun(modAbo)
# }
# 
# ## Compare constrained model with all unconstrained models
# anchorTest <- mxCompare(fitAbo, fitScalar)
# anchorOut <- data.frame(Name=paste0("Indicator",1:6), 
#                     X2=anchorTest$diffLL[seq(2,12,2)],
#                     df=anchorTest$diffdf[seq(2,12,2)],
#                     p=anchorTest$p[seq(2,12,2)])
# anchorOut
# 
# ## Select two indicators per factor with smallest X2 as anchor
# anchorOut[order(anchorOut$X2[1:4]),]
# anchorOut[4+order(anchorOut$X2[5:6]),]
# 
# ## Save anchors in object 
# anchors1 <- c(4,2)
# anchors2 <- c(5,6)
# 
# 
# ## -------------------------------
# ## Step 5: Test partial invariance
# ## -------------------------------
# 
# ## Specify which DIF effects can be estimated
# freeparT <- matrix(data=TRUE, nrow=1, ncol=6)
# freeparT[1,c(anchors1,anchors2)] <- FALSE
# 
# freeparL <- matrix(c(rep(c(TRUE,FALSE),4),rep(c(FALSE,TRUE),2)),nv,2, byrow=TRUE)
# freeparL[anchors1,1] <- FALSE
# freeparL[anchors2,2] <- FALSE
# 
# ## Specify matrices for unconstrained model
# matB1 <- mxMatrix(type="Full", nrow=1, ncol=nv,
#                   free=freeparT,
#                   values=0,
#                   name="matB1")
# matB2 <- mxMatrix(type="Full", nrow=1, ncol=nv,
#                   free=freeparT,
#                   values=0,
#                   name="matB2")
# matC1 <- mxMatrix(type="Full", nrow=nv, ncol=2,
#                   free=freeparL,
#                   values=0,
#                   name="matC1")
# matC2 <- mxMatrix(type="Full", nrow=nv, ncol=2,
#                   free=freeparL,
#                   values=0,
#                   name="matC2")
# 
# ## Make mxModel object and run the model
# modAnchors <- mxModel(model="AnchorsOnly", 
#                       matT, matT0, matB1, matB2,
#                       matL, matL0, matC1, matC2, 
#                       matE, matE0, matD1, matD2,
#                       matP, matP0, matH1, matH2,
#                       matA, matA0, matG1, matG2,  
#                       matIa, matIb, matV1, matV2, 
#                       matVar, matR, matCov, matM, matC, 
#                       expF, fitF, mxdata1)
# fitAnchors <- mxRun(modAnchors)
# 
# ## Run constrained model for each indicator (except the anchors)
# testIn <- c(1:14)[-c(anchors1,anchors2)]
# fitApo <- list()
# 
# for (i in testIn){
#   freeparTa <- freeparT
#   freeparLa <- freeparL
#   freeparTa[1,i] <- FALSE
#   freeparLa[i,ifelse(i<8,1,2)] <- FALSE
#   matB1 <- mxMatrix(type="Full", nrow=1, ncol=nv,
#                     free=freeparTa,
#                     values=0,
#                     name="matB1")
#   matB2 <- mxMatrix(type="Full", nrow=1, ncol=nv,
#                     free=freeparTa,
#                     values=0,
#                     name="matB2")
#   matC1 <- mxMatrix(type="Full", nrow=nv, ncol=2,
#                     free=freeparLa,
#                     values=0,
#                     name="matC1")
#   matC2 <- mxMatrix(type="Full", nrow=nv, ncol=2,
#                     free=freeparLa,
#                     values=0,
#                     name="matC2")
#   modApo <- mxModel(model=paste0("Anchors_plus_", i), 
#                     matT, matT0, matB1, matB2,
#                     matL, matL0, matC1, matC2, 
#                     matE, matE0, matD1, matD2,
#                     matP, matP0, matH1, matH2,
#                     matA, matA0, matG1, matG2,  
#                     matIa, matIb, matV1, matV2, 
#                     matVar, matR, matCov, matM, matC, 
#                     expF, fitF, mxdata1)
#   fitApo[[i]] <- mxRun(modApo)
# }
# 
# ## Compare fit of unconstrained model with all constrained models
# piTest <- mxCompare(fitAnchors, fitApo)
# piOut <- data.frame(Name=paste0("Indicator",testIn),
#                     X2=piTest$diffLL[2:11],
#                     df=piTest$diffdf[2:11],
#                     p=piTest$p[2:11],
#                     p.bonferroni=p.adjust(p=piTest$p[2:11], method="bonferroni"),
#                     p.BH=p.adjust(p=piTest$p[2:11], method="BH"))
# piOut
# 
# 
# ## -----------------------
# ## Final model with OpenMx
# ## -----------------------
# 
# ## Save start time
# start_OpenMx <- Sys.time()
# 
# ## Specify all non-DIF indicators
# finalIn <- c(1,2,3,4,5,6,8,9,10,11,12,14)
# 
# ## Specify free parameters
# freeparTb <- freeparT
# freeparLb <- freeparL
# 
# for(i in finalIn){
#   freeparTb[1,i] <- FALSE
#   freeparLb[i,ifelse(i<8,1,2)] <- FALSE
# }
# 
# ## Re-specify matrices
# matB1 <- mxMatrix(type="Full", nrow=1, ncol=nv,
#                   free=freeparTb,
#                   values=0,
#                   name="matB1")
# matB2 <- mxMatrix(type="Full", nrow=1, ncol=nv,
#                   free=freeparTb,
#                   values=0,
#                   name="matB2")
# matC1 <- mxMatrix(type="Full", nrow=nv, ncol=2,
#                   free=freeparLb,
#                   values=0,
#                   name="matC1")
# matC2 <- mxMatrix(type="Full", nrow=nv, ncol=2,
#                   free=freeparLb,
#                   values=0,
#                   name="matC2")
# 
# ## Fit final model
# modOpenmxPartial <- mxModel(model="PartialInvariance", 
#                       matT, matT0, matB1, matB2,
#                       matL, matL0, matC1, matC2, 
#                       matE, matE0, matD1, matD2,
#                       matP, matP0, matH1, matH2,
#                       matA, matA0, matG1, matG2,  
#                       matIa, matIb, matV1, matV2, 
#                       matVar, matR, matCov, matM, matC, 
#                       expF, fitF, mxdata1)
# fitOpenmxPartial <- mxRun(modOpenmxPartial)
# 
# ## Save end time
# end_OpenMx <- Sys.time()
# 
# 
# ## ---------------------------
# ## Create plots of final model
# ## ---------------------------
# 
# ## Functions for expected X7 and X13 scores
# FexpX7 <- function(theta, Male, Age){ 
#   fitOpenmxPartial$matT0$values[,7] + fitOpenmxPartial$matB1$values[,7]*Male + 
#   fitOpenmxPartial$matB2$values[,7]*Age + (fitOpenmxPartial$matL0$values[7,1] + 
#   fitOpenmxPartial$matC1$values[7,1]*Male + fitOpenmxPartial$matC2$values[7,1]*
#   Age)*theta
# }
# 
# FexpX13 <- function(theta, Male, Age){ 
#   fitOpenmxPartial$matT0$values[,13] + fitOpenmxPartial$matB1$values[,13]*Male + 
#   fitOpenmxPartial$matB2$values[,13]*Age + (fitOpenmxPartial$matL0$values[13,2] + 
#   fitOpenmxPartial$matC1$values[13,2]*Male + fitOpenmxPartial$matC2$values[13,2]*
#   Age)*theta
# }
# 
# ## Plot tracelines for DIF indicators
# library(ggplot2)
# myTheme <- theme(axis.title = element_text(size = 12),
#                  axis.text.y = element_text(size = 10),
#                  axis.text.x = element_text(size = 10),
#                  strip.text = element_text(size = 10),
#                  legend.text = element_text(size = 12),
#                  legend.title = element_text(size = 12),
#                  legend.key.size = grid::unit(5, units = "char"),
#                  legend.key = element_rect(fill = "white", colour = "white"),
#                  panel.background = element_rect(fill = "white", color = "grey"),
#                  panel.grid.minor = element_line(colour = "white"))
# 
# All <- expand.grid("Theta" = rep(seq(-4, 4, 0.01), 6), "Male" = 0:1, "Age" = c(-2,2), 
#                    "expX7" = NA, "expX13" = NA)
# for (i in 1:nrow(All)){All$expX7[i] <- FexpX7(All$Theta[i], All$Male[i], All$Age[i])
#                        All$expX13[i] <- FexpX13(All$Theta[i], All$Male[i], All$Age[i])
# }
# All$Male <- factor(All$Male, labels = c("Female","Male"))
# All$Age <- factor(All$Age, labels = c("Age = 38", "Age = 80"))
# # Age is converted from standardized score to actual age score:
# # Mean age = 58.66913, SD age = 10.54776
# # -2 --> -2*10.54776+58.66913 = 37.57361
# # 2 --> 2*10.54776+58.66913 = 79.76465
# 
# ggplot(All, aes(x = Theta, y = expX7, lty = Age, colour = Male)) + 
#   geom_line(size = 1) + 
#   scale_color_manual(values=c("grey80", "black")) +
#   scale_linetype_manual(values=c(1, 5)) +
#   labs(x = "SI", y = "Expected X7")  + 
#   myTheme 
# 
# ggplot(All, aes(x = Theta, y = expX13, lty = Age, colour = Male)) + 
#   geom_line(size = 1) + 
#   scale_color_manual(values=c("grey80", "black")) +
#   scale_linetype_manual(values=c(1, 5)) +
#   labs(x = "NA", y = "Expected X13")  + 
#   myTheme 
# 
# ## Functions for moderated factor parameters
# FacMean1 <- function(Male, Age){
#   fitOpenmxPartial$matA0$values[1,1] + fitOpenmxPartial$matG1$values[1,1]*Male +
#   fitOpenmxPartial$matG2$values[1,1]*Age
# }
# 
# FacMean2 <- function(Male, Age){
#   fitOpenmxPartial$matA0$values[2,1] + fitOpenmxPartial$matG1$values[2,1]*Male +
#   fitOpenmxPartial$matG2$values[2,1]*Age 
# }
# 
# FacVar1 <- function(Male, Age){
#   fitOpenmxPartial$matP0$values[1,1] * exp(fitOpenmxPartial$matH1$values[1,1]*Male +
#   fitOpenmxPartial$matH2$values[1,1]*Age) 
# }
# 
# FacVar2 <- function(Male, Age){
#   fitOpenmxPartial$matP0$values[2,2] * exp(fitOpenmxPartial$matH1$values[2,2]*Male +
#     fitOpenmxPartial$matH2$values[2,2]*Age) 
# }
# 
# FacCorr <- function(Male, Age){
#   (exp(2*(fitOpenmxPartial$matP0$values[1,2] + fitOpenmxPartial$matH1$values[1,2]*Male +
#   fitOpenmxPartial$matH2$values[1,2]*Age)) - 1)/(exp(2*(fitOpenmxPartial$matP0$values[1,2] + 
#   fitOpenmxPartial$matH1$values[1,2]*Male + fitOpenmxPartial$matH2$values[1,2]*Age)) + 1) 
# }
# 
# ## Plots factor parameters
# All2 <- expand.grid("Male" = 0:1, "Age" = seq(-3, 3, 0.01), "A_1" = NA, 
#                     "A_2" = NA, "P_1" = NA, "P_2" = NA, "Corr" = NA)
# for (i in 1:nrow(All2)){
#   All2$A_1[i] <- FacMean1(All2$Male[i], All2$Age[i])
#   All2$A_2[i] <- FacMean2(All2$Male[i], All2$Age[i])
#   All2$P_1[i] <- FacVar1(All2$Male[i], All2$Age[i])
#   All2$P_2[i] <- FacVar2(All2$Male[i], All2$Age[i])
#   All2$Corr[i] <- FacCorr(All2$Male[i], All2$Age[i])
# }
# All2$Male <- factor(All2$Male, labels = c("Female","Male"))
# 
# ggplot(All2, aes(x = Age, y = A_1, colour = Male)) + 
#   geom_line(size = 1) + 
#   scale_color_manual(values=c("grey80", "black")) +
#   scale_linetype_manual(values=c(1, 5)) +
#   labs(x = "Age", y = "SI mean")  + 
#   scale_x_continuous(breaks=c(-2, 0, 2), labels=c("38", "59","80")) +
#   myTheme 
# 
# ggplot(All2, aes(x = Age, y = A_2, colour = Male)) + 
#   geom_line(size = 1) + 
#   scale_color_manual(values=c("grey80", "black")) +
#   scale_linetype_manual(values=c(1, 5)) +
#   labs(x = "Age", y = "NA mean")  + 
#   scale_x_continuous(breaks=c(-2, 0, 2), labels=c("38", "59","80")) +
#   myTheme 
# 
# ggplot(All2, aes(x = Age, y = P_1, colour = Male)) + 
#   geom_line(size = 1) + 
#   scale_color_manual(values=c("grey80", "black")) +
#   scale_linetype_manual(values=c(1, 5)) +
#   labs(x = "Age", y = "SI variance")  + 
#   scale_x_continuous(breaks=c(-2, 0, 2), labels=c("38", "59","80")) +
#   myTheme 
# 
# ggplot(All2, aes(x = Age, y = P_2, colour = Male)) + 
#   geom_line(size = 1) + 
#   scale_color_manual(values=c("grey80", "black")) +
#   scale_linetype_manual(values=c(1, 5)) +
#   labs(x = "Age", y = "NA variance")  + 
#   scale_x_continuous(breaks=c(-2, 0, 2), labels=c("38", "59","80")) +
#   myTheme
# 
# ggplot(All2, aes(x = Age, y = Corr, colour = Male)) + 
#   geom_line(size = 1) + 
#   scale_color_manual(values=c("grey80", "black")) +
#   scale_linetype_manual(values=c(1, 5)) +
#   labs(x = "Age", y = "Factor correlation")  + 
#   scale_x_continuous(breaks=c(-2, 0, 2), labels=c("38", "59","80")) +
#   myTheme
# 
# 
# ## ------------------------
# ## Compare results to Mplus
# ## ------------------------
# 
# library(MplusAutomation)
# 
# ## Save start time
# start_Mplus <- Sys.time()
# 
# ## Create folder
# pathfix <- "~/FinalModel"
# dir.create(pathfix)
# 
# ## Store data in folder
# prepareMplusData(df=DS14, filename=paste0(pathfix, "/DS14dat.dat"))
# 
# ## Create and run partial-invariance model
# modMplusPartial <- '
#   DATA:
#   FILE = "DS14dat.dat";
#   
#   VARIABLE:
#   NAMES = Male Age x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14; 
#   CONSTRAINT = Male Age;
#   MISSING = .;
#   
#   ANALYSIS:
#   estimator is ML;
#   TOLERANCE = 1;
# 
#   MODEL:
#   SI by x1-x7* (l1-l7);
#   NA by x8-x14* (l8-l14);
#   SI (SIvar); 
#   NA (NAvar);
#   SI with NA (cov);
#   [SI] (SImean); 
#   [NA] (NAmean);
#   [x1-x14] (t1-t14);
#   x1-x14 (e1-e14);
#   Male @1;
#   Age @1;
#   [Male @0];
#   [Age @0];
#   
#   MODEL CONSTRAINT:
#   NEW (
#   h1_SI h1_NA g1_SI g1_NA
#   h2_SI h2_NA g2_SI g2_NA
#   p0 h1_cov h2_cov 
#   t7_0 b7_1 b7_2 t13_0 b13_1 b13_2
#   l7_0 c7_1 c7_2 l13_0 c13_1 c13_2
#   e1_0 e2_0 e3_0 e4_0 e5_0 e6_0 e7_0 e8_0 e9_0 e10_0 e11_0 e12_0 e13_0 e14_0
#   d1_1 d2_1 d3_1 d4_1 d5_1 d6_1 d7_1 d8_1 d9_1 d10_1 d11_1 d12_1 d13_1 d14_1
#   d1_2 d2_2 d3_2 d4_2 d5_2 d6_2 d7_2 d8_2 d9_2 d10_2 d11_2 d12_2 d13_2 d14_2);
#   
#   SIvar = 1 * EXP(h1_SI*Male + h2_SI*Age);
#   NAvar = 1 * EXP(h1_NA*Male + h2_NA*Age);
#   SImean = 0 + g1_SI*Male + g2_SI*Age;
#   NAmean = 0 + g1_NA*Male + g2_NA*Age;
#   
#   cov = SQRT(EXP(h1_SI*Male + h2_SI*Age))*
#         SQRT(EXP(h1_NA*Male + h2_NA*Age))*
#         (EXP(2*(p0 + h1_cov*Male + h2_cov*Age))-1)/
#         (EXP(2*(p0 + h1_cov*Male + h2_cov*Age))+1);
# 
#   e1 = e1_0 * EXP(d1_1*Male + d1_2*Age);
#   e2 = e2_0 * EXP(d2_1*Male + d2_2*Age);
#   e3 = e3_0 * EXP(d3_1*Male + d3_2*Age);
#   e4 = e4_0 * EXP(d4_1*Male + d4_2*Age);
#   e5 = e5_0 * EXP(d5_1*Male + d5_2*Age);
#   e6 = e6_0 * EXP(d6_1*Male + d6_2*Age);
#   e7 = e7_0 * EXP(d7_1*Male + d7_2*Age);
#   e8 = e8_0 * EXP(d8_1*Male + d8_2*Age);
#   e9 = e9_0 * EXP(d9_1*Male + d9_2*Age);
#   e10 = e10_0 * EXP(d10_1*Male + d10_2*Age);
#   e11 = e11_0 * EXP(d11_1*Male + d11_2*Age);
#   e12 = e12_0 * EXP(d12_1*Male + d12_2*Age);
#   e13 = e13_0 * EXP(d13_1*Male + d13_2*Age);
#   e14 = e14_0 * EXP(d14_1*Male + d14_2*Age);
#   t7 = t7_0 + b7_1*Male + b7_2*Age;
#   l7 = l7_0 + c7_1*Male + c7_2*Age;
#   t13 = t13_0 + b13_1*Male + b13_2*Age;
#   l13 = l13_0 + c13_1*Male + c13_2*Age;'
# 
# cat(modMplusPartial, file = paste0(pathfix, "/modPartial.inp", sep = ""))
# 
# ## Run model and save output
# runModels(pathfix)
# fitMplusPartial <- readModels(pathfix)
# 
# ## Save end time
# end_Mplus <- Sys.time()
# 
# ## Compare parameters estimates
# summary(fitOpenmxPartial)$parameters
# fitMplusPartial$parameters
# # NOTE: in the Mplus parameter list, the parameters being a function of
# # the background variables are also listed (can be ignored) and equal 999.000 
# # because they differ per individual
# 
# ## Compare times
# end_OpenMx - start_OpenMx
# end_Mplus - start_Mplus
