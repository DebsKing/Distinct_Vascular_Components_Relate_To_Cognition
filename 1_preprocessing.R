### 1. Preprocess data for analysis on vascular and cognitive factors

library(plyr) # count
library(psych) # describe

#-----------------------------------------------------------------------------#
### Preprocessing

# Raw data is imported as csv in this first script, 
# but it subsequently saved as R files to maintain
# precision (saving as csv repeatedly has a bug that looses dp).

## Load main data
# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Data/ExtractedData")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Data")
mydata <- read.csv("mydata.csv")

## Load medications data
# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Data/RawData")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Data")
drugdata <- read.csv("2021-08-09_CC700_medsBP_kat.csv")
# drugdata <- read.csv("drugdata.csv")

# ### Check datasets line up exactly
count(mydata$SubCCIDc == drugdata$SubCCIDc)
count(mydata$Age_Cognition == drugdata$Age)

#-----------------------------------------------------------------------------#
### Preprocessing medications data 

### Combine datasets, extracting drug data in numerical form:

# Four anti-hypertensives:
mydata$drugclass2ACE <- as.numeric(drugdata$drugclass2ACE)
mydata$drugclass2ARB  <- as.numeric(drugdata$drugclass2ARB)
mydata$drugclass2CCC  <- as.numeric(drugdata$drugclass2CCC)
mydata$drugclass2DiuThi  <- as.numeric(drugdata$drugclass2DiuThi)

# Other medications:
mydata$drugclass2BB  <- as.numeric(drugdata$drugclass2BB)
mydata$drugclass2Diu  <- as.numeric(drugdata$drugclass2Diu)
mydata$drugclass2Stat  <- as.numeric(drugdata$drugclass2Stat)

# Remove the medications data
rm(drugdata)

### Visualise data 
hist(mydata$drugclass2ACE)
qqnorm(mydata$drugclass2ACE)
describe(mydata$drugclass2ACE)
count(is.nan(mydata$drugclass2Stat)) # false, no NAs
count(mydata$drugclass2ACE == "NaN") # false, no NaNs

# Condense 4 antihypertensive medications into one summary var
mydata$AntiH_Total <- mydata$drugclass2ACE + mydata$drugclass2ARB + 
  mydata$drugclass2CCC + mydata$drugclass2DiuThi

plot(mydata$AntiH_Total)
range(mydata$AntiH_Total) # 0 to 2.

# Turn summary variable into binary for "on" / "off" medication status
mydata$AntiH_Binary [ mydata$AntiH_Total >= 1 ] <- 1
mydata$AntiH_Binary [ mydata$AntiH_Total < 1 ] <- 0

# Medications variables to use in subsequent analyses: 
# mydata$AntiH_Binary; mydata$drugclass2BB; mydata$drugclass2Diu; mydata$drugclass2Stat.

#-----------------------------------------------------------------------------#
### Preprocessing other variables

# Rename more usefully
mydata$Age <- mydata$Age_Cognition 
mydata$Sex <- mydata$GenderCodec_Cognition

# head(mydata$Genderc_Cognition) 
# head( mydata$GenderCodec_Cognition)
## Male = 1. Female = 2. 


### Split age into three groups, to later visualise age interactions.
# Young is 18-37. Mid is 38-67. Old is 69-88.

mydata$Age_Group [mydata$Age_Cognition<38] <- "Young"
count(mydata$Age_Group=="Young") # 164

mydata$Age_Group [37<mydata$Age_Cognition & mydata$Age_Cognition<68] <- "Middle"
count(mydata$Age_Group=="Middle") # 325

mydata$Age_Group [mydata$Age_Cognition>67] <- "Old"
count(mydata$Age_Group=="Old") # 219


### There is a data entry error for BP readings of participant participant SubCCIDc=CC510511.  
# Recode erroneous data to NAN. 
mydata$bp_dia1W1 [mydata$SubCCIDc=="CC510511" ]<- NaN
mydata$bp_dia2W1 [mydata$SubCCIDc=="CC510511" ]<- NaN
mydata$bp_dia3W1 [mydata$SubCCIDc=="CC510511" ]<- NaN
mydata$bp_sys1W1 [mydata$SubCCIDc=="CC510511" ]<- NaN
mydata$bp_sys2W1 [mydata$SubCCIDc=="CC510511" ]<- NaN
mydata$bp_sys3W1 [mydata$SubCCIDc=="CC510511" ]<- NaN


### Recode data on participants with flagged HRV data
mydata$ecgwin1_lf [mydata$flag_ecg1=="0"] <- NaN # Recode flag_ecg1=="0" to NaN 
mydata$ecgwin1_hf [mydata$flag_ecg1=="0"] <- NaN

#------------------------------------------------------------------------------#
### Save data

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

# Save as RDS for use in future scripts. 
saveRDS(mydata, file="1_mydata_preprocessed.rds")

# Save as CSV for viewing externally.
write.csv(mydata, "1_mydata_preprocessed.csv")

# # Check if CSV looses dp precision, which happens sometimes. 
# mydata <- read.csv("1_mydata_preprocessed.csv")
# # Check for decimal places in some vars:
# head(mydata$bmiW1)
# head(mydata$ecgwin1_hf)

rm(mydata)

#------------------------------------------------------------------------------#