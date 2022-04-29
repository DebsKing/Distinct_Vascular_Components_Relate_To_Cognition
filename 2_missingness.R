### 2. Investigating missingness and removing participants with too few recorded observations. 

#-----------------------------------------------------------------------------#

library(plyr)
library(utils)

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")
mydata <- readRDS("1_mydata_preprocessed.rds")

#------------------------------------------------------------------------------#
### Encode data as binary, based on missing or present

# 1 = present
# 0 = missing data
# Put into numeric format.
# Across all vascular & cognitive observed variables.

mydata$bp_dia1W1_mask [mydata$bp_dia1W1 != "NaN"] <- "1"
mydata$bp_dia1W1_mask [mydata$bp_dia1W1 == "NaN"] <- "0"
mydata$bp_dia1W1_mask <- as.numeric(unlist(mydata$bp_dia1W1_mask))

mydata$bp_dia2W1_mask [mydata$bp_dia2W1 != "NaN"] <- "1"
mydata$bp_dia2W1_mask [mydata$bp_dia2W1 == "NaN"] <- "0"
mydata$bp_dia2W1_mask <- as.numeric(unlist(mydata$bp_dia2W1_mask))

mydata$bp_dia3W1_mask [mydata$bp_dia3W1 != "NaN"] <- "1"
mydata$bp_dia3W1_mask [mydata$bp_dia3W1 == "NaN"] <- "0"
mydata$bp_dia3W1_mask <- as.numeric(unlist(mydata$bp_dia3W1_mask))

mydata$bp_sys1W1_mask [mydata$bp_sys1W1 != "NaN"] <- "1"
mydata$bp_sys1W1_mask [mydata$bp_sys1W1 == "NaN"] <- "0"
mydata$bp_sys1W1_mask <- as.numeric(unlist(mydata$bp_sys1W1_mask))

mydata$bp_sys2W1_mask [mydata$bp_sys2W1 != "NaN"] <- "1"
mydata$bp_sys2W1_mask [mydata$bp_sys2W1 == "NaN"] <- "0"
mydata$bp_sys2W1_mask <- as.numeric(unlist(mydata$bp_sys2W1_mask))

mydata$bp_sys3W1_mask [mydata$bp_sys3W1 != "NaN"] <- "1"
mydata$bp_sys3W1_mask [mydata$bp_sys3W1 == "NaN"] <- "0"
mydata$bp_sys3W1_mask <- as.numeric(unlist(mydata$bp_sys3W1_mask))

mydata$pulse1W1_mask [mydata$pulse1W1 != "NaN"] <- "1"
mydata$pulse1W1_mask [mydata$pulse1W1 == "NaN"] <- "0"
mydata$pulse1W1_mask <- as.numeric(unlist(mydata$pulse1W1_mask))

mydata$pulse2W1_mask [mydata$pulse2W1 != "NaN"] <- "1"
mydata$pulse2W1_mask [mydata$pulse2W1 == "NaN"] <- "0"
mydata$pulse2W1_mask <- as.numeric(unlist(mydata$pulse2W1_mask))

mydata$pulse3W1_mask [mydata$pulse3W1 != "NaN"] <- "1"
mydata$pulse3W1_mask [mydata$pulse3W1 == "NaN"] <- "0"
mydata$pulse3W1_mask <- as.numeric(unlist(mydata$pulse3W1_mask))

mydata$ecgwin1_lf_mask [mydata$ecgwin1_lf != "NaN"] <- "1"
mydata$ecgwin1_lf_mask [mydata$ecgwin1_lf == "NaN"] <- "0"
mydata$ecgwin1_lf_mask <- as.numeric(unlist(mydata$ecgwin1_lf_mask))

mydata$ecgwin1_hf_mask [mydata$ecgwin1_hf != "NaN"] <- "1"
mydata$ecgwin1_hf_mask [mydata$ecgwin1_hf == "NaN"] <- "0"
mydata$ecgwin1_hf_mask <- as.numeric(unlist(mydata$ecgwin1_hf_mask))

mydata$bmiW1_mask [mydata$bmiW1 != "NaN"] <- "1"
mydata$bmiW1_mask [mydata$bmiW1 == "NaN"] <- "0"
mydata$bmiW1_mask <- as.numeric(unlist(mydata$bmiW1_mask))

mydata$Cattell_SubScore1_mask [mydata$Cattell_SubScore1 != "NaN"] <- "1"
mydata$Cattell_SubScore1_mask [mydata$Cattell_SubScore1 == "NaN"] <- "0"
mydata$Cattell_SubScore1_mask <- as.numeric(unlist(mydata$Cattell_SubScore1_mask))

mydata$Cattell_SubScore2_mask [mydata$Cattell_SubScore2 != "NaN"] <- "1"
mydata$Cattell_SubScore2_mask [mydata$Cattell_SubScore2 == "NaN"] <- "0"
mydata$Cattell_SubScore2_mask <- as.numeric(unlist(mydata$Cattell_SubScore2_mask))

mydata$Cattell_SubScore3_mask [mydata$Cattell_SubScore3 != "NaN"] <- "1"
mydata$Cattell_SubScore3_mask [mydata$Cattell_SubScore3 == "NaN"] <- "0"
mydata$Cattell_SubScore3_mask <- as.numeric(unlist(mydata$Cattell_SubScore3_mask))

mydata$Cattell_SubScore4_mask [mydata$Cattell_SubScore4 != "NaN"] <- "1"
mydata$Cattell_SubScore4_mask [mydata$Cattell_SubScore4 == "NaN"] <- "0"
mydata$Cattell_SubScore4_mask <- as.numeric(unlist(mydata$Cattell_SubScore4_mask))

mydata$Proverbs_Score_mask [mydata$Proverbs_Score != "NaN"] <- "1"
mydata$Proverbs_Score_mask [mydata$Proverbs_Score == "NaN"] <- "0"
mydata$Proverbs_Score_mask <- as.numeric(unlist(mydata$Proverbs_Score_mask))

mydata$STW_total_mask [mydata$STW_total != "NaN"] <- "1"
mydata$STW_total_mask [mydata$STW_total == "NaN"] <- "0"
mydata$STW_total_mask <- as.numeric(unlist(mydata$STW_total_mask))

#-----------------------------------------------------------------------------#
### Quantify complete data per participant

### Loop calculates the total number of variables recorded per participant
for (i in 1:708) {
  mydata$Observations_all[i] <- (sum(
    # as.numeric(
    mydata$bp_dia1W1_mask[i], mydata$bp_dia2W1_mask[i], mydata$bp_dia3W1_mask[i],
    mydata$bp_sys1W1_mask[i],mydata$bp_sys2W1_mask[i],mydata$bp_sys3W1_mask[i],
    mydata$pulse1W1_mask[i], mydata$pulse2W1_mask[i], mydata$pulse3W1_mask[i],
    mydata$ecgwin1_lf_mask[i],  mydata$ecgwin1_hf_mask[i],
    mydata$bmiW1_mask[i],
    mydata$Cattell_SubScore1_mask[i], mydata$Cattell_SubScore2_mask[i],
    mydata$Cattell_SubScore3_mask[i], mydata$Cattell_SubScore4_mask[i],
    mydata$Proverbs_Score_mask[i], mydata$STW_total_mask[i] ))
}
range(as.numeric(mydata$Observations_all)) #0-18.

## Total VASCULAR variables recorded per participant
for (i in 1:708) {
  mydata$Observations_Vas[i] <- (sum
                                 (mydata$bp_dia1W1_mask[i], mydata$bp_dia2W1_mask[i], mydata$bp_dia3W1_mask[i],
                                   mydata$bp_sys1W1_mask[i],mydata$bp_sys2W1_mask[i],mydata$bp_sys3W1_mask[i],
                                   mydata$pulse1W1_mask[i], mydata$pulse2W1_mask[i], mydata$pulse3W1_mask[i],
                                   mydata$ecgwin1_lf_mask[i],  mydata$ecgwin1_hf_mask[i],
                                   mydata$bmiW1_mask[i]))
}
range(mydata$Observations_Vas) # 0-12
# hist(mydata$Observations_Vas, breaks=13)

## Total COGNITIVE variables recorded per participant
for (i in 1:708) {
  mydata$Observations_Cog[i] <- (sum
                                 (mydata$Cattell_SubScore1_mask[i], mydata$Cattell_SubScore2_mask[i],
                                   mydata$Cattell_SubScore3_mask[i], mydata$Cattell_SubScore4_mask[i],
                                   mydata$Proverbs_Score_mask[i], mydata$STW_total_mask[i] ))
}
rm(i)
range(mydata$Observations_Cog) # 0-6
# hist(mydata$Observations_Cog, breaks=7)

#-----------------------------------------------------------------------------#
### Summary on removing participants with low numbers of recorded variables

# Remove participants with <2/12 recorded vascular variables. 
# This gives a dataframe for vascular factor analyses, with n=668.

# Separately, remove participants with <2/6 recorded cognitive varaibles. 
# This gives a dataframe for cognitive factor analyses, with n=678. 

# Using two separate dataframes for the domain-specific factor analyses ensures
# that maximal variance is inputted to each analysis. 
# After the separate factor analyses, the two dataframes are merged based on the 
# 'subccidc' variables. This gives a final dataframe of n=655 for input to the
# linear regressions.

#-----------------------------------------------------------------------------#
### Removing participants based on vascular data

# Label participants to keep and remove  
mydata$Remove_vas [mydata$Observations_Vas >=2 ] <- "0" # keep
mydata$Remove_vas [mydata$Observations_Vas <2  ] <- "1" # remove

count(mydata$Remove_vas == 1) # Keep n=668. 

mydata_vas <-mydata 
# View(mydata_vas) # this adds 1st col of 1-708

# Remove participants coded '1' in Remove
for (i in 708:1) {                      # run 708:1, from bottom of df up, when removing rows.
  if (mydata_vas$Remove_vas[i] == 1){   # if labelled '1' for removal
    mydata_vas <- mydata_vas[-(i),]     # Remove that row
  }
}
rm(i)
dim(mydata_vas) # n=668xcols  

#-----------------------------------------------------------------------------#
### Save vascular

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

# Save as RDS for use in future scripts. 
saveRDS(mydata_vas, file="2_mydata_preprocessed_vascular_n668.rds")

# Save as CSV for viewing externally.
write.csv(mydata_vas, "2_mydata_preprocessed_vascular_n668.csv")

# # # Check if CSV looses dp precision, which happens sometimes. 
# mydata_vas <- read.csv("1_mydata_preprocessed.csv")
# # Check for decimal places in some vars:
# head(mydata$bmiW1)
# head(mydata$ecgwin1_hf)

rm(mydata_vas)
#-----------------------------------------------------------------------------#
### Removing participants based on cognitive data

# Label participants to keep and remove   
mydata$Remove_cog [mydata$Observations_Cog >=2 ] <- "0" # keep
mydata$Remove_cog [mydata$Observations_Cog <2  ] <- "1" # remove

count(mydata$Remove_cog == 1) # Keep n=678. 

mydata_cog <-mydata 

# View(mydata_cog) # this adds a 1st col of 1-708

## Remove participants coded '1' in Remove
for (i in 708:1) {                        # run from bottom of df up, when removing rows.
  if (mydata_cog$Remove_cog[i] == 1){     # remove iterations of '1'
    mydata_cog <- mydata_cog[-(i),]       # Remove that row
  }
}
rm(i)

dim(mydata_cog) # n=678xcols 

#-----------------------------------------------------------------------------#
### Save cognitive

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

# Save as RDS for use in future scripts. 
saveRDS(mydata_cog, file="2_mydata_preprocessed_cognitive_n678.rds")

# Save as CSV for viewing externally.
write.csv(mydata_cog, "2_mydata_preprocessed_cognitive_n678.csv")

# # Check if CSV looses dp precision, which happens sometimes.
# mydata_cog <- read.csv("2_mydata_preprocessed_cognitive_n678.csv")
# # Check for decimal places in some vars:
# head(mydata$bmiW1)
# head(mydata$ecgwin1_hf)

rm(mydata_cog)

#-----------------------------------------------------------------------------#

rm(mydata)
