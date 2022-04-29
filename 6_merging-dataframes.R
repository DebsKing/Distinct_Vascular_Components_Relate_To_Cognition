#### 6. Merging dataframes from vascular (n=668) and cognitive (n=678) factor analyses

#-----------------------------------------------------------------------------#
### Load vascular and cognitive dataframes

library(plyr)

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis/")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

mydata_vas <- readRDS("4_mydata_efa_vascular_output_n668.rds")
mydata_cog <- readRDS("5_mydata_discrepancy-score_output_n678.rds")

#------------------------------------------------------------------------------#
### Prepare to merge dataframes

# Trim vascular dataframe to the variables unique to that dataframe 
# (ie the newly created LVFs), and the participant IDs
mydata_vas <- subset (mydata_vas, select=c("LVF1", "LVF2", "LVF3", "SubCCIDc"))

head(mydata_vas)

# Dataframes will next be merged based on the ID variable "SubCCIDc".
# The newly created dataframe will only contain the "SubCCIDc"variable from 
# the first of the two input dataframes.

# To allow checks that participants are lined up exactly in the merged dataframe, 
# duplicate and rename "SubCCIDc" to "SubCCIDc_cog", in the cognitive dataframe.
# This allows the merged dataframe to retain both "SubCCIDc" and SubCCIDc_cog" 
# which can then be comapred.

mydata_cog$SubCCIDc_c<- mydata_cog$SubCCIDc
tail(colnames(mydata_cog))

#------------------------------------------------------------------------------#
### Merge dataframes 

# merge  based on the column they have in common, in this csae "SubCCIDc".
# Will only retain one copy of "SubCCIDc" column. 

mydata_combined <- merge(mydata_vas[,], mydata_cog[,])

# check initial columns are those from the trimmed down 'mydata_vas' dataframe.
head(colnames(mydata_combined)) # yes

# Check the final columns are recently created DS and SubCCIDc_c from 'mydata_cog' dataframe
tail(colnames(mydata_combined)) # yes

rm(mydata_vas, mydata_cog)

#-------------------------------------------------------------------------------#
### Double check merged dataframe

## Do the IDs match for all participants
sum(mydata_combined$SubCCIDc == mydata_combined$SubCCIDc_c) # yes 
count(mydata_combined$SubCCIDc == mydata_combined$SubCCIDc_c) # yes 


### Extra sanity checks

# Are any participants remaining who were coded for removal in script 2, 
# due to having <2 observations per domain?
count(mydata_combined$Remove_cog) # no
count(mydata_combined$Remove_vas) # no

# Check that the range of observations is >=2 threshold for inclusion.
range(mydata_combined$Observations_Vas) #2-12
range(mydata_combined$Observations_Cog) #2-6

#-------------------------------------------------------------------------------#
### Save  

# setwd("C:/Users/debs_/OneDrive - University of Cambridge/1_Vascular_Cognition/Analysis/")
setwd("/imaging/camcan/sandbox/do01/1_Vascular_Cognition/Analysis")

# Save as RDS for use in future scripts. 
saveRDS(mydata_combined, file="6_mydata_for_regressions_n655.rds")

# Save as CSV for viewing externally.
write.csv(mydata_combined, "6_mydata_for_regressions_n655.csv")

rm(mydata_combined)




