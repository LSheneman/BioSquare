library(ggplot2)
library(data.table)
library(scales)
library(plyr)
library(reshape2)
library(Rmisc)

setwd("/mnt/home/chjack/BioSquareGT/BRH/Data/SetL") # change

PopFiles = list.files(pattern = glob2rx("POP*.txt"))

#  =========== For the extinction graphs =====================

# For each run, it takes the number of rows where each column is 0, which indicates it is no longer present in the environment (extinct). It then binds with the run information tossing out the individual run data (which is needed for later).
allData <- lapply(PopFiles, function (x){
  dat <- read.csv(x, header = FALSE, col.names = c("Generation", "Type0", "Type1", "Type2", "Edit_PH", "Edit_PM", "Edit_MH", "Euclid_PH", "Euclid_PM","Euclid_MH", "Var_Plant", "Var_Herb","Var_Micro"))
  parts = strsplit(x, "_")[[1]]
  dat$run = as.factor(parts[3])
  newRun <- data.frame(cbind(x, Plant = length(dat[dat[,2]==0,][,1]), Herb =length(dat[dat[,3]==0,][,1]), Microbe = length(dat[dat[,4]==0,][,1]), Px = 0.01, Py = 0.99, Hx = 0.99, Hy = 0.01, Mx = 0.9, My = 0.1, MicrobeDist = "Far", Hypothesis = "BRH"))
 
  return(newRun)
})
BigDataFrame <- rbindlist(allData)


# This allows conversion of variables to numbers
BigDataFrame$Plant <- as.numeric(as.character(BigDataFrame$Plant))
BigDataFrame$Herb <- as.numeric(as.character(BigDataFrame$Herb))
BigDataFrame$Microbe <- as.numeric(as.character(BigDataFrame$Microbe))

BigDataFrame$Px <- as.numeric(as.character(BigDataFrame$Px))
BigDataFrame$Py <- as.numeric(as.character(BigDataFrame$Py))
BigDataFrame$Hx <- as.numeric(as.character(BigDataFrame$Hx))
BigDataFrame$Hy <- as.numeric(as.character(BigDataFrame$Hy))
BigDataFrame$Mx <- as.numeric(as.character(BigDataFrame$Mx))
BigDataFrame$My <- as.numeric(as.character(BigDataFrame$My))

# This translates the row extinct numbers to 0 if row number is greater than 0 or 1 if not. Later, the sums of the 1s will give us the probability of survival.
BigDataFrame$Plant2 <- ifelse(BigDataFrame$Plant >0, 0,1)
BigDataFrame$Herb2 <- ifelse(BigDataFrame$Herb >0, 0,1)
BigDataFrame$Microbe2 <- ifelse(BigDataFrame$Microbe >0, 0,1)

# Remember to change this before each run
write.table(BigDataFrame, file = "extinctionSetBRH_L.csv", sep = ",", row = F)

#  =========== For the Frequency graphs =====================

allData2 <- lapply(PopFiles, function (x){
  dat <- fread(x, header = FALSE, col.names = c("Generation", "Plant", "Herbivore", "Microbe", "Edit_PH", "Edit_PM", "Edit_MH", "Euclid_PH", "Euclid_PM","Euclid_MH", "Var_Plant", "Var_Herb","Var_Micro"))
  parts = strsplit(x, "_")[[1]]
  dat$run = as.factor(parts[3])
  
  return(dat)
}) 

BigDataFrame2 <- rbindlist(allData2)

PlantBDF <-summarySE(BigDataFrame2, measurevar = "Plant", groupvars = c("Generation"), na.rm = TRUE)
HerbBDF <-summarySE(BigDataFrame2, measurevar = "Herbivore", groupvars = c("Generation"), na.rm = TRUE)
MicrBDF <-summarySE(BigDataFrame2, measurevar = "Microbe", groupvars = c("Generation"), na.rm = TRUE)

PlantBDF <- rename(PlantBDF, c("Plant" = "Frequency"))
HerbBDF <- rename(HerbBDF, c("Herbivore" = "Frequency"))
MicrBDF <- rename(MicrBDF, c("Microbe" = "Frequency"))

PlantBDF$Organism <- rep("Plant", length(PlantBDF$Generation))
HerbBDF$Organism <- rep("Herbivore", length(HerbBDF$Generation))
MicrBDF$Organism <- rep("Microbe", length(MicrBDF$Generation))

FreqMarch18 <- rbind(PlantBDF, HerbBDF, MicrBDF)

# Remember to change!!!!!!- Based on the Plant-Microbe distance, Plant-Herbivore distance, and Hypothesis being tested.
 
FreqMarch18$MicrobeDist <- rep("Near", length(FreqMarch18$Generation))
FreqMarch18$InitDistPH <- rep("0.52", length(FreqMarch18$Generation))
FreqMarch18$Hypothesis <- rep("ERH", length(FreqMarch18$Generation))

write.table(FreqMarch18, file = "FreqERH_C.csv", sep = ",", row = F)

# Distances Data

BigDataFrameD <- BigDataFrame[rowSums(!is.na(TestDF)) > 8,]

PH_Edit <- summarySE(BigDataFrameD, measurevar = "Edit_PH", groupvars = "Generation", na.rm = T)
PM_Edit <- summarySE(BigDataFrameD, measurevar = "Edit_PM", groupvars = "Generation", na.rm = T)
MH_Edit <- summarySE(BigDataFrameD, measurevar = "Edit_MH", groupvars = "Generation", na.rm = T)

PH_Euclid <- summarySE(BigDataFrameD, measurevar = "Euclid_PH", groupvars = "Generation", na.rm = T)
PM_Euclid <- summarySE(BigDataFrameD, measurevar = "Euclid_PM", groupvars = "Generation", na.rm = T)
MH_Euclid <- summarySE(BigDataFrameD, measurevar = "Euclid_MH", groupvars = "Generation", na.rm = T)


myds_PH_Euclid.1 <- rename(PH_Euclid, c("Euclid_PH" = "Euclid"))
myds_PM_Euclid.1 <- rename(PM_Euclid, c("Euclid_PM" = "Euclid"))
myds_MH_Euclid.1 <- rename(MH_Euclid, c("Euclid_MH" = "Euclid"))

myds_PH_Euclid.1$Btw <- rep("PH", length(myds_PH_Euclid.1$Generation))
myds_PM_Euclid.1$Btw <- rep("PM", length(myds_PM_Euclid.1$Generation))
myds_MH_Euclid.1$Btw <- rep("MH", length(myds_MH_Euclid.1$Generation))

mydsFullEuclid <- rbind(myds_PH_Euclid.1, myds_PM_Euclid.1, myds_MH_Euclid.1)

myds_PH_Edit.1 <- rename(PH_Edit, c("Edit_PH" = "Edit"))
myds_PM_Edit.1 <- rename(PM_Edit, c("Edit_PM" = "Edit"))
myds_MH_Edit.1 <- rename(MH_Edit, c("Edit_MH" = "Edit"))

myds_PH_Edit.1$Btw <- rep("PH", length(myds_PH_Edit.1$Generation))
myds_PM_Edit.1$Btw <- rep("PM", length(myds_PM_Edit.1$Generation))
myds_MH_Edit.1$Btw <- rep("MH", length(myds_MH_Edit.1$Generation))

mydsFullEdit <- rbind(myds_PH_Edit.1, myds_PM_Edit.1, myds_MH_Edit.1)


mydsFull <- merge(mydsFullEuclid, mydsFullEdit, by = c("Generation", "Btw", "N"))

mydsFull <- rename(mydsFull, c("sd.x" = "EuclidSD"))	
mydsFull <- rename(mydsFull, c("se.x" = "EuclidSE"))						  
mydsFull <- rename(mydsFull, c("ci.x" = "EuclidCI"))							  

mydsFull <- rename(mydsFull, c("sd.y" = "EditSD"))	
mydsFull <- rename(mydsFull, c("se.y" = "EditSE"))						  
mydsFull <- rename(mydsFull, c("ci.y" = "EditCI"))

mydsFull$InitDist <- rep(0.13, length(mydsFull$Generation))

write.table(mydsFull, file = "Distances_.csv", sep = ",", row = F) # Change

mydsFullSub <- subset(mydsFull, Generation == 1023| Generation == 55295 | Generation == 111615 | Generation == 166911 | Generation == 222207 | Generation == 277503 | Generation == 333823 | Generation == 388095 | Generation == 444415 | Generation == 499711)

write.table(mydsFullSub, file = "DistancesSub_.csv", sep = ",", row = F) # Change


# Variances Data
Var_PlantDF <- summarySE(BigDataFrameD, measurevar = "Var_Plant", groupvars = "Generation", na.rm = T)
Var_HerbDF <- summarySE(BigDataFrameD, measurevar = "Var_Herb", groupvars = "Generation", na.rm = T)
Var_MicroDF <- summarySE(BigDataFrameD, measurevar = "Var_Micro", groupvars = "Generation", na.rm = T)



myds_Var_PlantDF.1 <- rename(Var_PlantDF, c("Var_Plant" = "Variance"))
myds_Var_HerbDF.1 <- rename(Var_HerbDF, c("Var_Herb" = "Variance"))
myds_Var_MicroDF.1 <- rename(Var_MicroDF, c("Var_Micro" = "Variance"))

myds_Var_PlantDF.1$Btw <- rep("PH", length(myds_Var_PlantDF.1$Generation))
myds_Var_HerbDF.1$Btw <- rep("PM", length(myds_Var_HerbDF.1$Generation))
myds_Var_MicroDF.1$Btw <- rep("MH", length(myds_Var_MicroDF.1$Generation))

mydsFullVariance <- rbind(myds_Var_PlantDF.1, myds_Var_HerbDF.1, myds_Var_MicroDF.1)

write.table(mydsFullVariance, file = "Variances_.csv", sep = ",", row = F) # Change

mydsFullVarianceSub <- subset(mydsFullVariance, Generation == 1023| Generation == 55295 | Generation == 111615 | Generation == 166911 | Generation == 222207 | Generation == 277503 | Generation == 333823 | Generation == 388095 | Generation == 444415 | Generation == 499711)

write.table(mydsFullVarianceSub, file = "VariancesSub_.csv", sep = ",", row = F) # Change
