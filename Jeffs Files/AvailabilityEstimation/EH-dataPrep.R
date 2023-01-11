# 2021 April 8, Jeff Moore
# Called by EH-runBUGS.R
# First section of this file reads in Encounter History file (from DASBR data) and preps it for BUGS analysis to estimate availability time (creates BUGS data object)
# Second section of this generates a suite of data summary plots
# For an analysis of new data, code will have to be carefully reviewed line by line, because some steps are dataset specific
# Notes:
    # Distances are in km unless otherwise specified
    # for finding average detection angles or distances (when there are multiple angles for a group), we average on the angle tangents, i.e., arctan{average[tan(angle*pi/180)]} * 180/pi
    # in input .csv file, angles range from ~ 90 (horizontal) to ~ 180 (straight down).
    # in the code below, we use the transformation (180 - measured angle), so that 90 deg = horizontal and 0 deg = straight down
# Citations:
    # Barlow et al. 2020?

######################################################
# Input and manipulate data, create BUGS data object #
######################################################

##### set working directory #####
setwd("C:/jeff/NOAA/PASCAL/analysis/finalAnalysisDataAndCode-March2021/AvailabilityEstimation")

##### Ancillary info and specifications #####
mean.depth = 1.217  # mean dive depth (km), from Barlow et al. 2020
station.depth = 0.113  # mean depth of hydrophone array (km)
TRUNC = 4  # truncation distance (km)

##### Read in data file #####
ehdata = read.csv("DetectionHistoryWithAngles.csv")

##### Data manipulation #####
# Reduce entire dataframe to matrix with unnecessary fields and records removed
# Matrix will consist of detection angles for detected animals, 0s for non-detections where detectors were on, and NAs where detectors were off
  intervalRows = 6:36  # indicate which fields in the data spreadsheet are the actual encounter data (keep fields labeled x1 to x31 and exclude covariates)
  intMatrix = as.matrix(ehdata[,intervalRows])  # remove the covariates, leaving just the encounter history matrix
  intMatrix[intMatrix==-1] <- NA  # if duty-cycle-off was indicated by -1 in the csv, this line changes those to NAs (we want NAs for occasions where detectors were off)
  intMatrix = intMatrix[!ehdata$Station==6 & !ehdata$Station==15,]  # Remove rows for DASBRS 6 and 15, due to duty cycling errors on these DASBRs

##### Create other matrrix objects ##### 
# empty matrices of correct dimension
  detMatrix = effMatrix = posAngleMatrix = intMatrix  
# binary detection matrix (changes detection angles to 1s, changes NAs to 0s)
  detMatrix[intMatrix>0]=1;detMatrix[is.na(intMatrix)]=0  # if duty-cycle off are NAs in the csv:  matrix that changes angles to 1 (binary indicator of detection) and 0s otherwise
# effort matrix (1 if detector on, 0 if off due to duty cycling)
  effMatrix[!is.na(effMatrix)]=1; effMatrix[is.na(effMatrix)]=0  # (if duty-cycle off are NAs in the csv) matrix of 1s for duty-cyle on, 0s for duty-cycle off
# matrix of positive measure angles (matrix has angles to positive detections and NAs otherwise)
  posAngleMatrix[posAngleMatrix==0] = NA  
# in a new object, replace 1s from detMatrix with the interval of detection (used for calculating covariate "lastInt" below)
  detMatrix2 = detMatrix
   for(i in 1:nrow(detMatrix)){
     detMatrix2[i,] = detMatrix[i,]*c(1:ncol(detMatrix))
   }

##### Derived data and covariate fields #####
# Estimated distance to initial detection (based on angle to initial encounter), assuming animal is at its mean dive depth
  DetDist0.km = (mean.depth-station.depth) * tan((180-intMatrix[,1])*pi/180)
# Mean detection distance (average of all detection distances in the encounter)
  DetDistMean.km = (mean.depth-station.depth) * rowMeans(tan((180-posAngleMatrix)*pi/180),na.rm=T)
# Difference between initial and mean detection
  DistDiff = DetDist0.km - DetDistMean.km
# Number of detections per encounter history
  numDet = rowSums(detMatrix)
# Number of intervals during which the detector was on (duty cycling), a measure of effort
  numInts = rowSums(effMatrix)
# crude measure of detection rate: proportion of intervals detected out of the number of on-intervals
  dpue = numDet/numInts
# last interval in which animal was detection (a measure of encounter length)
  lastInt = apply(detMatrix2,1,max)
# alternate crude measure of detection rate: proportion of intervals detected within the period that it was being heard
  dpue2 = numDet/lastInt

##### Create BUGS data variables and object #####
  # only retain records for which initial detection distance is within the truncation distance
  # but note that covariate for analysis is the mean detection distance (for those records where initial distance was < 4km)
    # sample size (number of encounter histories)
      ndives = nrow(detMatrix[DetDist0.km<TRUNC,])
    # number of encounter history occasions
      nocc = ncol(intMatrix[DetDist0.km < TRUNC,])
    # index whether DASRR is on 100% of time (dutyidx=1) or duty-cycled (dutyidx=0)
      dutyidx = rep(NA,ndives)
        for(i in 1:ndives){
          dutyidx[i] = 1 - (sum(effMatrix[DetDist0.km < TRUNC,][i, ]==0)>0)
        }
    # BUGS data object
      bugs.data.EH = list(nocc=nocc, nind=ndives, eh=detMatrix[DetDist0.km<TRUNC,], dcyc=effMatrix[DetDist0.km<TRUNC,])
  

#####################################################
# Generate (and save) some plots and data summaries #
#####################################################     
      
##### set working directory for where plots will be saved #####
setwd("C:/jeff/NOAA/PASCAL/analysis/finalAnalysisDataAndCode-March2021/dataSummaryPlots")

# Plot initial vs. mean detection distances, all data
x11()
par(mfrow=c(2,2))
hist(DetDist0.km, breaks=50, main=paste("n = ", length(DetDist0.km),sep="")); hist(DetDistMean.km,breaks=50, main=paste("n = ", length(DetDist0.km),sep="")); plot(DetDist0.km,DetDistMean.km); abline(a=0,b=1); hist(DistDiff, breaks=50)
savePlot(filename="Detection Distance Summaries.tiff", type="tiff")
dev.off()

# Plot initial vs. mean detection distances, data truncated to MEAN detection distance < TRUNC distance
x11()
par(mfrow=c(2,2))
hist(DetDist0.km[DetDistMean.km<TRUNC], breaks=25, main=paste("n = ", length(DetDist0.km[DetDistMean.km<TRUNC]),sep="")); hist(DetDistMean.km[DetDistMean.km<TRUNC],breaks=25, main=paste("n = ", length(DetDistMean.km[DetDistMean.km<TRUNC]),sep="")); plot(DetDist0.km[DetDistMean.km<TRUNC],DetDistMean.km[DetDistMean.km<TRUNC], main=""); abline(a=0,b=1); hist(DistDiff[DetDistMean.km<TRUNC], breaks=50)
savePlot(filename="Detection Distance Summaries - trunc meanDist 4km.tiff", type="tiff")
dev.off()

# Plot initial vs. mean detection distances, data truncated to INITIAL detection distance < TRUNC distance
x11()
par(mfrow=c(2,2))
hist(DetDist0.km[DetDist0.km<TRUNC], breaks=25, main=paste("n = ", length(DetDist0.km[DetDist0.km<TRUNC]),sep="")); hist(DetDistMean.km[DetDist0.km<TRUNC],breaks=25, main=paste("n = ", length(DetDistMean.km[DetDist0.km<TRUNC]),sep="")); plot(DetDist0.km[DetDist0.km<TRUNC],DetDistMean.km[DetDist0.km<TRUNC], main=""); abline(a=0,b=1); hist(DistDiff[DetDist0.km<TRUNC], breaks=50)
savePlot(filename="Detection Distance Summaries - trunc initialDist 4km.tiff", type="tiff")
dev.off()

### Repeat some of the above plots, just for detectors that were never duty cycled

# Plot initial vs. mean detection distances, data truncated to MEAN detection distance < TRUNC distance
x11()
par(mfrow=c(2,2))
hist(DetDist0.km[DetDistMean.km<TRUNC & numInts==max(numInts) ], breaks=25, xlab="Initial Distance (km)", main=paste("n = ", length(DetDist0.km[DetDistMean.km<TRUNC & numInts==max(numInts)]),sep="")); hist(DetDistMean.km[DetDistMean.km<TRUNC & numInts==max(numInts)],breaks=25, xlab="Mean Distance (km)", main=paste("n = ", length(DetDistMean.km[DetDistMean.km<TRUNC & numInts==max(numInts)]),sep="")); plot(DetDist0.km[DetDistMean.km<TRUNC & numInts==max(numInts)],DetDistMean.km[DetDistMean.km<TRUNC & numInts==max(numInts)], xlab="Initial detection km", ylab="Mean detection km", main=""); abline(a=0,b=1); hist(DistDiff[DetDistMean.km<TRUNC & numInts==max(numInts)], breaks=50, xlab="Initial distance - Mean distance", main="")
savePlot(filename="Detection Distance Summaries - trunc meanDist 4km - noDutyCycle.tiff", type="tiff")
dev.off()

# Plot initial vs. mean detection distances, data truncated to INITIAL detection distance < TRUNC distance
x11()
par(mfrow=c(2,2))
hist(DetDist0.km[DetDist0.km<TRUNC & numInts==max(numInts) ], breaks=25, xlab="Initial Distance (km)", main=paste("n = ", length(DetDist0.km[DetDist0.km<TRUNC & numInts==max(numInts)]),sep="")); hist(DetDistMean.km[DetDist0.km<TRUNC & numInts==max(numInts)],breaks=25, xlab="Mean Distance (km)", main=paste("n = ", length(DetDistMean.km[DetDist0.km<TRUNC & numInts==max(numInts)]),sep="")); plot(DetDist0.km[DetDist0.km<TRUNC & numInts==max(numInts)],DetDistMean.km[DetDist0.km<TRUNC & numInts==max(numInts)], xlab="Initial detection km", ylab="Mean detection km", main=""); abline(a=0,b=1); hist(DistDiff[DetDist0.km<TRUNC & numInts==max(numInts)], breaks=50, xlab="Initial distance - Mean distance", main="")
savePlot(filename="Detection Distance Summaries - trunc initDist 4km - noDutyCycle.tiff", type="tiff")
dev.off()

### Plot mean detection distance vs. some encounter rate metrics
axisSize = 1.5
labSize = 1.75
fontType = 2 # bold
#x11()
tiff(filename = "Encounter rate metrics vs mean distance.tif", width=10, height=10,units="in",res=300)
par(mfrow=c(2,2), mar=c(5,5,4,2))
# Set 1: all data within truncation range
plot(DetDistMean.km[DetDist0.km<TRUNC], dpue[DetDist0.km<TRUNC], main="", xlab="Mean detection distance (km)", ylab="Detections per on-interval", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
plot(DetDistMean.km[DetDist0.km<TRUNC], lastInt[DetDist0.km<TRUNC], main="", xlab="Mean detection distance (km)", ylab="Last interval detected", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
# Set 2: only for non duty-cycled detectors
plot(DetDistMean.km[DetDist0.km<TRUNC & numInts==max(numInts)], dpue[DetDist0.km<TRUNC & numInts==max(numInts)], main="", xlab="Mean detection distance (km)", ylab="Detections per on-interval", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
plot(DetDistMean.km[DetDist0.km<TRUNC & numInts==max(numInts)], lastInt[DetDist0.km<TRUNC & numInts==max(numInts)], main="", xlab="Mean detection distance (km)", ylab="Last interval detected", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
mtext(paste("All data within truncation distance (n = ", length(DetDistMean.km[DetDist0.km<TRUNC]), ")", sep=""), side=3, outer=T,cex=2, line= -3, font=2)
mtext(paste("Only from detectors not duty-cycled (n = ", sum(DetDist0.km<TRUNC & numInts==max(numInts)), ")", sep=""), side=3, outer=T, line=-33,cex=2,font=2)
#savePlot(filename = "Encounter rate metrics vs mean distance.tif", type="tiff")
dev.off()

# Plot initial detection distance vs. some encounter rate metrics
x11()
par(mfrow=c(2,2))
# Set 1: all data within truncation range
plot(DetDist0.km[DetDist0.km<TRUNC], dpue[DetDist0.km<TRUNC], main=paste("n = ", length(DetDist0.km[DetDist0.km<TRUNC]),sep=""), xlab="Initial detection distance", ylab="Detections per on-interval")
plot(DetDist0.km[DetDist0.km<TRUNC], lastInt[DetDist0.km<TRUNC], main="All data within truncation distance", xlab="Initial detection distance", ylab="Last interval detected")
# Set 2: only for non duty-cycled detectors
plot(DetDist0.km[DetDist0.km<TRUNC & numInts==max(numInts)], dpue[DetDist0.km<TRUNC & numInts==max(numInts)], main=paste("n = ", sum(DetDist0.km<TRUNC & numInts==max(numInts)),sep=""), xlab="Initial detection distance", ylab="Detections per on-interval")
plot(DetDist0.km[DetDist0.km<TRUNC & numInts==max(numInts)], lastInt[DetDist0.km<TRUNC & numInts==max(numInts)], main="Only from detectors not duty-cycled", xlab="Initial detection distance", ylab="Last interval detected")
savePlot(filename = "Encounter rate metrics vs initial distance.tif", type="tiff")
dev.off()

# duration length (last interval detected) on non duty-cycled detectors
mean(lastInt[DetDist0.km<TRUNC & numInts==max(numInts)])  # mean value for "last Interval", for non-duty cycled detectors, within truncation distance
# 7.57 for 4k trunction
length(lastInt[DetDist0.km<TRUNC & numInts==max(numInts)])  # sample size
# n=44 for 4k truncation


