# 2021 April 8, Jeff Moore
# Called by encounterRate-runBUGS.R
# First section of this file reads in encounter rate data (from DASBRs) and preps it for BUGS analysis to estimate availability time (creates BUGS data object)
# Second section of this generates some data summary plots
# For an analysis of new data, code will have to be carefully reviewed line by line, because some steps are dataset specific
# Notes:
# Distances are in km unless otherwise specified
# in input .csv file, angles range from ~ 90 (horizontal) to ~ 180 (straight down).
# in the code below, we use the transformation (180 - measured angle), so that 90 deg = horizontal and 0 deg = straight down

######################################################
# Input and manipulate data, create BUGS data object #
######################################################

##### set working directory #####
setwd("C:/jeff/NOAA/PASCAL/analysis/finalAnalysisDataAndCode-March2021/DensityEstimation")

##### Set some inputs #####
W.trunc = 4         # 4km truncation distance for including in analysis
mean.depth = 1.217  # mean dive depth (km). Enter same value here as in EH-dataPrep.R
#A = 1142000        # size of CA Current study area (km^2) out to 300 nmi
A = 1057925         # size of study area > 500m deep

##### Effort data #####
#  number of 2-min intervals on each DASBR (stations 1 - 22)
#  from email from Jay dated Sept 8, 2017
k.raw = c(1629, 1630, 1762, 1841, 1918, 3078, 13669, 6807, 2812, 13996, 7073, 2801, 12133, 0, 3135, 2800, 14250, 6934, 2707, 3329, 3087, 3376)
k = k.raw[k.raw>0]            # remove station 14 since there was no effort there
k.tot = sum(k)                # total number of 2-min intervals in the dataset
n.sites = length(k)           # number of stations with effort data
stationNames = c(1:13,15:22)  # stations with data

##### DASBR depths #####
# vector of station depths
stationDepth = rep(0.115,length(k.raw))  # depth (km) of ST4300 DASBRS
stationDepth[c(8,11,18)] = 0.091         # depth (km) of SM2Bat DASBRs
stationDepth[c(7,10,13,17)] = 0.113      # depth (km) of SM3M DASBRS
stationDepth = stationDepth[k.raw>0]     # remove DASBR 14
# label station types for the 21 functioning detectors
#stationTypes = c("ST4300", "ST4300","ST4300","ST4300", "ST4300","ST4300","SM3M", "SM2Bat","ST4300", "SM3M","SM2Bat","ST4300","SM3M","ST4300","ST4300","SM3M","SM2Bat","ST4300","ST4300","ST4300", "ST4300")
stationTypes = c(3,3,3,3,3,3,2,1,3,2,1,3,2,3,3,2,1,3,3,3,3)  # numeric label instead of text names
# find mean array depth
#mean.arrayDepth = mean(stationDepth)    # average depth of the 21 functioning DASBRs (km)
mean.arrayDepth = 0.113                  # Just use this value instead of station-specific depths.  This approximation is good (approximation error is trivial relative to other sources of error)

##### READ IN THE ENCOUNTER RATE DATA #####
ziph.data.in = read.csv(here('Jeffs Files', 'DensityEstimation', "ZcEventFilesWithAngles_drifts_1-22 - CorrectedAngles v2.csv"))

##### DATA MANIPULATION #####
ziph.angles1 = 180 - ziph.data.in$MeanAngle                         # new variable: express angles so that 0 = animal below hydrophones, 90 = animal at surface (level w phones)
recorder.idx = rep(1,nrow(ziph.data.in))                            # new variable: indicate the recorder type for each row using a numeric label
recorder.idx[ziph.data.in$Recorder=="SM3M"]=2
recorder.idx[ziph.data.in$Recorder=="ST4300"]=3
ziph.data = cbind(ziph.data.in, ziph.angles1, recorder.idx)         # new array: bind the input array with the two new variables above
ziph.data2 = ziph.data[ziph.data$nClicks>=3,]                       # filter out detections based on fewer than 3 clicks (because detections with fewer clicks could not confidently be identified as Ziphius... this did not eliminate much data)
trunc.angle = atan(W.trunc/(mean.depth-mean.arrayDepth)) *180/pi    # find truncation angle (in radians), given specified truncation distance above (W.trunc)
ziph.data.trunc = ziph.data2[ziph.data2$ziph.angles1 <= trunc.angle & !is.na(ziph.data2$ziph.angles1),]     # filter out NAs detections beyond truncation angle
n.per.station = table(ziph.data.trunc$Station)                      # find the number of detections per station
detPerStation = cbind(stationNames, rep(0,n.sites))                 # new matrix:  column 1 = station number (1-21), column 2 = number detections per station
detPerStation[stationNames %in% names(n.per.station),2] = n.per.station
colnames(detPerStation) = c("station", "detections")

##### Create BUGS data variables and object #####
#importtant thing is stationnames & k are in same order. Originally
# sn is used as an index, cant do that going forward. Just have to enforce
# same ordering of n, k
bugs.data.ER = list(n.sites=n.sites, n=detPerStation[,2], A=A, k=k, stationNames=stationNames)

#####################################################
# Generate (and save) some plots and data summaries #
#####################################################

# histogram of angles for retained detections
hist(ziph.data.trunc$ziph.angles1, breaks=50)

# histogram of distances
hist((mean.depth-mean.arrayDepth)*tan(ziph.data.trunc$ziph.angles1*pi/180),breaks=50, xlab="RANGE (KM)", main="DETECTION DISTANCES IF ANIMAL DEPTH = MEAN DEPTH")
n.trunc = nrow(ziph.data.trunc); n.trunc # sample size for detections (retained)

# histogram of detection angles by recorder type
x11()
par(mfrow = c(1,3))
hist(ziph.data.trunc$ziph.angles1[ziph.data.trunc$Recorder=="SM2Bat"],breaks=seq(0,90,10), xlab="Angle", main="SM2Bat")
hist(ziph.data.trunc$ziph.angles1[ziph.data.trunc$Recorder=="SM3M"],breaks=seq(0,90,10), xlab="Angle", main="SM3M")
hist(ziph.data.trunc$ziph.angles1[ziph.data.trunc$Recorder=="ST4300"],breaks=seq(0,90,10), xlab="Angle", main="ST4300")

