#2021 April 8, Jeff Moore
# Goodness of fit analysis for encounter histories

##### libraries #####
library(Hmisc)  # for errbar

##### OVERVIEW OF SIMULAITON PROCESS #####
# GOF analysis, like EH analysis, is conditioned on first capture
# Option 1: simulate n=44 encounter histories.  This is the subset of encounter histories recorded on DASBRs that were not duty-cycled
# Option 2: simulate n=247 encounter histories.  This is the number of encounter histories in real dataset where the initial detection distance was less than 4km
# The two options are coded in separate sections below
# RUN CODE FOR **ONE** OF THESE OPTIONS, THEN SCROLL BELOW AND RUN CODE FOR PLOTTING (REQUIRES SOME FILE NAME CHANGES IN THE PLOT SECTION)

############################################################
##### OPTION 1: SIMULATE FOR DETECTORS NOT DUTY CYCLED #####
############################################################

##### INPUTS #####
nEH = 44  # number of encounter histories to generate (n=247 or 44) NROW OF MATRIX
occ = 31  # number of capture occasions, same as real data NCOL OF MATRIX
occ.plot = 20  # plotting parameter
phi = 0.892  # probability of remaining 'alive' (available) at each capture occasion (this is the model estimate)
p = 0.43 # probability of being detected, given available (this is the model estimate)

##### BRING IN THE REAL DATA #####
# first, read in relevant lines of code from EH-dataPrep.R if objects generated from that code are not already in workspace
  # i.e., read in lines to create detMatrix, DetDist0.km, and dutyidx
# then...
ehData = detMatrix[DetDist0.km<TRUNC,][dutyidx==1,]  # data for non-duty-cycled detectors (n=44)

##### SIMULATION #####

# Empty objects
  availMatrix.s = matrix(0,nEH,occ)              # availability matrix, indexes whether animal is 'available' at time t
  ehMatrix.s = detMatrix2.s = matrix(0,nEH,occ)  # encounter history matrix, indexes whether an animal is detected at time t
  meanDetections.s = rep(NA,1000)
  fd.numDet.s = matrix(0,1000,occ)
  fd.dur.s = matrix(0,1000,occ)
  lastInt.s = matrix(0,1000,nEH)

# Loop through 1000 simulated datasets
for(i in 1:1000){
  availMatrix.s[,1] = ehMatrix.s[,1] = 1                         # first column is initial capture, so this column is all 1's
  for(j in 1:nEH){                                               # for each encounter history
    for(t in 2:occ){                                             # for each t within j
      availMatrix.s[j,t] = rbinom(1,1,availMatrix.s[j,t-1]*phi)  # availability at t is a bernouli variable, given phi and availability at t-1
      ehMatrix.s[j,t] = availMatrix.s[j,t] * rbinom(1,1,p)       # detection at t is bernouli variable, given availbility at t
    } # t
    detMatrix2.s[j,] = ehMatrix.s[j,]*c(1:occ)                   # for finding encounter duration, record the interval in which each detection occurs
  } # j

  # simulated data summaries
  meanDetections.s[i] = mean(rowSums(ehMatrix.s))                  # avg number of detections per encounter history for simulated data
  fd.numDet.s[i,] = hist(rowSums(ehMatrix.s),breaks=0:occ)$counts  # freq distn of number detections
  lastInt.s[i,] = apply(detMatrix2.s,1,max)                        # encounter lengths
  fd.dur.s[i,] = hist(lastInt.s[i,],breaks=0:occ)$counts           # freq distn for encounter lengths

} # i


################################################
##### OPTION 2: SIMULATE FOR ALL DETECTORS #####
################################################

##### INPUTS #####
nEH = 247  # number of encounter histories to generate (n=247 or 44)
occ = 31  # number of capture occasions, same as real data
occ.plot = 31  # plotting parameter
phi = 0.892  # probability of remaining 'alive' (available) at each capture occasion (this is the model estimate)
p = 0.43 # probability of being detected, given available (this is the model estimate)

##### BRING IN THE REAL DATA #####
# first, read in relevant lines of code from EH-dataPrep.R if objects generated from that code are not already in workspace
# i.e., read in lines to create detMatrix, DetDist0.km, and dutyidx
# then...
ehData = detMatrix[DetDist0.km<TRUNC,]  # data for all detectors (n=247)

##### SIMULATION #####

# Empty objects
  availMatrix.s = matrix(0,nEH,occ)              # availability matrix, indexes whether animal is 'available' at time t
  ehMatrix.s = detMatrix2.s = matrix(0,nEH,occ)  # encounter history matrix, indexes whether an animal is detected at time t
  meanDetections.s = rep(NA,1000)
  fd.numDet.s = matrix(0,1000,occ)
  fd.dur.s = matrix(0,1000,occ)
  lastInt.s = matrix(0,1000,nEH)
  effMatrix.s = effMatrix[DetDist0.km < TRUNC, ]  # indicates where duty cycling occurs; requires effMatrix to have been defined by running through EHdataPrep.R

# Loop through 1000 simulated datasets
for(i in 1:1000){
  availMatrix.s[,1] = ehMatrix.s[,1] = 1                                       # first column is initial capture, so this column is all 1's
  for(j in 1:nEH){                                                             # for each encounter history
    for(t in 2:occ){                                                           # for each t within j
      availMatrix.s[j,t] = rbinom(1,1,availMatrix.s[j,t-1]*phi)                # availability at t is a bernouli variable, given phi and availability at t-1
      ehMatrix.s[j,t] = availMatrix.s[j,t] * effMatrix.s[j,t] * rbinom(1,1,p)  # detection at t is bernouli variable, given availbility at t. Note effMatrix.s indicates duty cycling
    } # t
    detMatrix2.s[j,] = ehMatrix.s[j,]*c(1:occ)  # for finding encounter duration
  } # j

  # simulated data summaries
  meanDetections.s[i] = mean(rowSums(ehMatrix.s))                  # avg number of detections per encounter history for simulated data
  fd.numDet.s[i,] = hist(rowSums(ehMatrix.s),breaks=0:occ)$counts  # freq distn of number detections
  lastInt.s[i,] = apply(detMatrix2.s,1,max)                        # encounter lengths
  fd.dur.s[i,] = hist(lastInt.s[i,],breaks=0:occ)$counts           # freq distn for encounter lengths

} # i


##########################################################
##### PLOTS, FOR OPTION 1 OR OPTION 2. READ NOTES!!! #####
##########################################################
# NOTE: Depending on option 1 or option 2, some code has to be changed for a couple of the plots below.  Specifically:
     # Change the ymax value
     # For the 5th plot below, pound out the lines not being used (these are annotated)
# ALSO NOTE: plots are saved in "...out/GOF", and output names for the two Options are not differentiated.  Thus, copy outputs to a new subfolder (within GOF) that is named to indicate Option 1 or 2

# set this for plotting
ymax = 150    # 20 for n=44 (Option 1), 150 for n=247 (Option 2)

# 1. compare mean number of detections per encounter history, for real data vs. simulated data
#x11()
tiff(filename="out/GOF/numberDetectionsPerEH-mean.tif", width=4, height=4, units="in", res=300)
meanDetections.r =  mean(rowSums(ehData)) ;meanDetections.r  # real data
hist(meanDetections.s, xlab="Number of detections per encounter history", main="", font.lab=2); abline(v=meanDetections.r, col="red")  # distribution of means for the simulated datasets
#savePlot(filename="out/GOF/numberDetectionsPerEH-mean.jpg",type="jpg")
prob.Higher.Val = sum(meanDetections.s > meanDetections.r)/1000; prob.Higher.Val  # what percentage of simulated means are greater than the observed value?
mean(meanDetections.s)
dev.off()

# 2. compare frequency distribution for number of detections (real data vs. simulated data)
#x11()
tiff(filename="out/GOF/numberDetectionsPerEH-hist.tif", width=4, height=4, units="in", res=300)
fd.numDet.r = hist(rowSums(ehData), breaks=0:occ)$counts; fd.numDet.r
# Plot frequency distributions
errbar(1:occ.plot, y=colMeans(fd.numDet.s)[1:occ.plot], yplus=apply(fd.numDet.s,2,quantile,0.975)[1:occ.plot], yminus=apply(fd.numDet.s,2,quantile,0.025)[1:occ.plot], ylim=c(0,ymax), xlab="Number of occasions detected", ylab="Frequency", errbar.col="gray", col="gray", font.lab=2, cex=0.75)
points(1:occ.plot,fd.numDet.r[1:occ.plot], pch=19, cex=0.75)
#savePlot(filename="out/GOF/numberDetectionsPerEH-hist.jpg", type="jpg")
dev.off()

# 3. compare frequency distn of data with that of some random simulations
#x11()
tiff(filename="out/GOF/numberDetectionsPerEH-examples.tif", width=10, height=10, units="in", res=300)
par(mar=c(3,3,3,3))
par(mfrow=c(3,3))
plot(1:occ.plot, fd.numDet.r[1:occ.plot], pch=19, ylim=c(0,ymax),xlab="", ylab="", cex=2, cex.axis=2)
for(i in 1:8){
  plot(1:occ.plot,fd.numDet.s[i,1:occ.plot], ylim=c(0,ymax), xlab="", ylab="", cex=2, cex.axis=2)
}
#savePlot(filename="out/GOF/numberDetectionsPerEH-examples.jpg", type="jpg")
dev.off()

# 4. compare encounter duration (time of last detection), for real vs. simulated data
#x11()
tiff(filename="out/GOF/encounterDuration-mean.tif", width=4, height=4, units="in", res=300)
meanDuration.r = mean(lastInt[DetDist0.km<TRUNC]); meanDuration.r # mean encounter duration (interval of last detection) for the real data
mean(rowMeans(lastInt.s))
hist(rowMeans(lastInt.s), xlab="Last interval detected", main="", font.lab=2); abline(v=meanDuration.r, col="red")  # distribution of means for the simulated datasets
#savePlot(filename="out/GOF/encounterDuration-mean.jpg", type="jpg")
prob.higher.value = (sum(rowMeans(lastInt.s) > meanDuration.r)/1000); prob.higher.value
dev.off()

# 5. compare frequency distribution for number of detections (real data vs. simulated data)
#x11()
tiff(filename="out/GOF/encounterDuration-hist.tif", width=4, height=4, units="in", res=300)
# for n=44 case (Option 1)
#     fd.duration.r = hist(lastInt[DetDist0.km<TRUNC & numInts==max(numInts)], breaks=c(0:31), main="", xlab="Last interval detected")$counts; fd.duration.r
# for n=247 case (Option 2)
     fd.duration.r = hist(lastInt[DetDist0.km<TRUNC], breaks=c(0:31), main="", xlab="Last interval detected")$counts; fd.duration.r
# Plot freq distns
errbar(1:occ, y=colMeans(fd.dur.s)[1:occ], yplus=apply(fd.dur.s,2,quantile,0.975)[1:occ], yminus=apply(fd.dur.s,2,quantile,0.025)[1:occ], ylim=c(0,ymax), xlab="Last interval detected", ylab="Frequency", errbar.col="gray", col="gray", font.lab=2, cex=0.75 )
points(1:occ,fd.duration.r[1:occ], pch=19, cex=0.75)
#savePlot(filename="out/GOF/encounterDuration-hist.jpg", type="jpg")
dev.off()

# 6. compare frequency distn of duration-data with that of some random simulations
#x11()
tiff(filename="out/GOF/encounterDuration-examples.tif", width=10, height=10, units="in", res=300)
par(mar=c(3,3,3,3))
par(mfrow=c(3,3))
plot(1:occ, fd.duration.r[1:occ], pch=19, ylim=c(0,ymax),xlab="", ylab="", cex=2, cex.axis=2)
for(i in 1:8){
  plot(1:occ,fd.dur.s[i,1:occ], ylim=c(0,ymax), xlab="", ylab="", cex=2, cex.axis=2)
}
#savePlot(filename="out/GOF/encounterDuration-examples.jpg", type="jpg")
dev.off()
