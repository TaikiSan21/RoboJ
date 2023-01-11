# create detection history for mark-recapture analysis of detection probability
# NOTE this program treats Drifts 6 & 15 (which had irregular duty cycles) as if they were not 
#     duty-cycled.  This approximation may be wrong for some detections.  May not want to include
#     these two in the mark-recapture analysis.
par(xaxs="i",yaxs="i",cex.axis=1.3,cex.lab=1.5,font=2,font.lab=2)  #format plots

# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
setwd("C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/Files for Submission/Supplemental R-code & data")

  MaxTime= 61.999  #maximum encounter time in minutes
  meanDepth<- 1217    #mean depth from Barlow et al. 2020 tracking study in waters >2000m
  sdDepth<- 354      #sd depth from same source
  HP_depth<- 110     #depth at mid-point of hydrohone pair
  truncDist<- 4000   #nominal horizontal truncation distance
  TruncAngle<- atan(truncDist/(meanDepth-HP_depth))*180/pi  #actual truncation angle based on nominal truncation dist. and mean depth
  
  
  AllEventFilesWithAngle= read.csv(file="AllEventAngles_drifts_1-22.csv")
  AllEventFilesWithAngle$UTCtime= strptime(as.character(AllEventFilesWithAngle$FileDateTime), "%Y-%m-%d %H:%M:%S", tz="GMT")
  AllEventFilesWithAngle$MeanAngle= 180 - AllEventFilesWithAngle$MeanAngle
  AllEventFilesWithAngle= AllEventFilesWithAngle[AllEventFilesWithAngle$nClicks >= 3,]
  AllEventFilesWithAngle= AllEventFilesWithAngle[!is.na(AllEventFilesWithAngle$MeanAngle),]
  summary(as.factor(AllEventFilesWithAngle$eventType))
  BeakedWhales= c("ZC","MS","BWunid","BB","BW39V","BW43","BW","BW26-47","BW38")
  BWEventFilesWithAngle= AllEventFilesWithAngle[AllEventFilesWithAngle$eventType %in% BeakedWhales,]
  summary(as.factor(BWEventFilesWithAngle$eventType))
  sum(summary(as.factor(BWEventFilesWithAngle$eventType)))
  
# limit sample to Ziphius
  ZCevents= AllEventFilesWithAngle[AllEventFilesWithAngle$eventType=="ZC",]
  ZCevents$DeltaMin[is.na(ZCevents$DeltaMin)]= -9999
  nlines= length(ZCevents$DeltaMin)


# eliminate drifts on seamount ... need to do more QA-QC on those.  This program is not working well in defining dives.
  ZCevents= ZCevents[ZCevents$Station<=22,]

# plot distribution of detection angles and truncation angle
  hist(180-ZCevents$MeanAngle,xlab="Detection Angle (deg)",ylim=c(0,160),breaks=seq(0,100,5),main=NULL)
  lines(c(TruncAngle,TruncAngle),c(0,200),lty="dashed")
  
# write ZC event file for input to density analysis
  # write.csv(ZCevents,file="ZcEventFilesWithAngles_drifts_1-22 - CorrectedAngles v2.csv")
  
# limit samples to those that will be used for abundance estimation (truncation distance and minimum click count)
  ZCevents_trunc= ZCevents[ZCevents$MeanAngle>=(180-TruncAngle),]
  # estimate the number of snapshots with detections for each drift w/in trunc angle
  NumEventsPerDrift= array()
  for (iDrift in 1:22) {
    NumEventsPerDrift[iDrift]= sum(ZCevents_trunc$Station==iDrift)
  }
  ZCeventsPerDrift= data.frame(Drift=1:22,ZCevents=NumEventsPerDrift)
  write.csv(ZCeventsPerDrift,file="ZCeventsPerDrift.csv")
  
  
# Search through Ziphius events and assign a unique number for each dive (using MaxTime limit)
  ZCevents$DiveNum= NA
  ZCevents$DiveNum[1]= 1
  ZCevents$DiveTime= NA
  ZCevents$DiveTime[1]= 0
  DiveStart= ZCevents$UTCtime[1]
  Dnumber= 1
  for (i in 2:length(ZCevents$DiveNum)) {
    ZCevents$DiveTime[i]= difftime(ZCevents$UTCtime[i],DiveStart,units="mins")
    if ((ZCevents$DiveTime[i] > MaxTime)| (ZCevents$DiveTime[i] < 0)) {
      Dnumber= Dnumber + 1
      ZCevents$DiveTime[i]= 0
      DiveStart= ZCevents$UTCtime[i]
    }
    ZCevents$DiveNum[i]= Dnumber  
  }
  cat(" total number of dives:",Dnumber,"\n")

# Create Detection History dataframe
  ndives= max(ZCevents$DiveNum)
  DetectionHistory= array(0,dim=c(ndives,1+MaxTime/2))
  Recorder= rep(NA,ndives)
  Station= rep(NA,ndives)
  Time= rep(NA,ndives)
  EventId= rep(NA,ndives)
  ZCevents$FilesPerDive= NA
  DeltaTime= rep(NA,ndives)
  StartAngle= rep(NA,ndives); EndAngle= rep(NA,ndives)
  nAngles= 0; Angles=rep(0,2); AngleTimes= rep(0,2)
# create dectection history for each recognized dive
  for (idive in 1:ndives) {  
    DiveEvents= ZCevents[ZCevents$DiveNum == idive,]
    n= length(DiveEvents$MeanAngle)
    Recorder[idive]= as.character(DiveEvents$Recorder[1])
    Station[idive]= DiveEvents$Station[1]
    EventId[idive]= as.character(DiveEvents$EventId[1])
    Time[idive]= as.character(DiveEvents$UTCtime[1])
    DeltaTime[idive]= difftime(DiveEvents$UTCtime[n],DiveEvents$UTCtime[1],unit='min')
    StartAngle[idive]= DiveEvents$MeanAngle[1]
    EndAngle[idive]=   DiveEvents$MeanAngle[n]
    for (i in 1:n) {
      elapsedMin= DiveEvents$DiveTime[i]
      if (elapsedMin < MaxTime) {
          DetectionHistory[idive,1+(elapsedMin/2)]= DiveEvents$MeanAngle[i]
          nAngles= nAngles + 1
          Angles[nAngles]= DiveEvents$MeanAngle[i]
          AngleTimes[nAngles]= difftime(DiveEvents$UTCtime[i],DiveEvents$UTCtime[1],unit='min')
      }
    }
  }

# for instruments on duty cycle, indicate when detections were not possible with -1
  DetectionHistory[Recorder=="SM2Bat",seq(2,1+MaxTime/2,2)]= -1
  DetectionHistory[(Recorder=="ST4300")&(Station!=6)&(Station!=15)&(Station<=22),seq(2,1+MaxTime/2,5)]= -1
  DetectionHistory[(Recorder=="ST4300")&(Station!=6)&(Station!=15)&(Station<=22),seq(3,1+MaxTime/2,5)]= -1
  DetectionHistory[(Recorder=="ST4300")&(Station!=6)&(Station!=15)&(Station<=22),seq(4,1+MaxTime/2,5)]= -1
  DetectionHistory[(Recorder=="ST4300")&(Station!=6)&(Station!=15)&(Station<=22),seq(5,1+MaxTime/2,5)]= -1
  
# create dataframe from detection history matrix and save as csv file
  DetectionHistoryOut= data.frame(Station,Time,EventId,Recorder,DetectionHistory)
  DetectionHistoryOut= DetectionHistoryOut[!is.na(DetectionHistoryOut$Recorder),]
  write.csv(DetectionHistoryOut,file="DetectionHistory.csv")

# determine duration of encounter w/o angle truncation
  DetectDuration= rep(0,ndives)
  for (idive in 1:ndives) {
    DetectDuration[idive]= (max(which(DetectionHistory[idive,]>0)) - min(which(DetectionHistory[idive,]>0)) + 1) * 2
  }
  hist(DetectDuration[Recorder=="SM3M"],breaks=seq(0,80,2),xlab="SM3M Encounter Duration (min)",
       main=paste("Mean =",mean(DetectDuration[Recorder=="SM3M"]),"  SD=",sd(DetectDuration[Recorder=="SM3M"])))
  hist(DetectDuration[Recorder=="SM2Bat"],breaks=seq(0,80,2),xlab="SM2Bat Encounter Duration (min)",
       main=paste("Mean =",mean(DetectDuration[Recorder=="SM2Bat"]),"  SD=",sd(DetectDuration[Recorder=="SM2Bat"])))
  hist(DetectDuration[Recorder=="ST4300"],breaks=seq(0,80,2),xlab="ST4300 Encounter Duration (min)",
       main=paste("Mean =",mean(DetectDuration[Recorder=="ST4300"]),"  SD=",sd(DetectDuration[Recorder=="ST4300"])))
  
# determine duration of encounter w/ angle truncation
  DetectDurationTrunc= rep(0,ndives)
  DetectionHistoryTrunc= DetectionHistory
  DetectionHistoryTrunc[(180-DetectionHistoryTrunc) > TruncAngle]= 0
  for (idive in 1:ndives) {
    ndet= sum(DetectionHistoryTrunc[idive,]>0)
    if (ndet > 0) {
      DetectDurationTrunc[idive]= duration= (max(which(DetectionHistoryTrunc[idive,]>0)) - min(which(DetectionHistoryTrunc[idive,]>0)) + 1) * 2
    }
  }
  hist(DetectDurationTrunc[Recorder=="SM3M"&DetectDurationTrunc>0],breaks=seq(0,80,2),xlab="SM3M Encounter Duration (min)",
       main=paste("Mean =",mean(DetectDurationTrunc[Recorder=="SM3M"&DetectDurationTrunc>0],na.rm=T),"  SD=",
       sd(DetectDurationTrunc[Recorder=="SM3M"&DetectDurationTrunc>0],na.rm=T)))
  hist(DetectDurationTrunc[Recorder=="SM2Bat"&DetectDurationTrunc>0],breaks=seq(0,80,2),xlab="SM2Bat Encounter Duration (min)",
       main=paste("Mean =",mean(DetectDurationTrunc[Recorder=="SM2Bat"&DetectDurationTrunc>0],na.rm=T),"  SD=",
       sd(DetectDurationTrunc[Recorder=="SM2Bat"&DetectDurationTrunc>0],na.rm=T)))
  hist(DetectDurationTrunc[Recorder=="ST4300"&DetectDurationTrunc>0],breaks=seq(0,80,2),xlab="ST4300 Encounter Duration (min)",
       main=paste("Mean =",mean(DetectDurationTrunc[Recorder=="ST4300"&DetectDurationTrunc>0],na.rm=T),"  SD=",
       sd(DetectDurationTrunc[Recorder=="ST4300"&DetectDurationTrunc>0],na.rm=T)))
  
  plot(DeltaTime,(StartAngle-EndAngle),type="p")
  abline(h=0)
  plot(StartAngle,EndAngle,type="p",xlim=c(0,90),ylim=c(0,90))
  lines(c(0,90),c(0,90),col="red")
  plot(StartAngle,(StartAngle-EndAngle),type="p")
  abline(h=0)
  
  # Angles= 180 - Angles
  # plot(Angles[1:(nAngles-1)],Angles[2:nAngles],type="l")  
  # plot(AngleTimes,Angles,type="l",xlab="Elapsed Time (min)")
  
  