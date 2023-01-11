# create detection history for mark-recapture analysis of detection probability
# NOTE this program treats Drifts 6 & 15 (which had irregular duty cycles) as if they were not
#     duty-cycled.  This approximation may be wrong for some detections.  May not want to include
#     these two in the mark-recapture analysis.


# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
createDetectionHistory <- function(ea, MaxTime=61.999, meanDepth =1217, sdDepth=354,
                                   HP_depth = 110, truncDist = 4000, plot=FALSE, snapshot=2,
                                   dutyCycle=NULL) {
    TruncAngle <- atan(truncDist/(meanDepth-HP_depth))*180/pi  #actual truncation angle based on nominal truncation dist. and mean depth

    if(is.character(ea$FileDateTime)) {
        ea$UTCtime <- as.POSIXct(ea$FileDateTime, "%Y-%m-%d %H:%M:%S", tz='UTC')
    } else {
        ea$UTCtime <- ea$FileDateTime
    }
    # Why is this done and then undone in Jeff's code
    ea$MeanAngle <- 180 - ea$MeanAngle
    ea <- ea[ea$nClicks >= 3,]
    ea <- ea[!is.na(ea$MeanAngle),]

    # limit sample to Ziphius
    ea$DeltaMin[is.na(ea$DeltaMin)]= -9999
    nlines <- length(ea$DeltaMin)

    # plot distribution of detection angles and truncation angle
    if(plot) {
        op <- par(xaxs="i",yaxs="i",cex.axis=1.3,cex.lab=1.5,font=2,font.lab=2)  #format plots
        on.exit(par(op))
        hist(180-ea$MeanAngle,xlab="Detection Angle (deg)",ylim=c(0,160),breaks=seq(0,100,5),main=NULL)
        lines(c(TruncAngle,TruncAngle),c(0,200),lty="dashed")
    }

    # write ZC event file for input to density analysis
    ##### INPUT FOR JEFFS FUNCTiON encounterRate-DataPrep####
    # write.csv(ea,file="ZcEventFilesWithAngles_drifts_1-22 - CorrectedAngles v2.csv", row.names=FALSE)

    # limit samples to those that will be used for abundance estimation (truncation distance and minimum click count)
    ea_trunc <- ea[ea$MeanAngle>=(180-TruncAngle),]
    # estimate the number of snapshots with detections for each drift w/in trunc angle
    eventsPerStation <- group_by(ea_trunc, Station) %>%
        summarise(numEvents = n())

    # Search through Ziphius events and assign a unique number for each dive (using MaxTime limit)
    ea <- arrange(ea, Station, UTCtime)
    ea$DiveNum= NA
    ea$DiveNum[1]= 1
    ea$DiveTime= NA
    ea$DiveTime[1]= 0
    DiveStart= ea$UTCtime[1]
    Dnumber= 1
    for(i in 2:length(ea$DiveNum)) {
        ea$DiveTime[i]= as.numeric(difftime(ea$UTCtime[i],DiveStart,units="mins"))
        if((ea$Station[i] != ea$Station[i-1]) ||
           (ea$DiveTime[i] > MaxTime) ||
           (ea$DiveTime[i] < 0)) {
            Dnumber= Dnumber + 1
            ea$DiveTime[i]= 0
            DiveStart= ea$UTCtime[i]
        }
        ea$DiveNum[i]= Dnumber
    }
    cat(" total number of dives:",Dnumber,"\n")

    # Create Detection History dataframe
    ndives= max(ea$DiveNum)
    # is this div by 2 because 2 minutes?? And then is max time 61.999 so that div/2 is lower?
    DetectionHistory= array(0,dim=c(ndives,1+MaxTime/snapshot))
    Recorder= rep(NA,ndives)
    Station= rep(NA,ndives)
    Time= rep(NA,ndives)
    EventId= rep(NA,ndives)
    ea$FilesPerDive= NA
    DeltaTime= rep(NA,ndives)
    StartAngle= rep(NA,ndives); EndAngle= rep(NA,ndives)
    nAngles= 0; Angles=rep(0,2); AngleTimes= rep(0,2)
    # create dectection history for each recognized dive
    for (idive in 1:ndives) {
        DiveEvents <- ea[ea$DiveNum == idive,]
        n <- length(DiveEvents$MeanAngle)
        Recorder[idive] <- as.character(DiveEvents$Recorder[1])
        Station[idive] <- DiveEvents$Station[1]
        EventId[idive] <- as.character(DiveEvents$EventId[1])
        Time[idive] <- as.character(DiveEvents$UTCtime[1])
        DeltaTime[idive] <- difftime(DiveEvents$UTCtime[n],DiveEvents$UTCtime[1],units='min')
        StartAngle[idive] <- DiveEvents$MeanAngle[1]
        EndAngle[idive] <- DiveEvents$MeanAngle[n]
        for (i in 1:n) {
            elapsedMin= DiveEvents$DiveTime[i]
            if (elapsedMin < MaxTime) {
                DetectionHistory[idive,1+(elapsedMin/snapshot)]= DiveEvents$MeanAngle[i]
                nAngles= nAngles + 1
                Angles[nAngles]= DiveEvents$MeanAngle[i]
                AngleTimes[nAngles]= difftime(DiveEvents$UTCtime[i],DiveEvents$UTCtime[1],units='min')
            }
        }
    }

    # for instruments on duty cycle, indicate when detections were not possible with -1
    # so ST4300 were on 2/10 duty cycle, SM2Bat was on 2/4 duty cycle?
    # DetectionHistory[Recorder=="SM2Bat",seq(2,1+MaxTime/2,2)]= -1
    # DetectionHistory[(Recorder=="ST4300")&(Station!=6)&(Station!=15)&(Station<=22),seq(2,1+MaxTime/2,5)]= -1
    # DetectionHistory[(Recorder=="ST4300")&(Station!=6)&(Station!=15)&(Station<=22),seq(3,1+MaxTime/2,5)]= -1
    # DetectionHistory[(Recorder=="ST4300")&(Station!=6)&(Station!=15)&(Station<=22),seq(4,1+MaxTime/2,5)]= -1
    # DetectionHistory[(Recorder=="ST4300")&(Station!=6)&(Station!=15)&(Station<=22),seq(5,1+MaxTime/2,5)]= -1

    # x$Recorder <- 'ST4300'
    # x$Recorder[x$DriftNum %in% c(7, 10, 13, 17)] <- 'SM3M'
    # x$Recorder[x$DriftNum %in% c(8, 11, 14, 18)] <- 'SM2Bat'

    if(is.null(dutyCycle)) {
        dutyCycle <- tribble(
            ~Station, ~dutyCycle,
            1, '2/8',
            2, '2/8',
            3, '2/8',
            4, '2/8',
            5, '2/8',
            6, '2/0',
            7, '2/0',
            8, '2/2',
            9, '2/8',
            10, '2/0',
            11, '2/2',
            12, '2/8',
            13, '2/0',
            14, '2/2',
            15, '2/0',
            16, '2/8',
            17, '2/0',
            18, '2/2',
            19, '2/8',
            20, '2/8',
            21, '2/8',
            22, '2/8'
        )
    }
    for(s in unique(Station)) {
        stationOff <- markOffCycle(dutyCycle$dutyCycle[dutyCycle$Station == s])
        DetectionHistory[Station == s, stationOff] <- -1
    }
    # create dataframe from detection history matrix and save as csv file
    DetectionHistoryOut= data.frame(Station,Time,EventId,Recorder,DetectionHistory)
    DetectionHistoryOut= DetectionHistoryOut[!is.na(DetectionHistoryOut$Recorder),]
    #### OUTPUT FOR JEFF EH-dataPrep ####
    # write.csv(DetectionHistoryOut,file="DetectionHistory.csv", row.names=FALSE)

    if(plot) {
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
    }
    # Angles= 180 - Angles
    # plot(Angles[1:(nAngles-1)],Angles[2:nAngles],type="l")
    # plot(AngleTimes,Angles,type="l",xlab="Elapsed Time (min)")
    list(ea=ea,
         detHistory=DetectionHistoryOut)
}
