# RoboJeff

ehDataPrep <- function(detHist, meanDepth=1217, hpDepth=110, truncDist=4e3, plot=FALSE, plotDir='.', maxTime=61.999) {
    intMatrix <- detHist$dhMat
    intMatrix[intMatrix == -1] <- NA  # if duty-cycle-off was indicated by -1 in the csv, this line changes those to NAs (we want NAs for occasions where detectors were off)
    # Create interval length matrix
    lenMatrix <- matrix(0, nrow=nrow(intMatrix), ncol=ncol(intMatrix))
    for(i in 1:nrow(lenMatrix)) {
        lenMatrix[i, ] <- detHist$dcMap[[detHist$dhDf$station[i]]]
    }
    ##### Create other matrrix objects #####
    # empty matrices of correct dimension
    detMatrix <- effMatrix <- posAngleMatrix <- intMatrix

    # binary detection matrix (changes detection angles to 1s, changes NAs to 0s)
    # if duty-cycle off are NAs in the csv:  matrix that changes angles to 1 (binary indicator of detection) and 0s otherwise
    detMatrix[intMatrix>0] <- 1
    detMatrix[is.na(intMatrix)] <- 0

    # effort matrix (1 if detector on, 0 if off due to duty cycling)
    # (if duty-cycle off are NAs in the csv) matrix of 1s for duty-cyle on, 0s for duty-cycle off
    effMatrix[!is.na(effMatrix)] <- 1
    effMatrix[is.na(effMatrix)] <- 0

    # matrix of positive measure angles (matrix has angles to positive detections and NAs otherwise)
    posAngleMatrix[posAngleMatrix == 0] <- NA

    ##### Derived data and covariate fields #####
    # Estimated distance to initial detection (based on angle to initial encounter), assuming animal is at its mean dive depth
    detDist0_km <- (meanDepth-hpDepth) * tan((180-intMatrix[,1])*pi/180) / 1e3
    # Mean detection distance (average of all detection distances in the encounter)
    detDistMean_km <- (meanDepth-hpDepth) * rowMeans(tan((180-posAngleMatrix)*pi/180),na.rm=T) / 1e3
    # Difference between initial and mean detection
    distDiff <- detDist0_km - detDistMean_km
    # Number of detections per encounter history
    numDet <- rowSums(detMatrix)
    # Number of intervals during which the detector was on (duty cycling), a measure of effort
    numInts <- rowSums(effMatrix)
    # crude measure of detection rate: proportion of intervals detected out of the number of on-intervals
    dpue <- numDet/numInts
    # last interval in which animal was detection (a measure of encounter length)
    lastInt <- apply(detMatrix, 1, function(x) {
        isDet <- x > 0
        if(!any(isDet)) {
            return(0)
        }
        max(which(isDet))
    })
    # alternate crude measure of detection rate: proportion of intervals detected within the period that it was being heard
    dpue2 <- numDet/lastInt

    ##### Create BUGS data variables and object #####
    # only retain records for which initial detection distance is within the truncation distance
    # but note that covariate for analysis is the mean detection distance (for those records where initial distance was < 4km)
    # sample size (number of encounter histories)
    distOK <- detDist0_km < truncDist
    ndives <- nrow(detMatrix[distOK, ])
    # number of encounter history occasions
    nocc <- ncol(intMatrix[distOK, ])
    # index whether DASRR is on 100% of time (dutyidx=1) or duty-cycled (dutyidx=0)
    dutyidx <- rep(NA, ndives)
    for(i in 1:ndives){
        dutyidx[i] = 1 - (sum(effMatrix[distOK,][i, ]==0)>0)
    }
    # BUGS data object
    bugsDataEH <- list(nOcc=nocc,
                       nInd=ndives,
                       eh=detMatrix[distOK, ],
                       dCyc=effMatrix[distOK, ],
                       intLen = lenMatrix[distOK, ],
                       nMin = ceiling(maxTime))
    if(plot) {
        if(!dir.exists(plotDir)) {
            dir.create(plotDir)
        }
        # Plot initial vs. mean detection distances, all data
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km, breaks=50, main=paste("n = ", length(detDist0_km),sep=""))
        hist(detDistMean_km,breaks=50, main=paste("n = ", length(detDist0_km),sep=""))
        plot(detDist0_km,detDistMean_km); abline(a=0,b=1); hist(distDiff, breaks=50)
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries.tiff"), type="tiff")
        dev.off()

        # Plot initial vs. mean detection distances, data truncated to MEAN detection distance < truncDist distance
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km[detDistMean_km<truncDist], breaks=25, main=paste("n = ", length(detDist0_km[detDistMean_km<truncDist]),sep=""))
        hist(detDistMean_km[detDistMean_km<truncDist],breaks=25, main=paste("n = ", length(detDistMean_km[detDistMean_km<truncDist]),sep=""))
        plot(detDist0_km[detDistMean_km<truncDist],detDistMean_km[detDistMean_km<truncDist], main="")
        abline(a=0,b=1)
        hist(distDiff[detDistMean_km<truncDist], breaks=50)
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries - trunc meanDist 4km.tiff"), type="tiff")
        dev.off()

        # Plot initial vs. mean detection distances, data truncated to INITIAL detection distance < truncDist distance
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km[detDist0_km<truncDist], breaks=25, main=paste("n = ", length(detDist0_km[detDist0_km<truncDist]),sep=""))
        hist(detDistMean_km[detDist0_km<truncDist],breaks=25, main=paste("n = ", length(detDistMean_km[detDist0_km<truncDist]),sep=""))
        plot(detDist0_km[detDist0_km<truncDist],detDistMean_km[detDist0_km<truncDist], main="")
        abline(a=0,b=1)
        hist(distDiff[detDist0_km<truncDist], breaks=50)
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries - trunc initialDist 4km.tiff"), type="tiff")
        dev.off()

        ### Repeat some of the above plots, just for detectors that were never duty cycled

        # Plot initial vs. mean detection distances, data truncated to MEAN detection distance < truncDist distance
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km[detDistMean_km<truncDist & numInts==max(numInts) ], breaks=25, xlab="Initial Distance (km)",
             main=paste("n = ", length(detDist0_km[detDistMean_km<truncDist & numInts==max(numInts)]),sep=""))
        hist(detDistMean_km[detDistMean_km<truncDist & numInts==max(numInts)],breaks=25, xlab="Mean Distance (km)",
             main=paste("n = ", length(detDistMean_km[detDistMean_km<truncDist & numInts==max(numInts)]),sep=""))
        plot(detDist0_km[detDistMean_km<truncDist & numInts==max(numInts)],detDistMean_km[detDistMean_km<truncDist & numInts==max(numInts)],
             xlab="Initial detection km", ylab="Mean detection km", main="")
        abline(a=0,b=1)
        hist(distDiff[detDistMean_km<truncDist & numInts==max(numInts)], breaks=50, xlab="Initial distance - Mean distance", main="")
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries - trunc meanDist 4km - noDutyCycle.tiff"), type="tiff")
        dev.off()

        # Plot initial vs. mean detection distances, data truncated to INITIAL detection distance < truncDist distance
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km[detDist0_km<truncDist & numInts==max(numInts) ], breaks=25, xlab="Initial Distance (km)",
             main=paste("n = ", length(detDist0_km[detDist0_km<truncDist & numInts==max(numInts)]),sep=""))
        hist(detDistMean_km[detDist0_km<truncDist & numInts==max(numInts)],breaks=25, xlab="Mean Distance (km)",
             main=paste("n = ", length(detDistMean_km[detDist0_km<truncDist & numInts==max(numInts)]),sep=""))
        plot(detDist0_km[detDist0_km<truncDist & numInts==max(numInts)],detDistMean_km[detDist0_km<truncDist & numInts==max(numInts)],
             xlab="Initial detection km", ylab="Mean detection km", main="")
        abline(a=0,b=1)
        hist(distDiff[detDist0_km<truncDist & numInts==max(numInts)], breaks=50, xlab="Initial distance - Mean distance", main="")
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries - trunc initDist 4km - noDutyCycle.tiff"), type="tiff")
        dev.off()

        ### Plot mean detection distance vs. some encounter rate metrics
        axisSize = 1.5
        labSize = 1.75
        fontType = 2 # bold
        #x11()
        tiff(filename = file.path(plotDir, "Encounter rate metrics vs mean distance.tif"), width=10, height=10,units="in",res=300)
        par(mfrow=c(2,2), mar=c(5,5,4,2))
        # Set 1: all data within truncation range
        plot(detDistMean_km[detDist0_km<truncDist], dpue[detDist0_km<truncDist], main="",
             xlab="Mean detection distance (km)", ylab="Detections per on-interval", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
        plot(detDistMean_km[detDist0_km<truncDist], lastInt[detDist0_km<truncDist], main="",
             xlab="Mean detection distance (km)", ylab="Last interval detected", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
        # Set 2: only for non duty-cycled detectors
        plot(detDistMean_km[detDist0_km<truncDist & numInts==max(numInts)], dpue[detDist0_km<truncDist & numInts==max(numInts)],
             main="", xlab="Mean detection distance (km)", ylab="Detections per on-interval", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
        plot(detDistMean_km[detDist0_km<truncDist & numInts==max(numInts)], lastInt[detDist0_km<truncDist & numInts==max(numInts)],
             main="", xlab="Mean detection distance (km)", ylab="Last interval detected", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
        mtext(paste("All data within truncation distance (n = ", length(detDistMean_km[detDist0_km<truncDist]), ")", sep=""), side=3, outer=T,cex=2, line= -3, font=2)
        mtext(paste("Only from detectors not duty-cycled (n = ", sum(detDist0_km<truncDist & numInts==max(numInts)), ")", sep=""), side=3, outer=T, line=-33,cex=2,font=2)
        #savePlot(filename = "Encounter rate metrics vs mean distance.tif", type="tiff")
        dev.off()

        # Plot initial detection distance vs. some encounter rate metrics
        x11()
        par(mfrow=c(2,2))
        # Set 1: all data within truncation range
        plot(detDist0_km[detDist0_km<truncDist], dpue[detDist0_km<truncDist], main=paste("n = ", length(detDist0_km[detDist0_km<truncDist]),sep=""),
             xlab="Initial detection distance", ylab="Detections per on-interval")
        plot(detDist0_km[detDist0_km<truncDist], lastInt[detDist0_km<truncDist], main="All data within truncation distance",
             xlab="Initial detection distance", ylab="Last interval detected")
        # Set 2: only for non duty-cycled detectors
        plot(detDist0_km[detDist0_km<truncDist & numInts==max(numInts)], dpue[detDist0_km<truncDist & numInts==max(numInts)],
             main=paste("n = ", sum(detDist0_km<truncDist & numInts==max(numInts)),sep=""), xlab="Initial detection distance", ylab="Detections per on-interval")
        plot(detDist0_km[detDist0_km<truncDist & numInts==max(numInts)], lastInt[detDist0_km<truncDist & numInts==max(numInts)],
             main="Only from detectors not duty-cycled", xlab="Initial detection distance", ylab="Last interval detected")
        savePlot(filename = file.path(plotDir, "Encounter rate metrics vs initial distance.tiff"), type="tiff")
        dev.off()

        # duration length (last interval detected) on non duty-cycled detectors
        mean(lastInt[detDist0_km<truncDist & numInts==max(numInts)])  # mean value for "last Interval", for non-duty cycled detectors, within truncation distance
        # 7.57 for 4k trunction
        length(lastInt[detDist0_km<truncDist & numInts==max(numInts)])  # sample size
        # n=44 for 4k truncation
    }
    bugsDataEH
}

erDataPrep <- function(es, snaps, meanDepth=1217, hpDepth=110, truncDist=4000, A=1057925) {
    truncAngle <- atan(truncDist / (meanDepth-hpDepth)) * 180 / pi
    es <- filter(es, !is.na(meanAngle),
               meanAngle <= truncAngle,
               nClicks >= 3)
    es$station <- as.character(es$station)
    nDets <- group_by(es, station) %>%
        summarise(nDets = n())
    nDets <- left_join(nDets, snaps, by='station')
    nDets <- filter(nDets,
                    !is.na(nSnaps),
                    nDets > 0)
    # ziph.angles1 = 180 - ziph.data.in$MeanAngle                         # new variable: express angles so that 0 = animal below hydrophones, 90 = animal at surface (level w phones)
    # recorder.idx = rep(1,nrow(ziph.data.in))                            # new variable: indicate the recorder type for each row using a numeric label
    # recorder.idx[ziph.data.in$Recorder=="SM3M"]=2
    # recorder.idx[ziph.data.in$Recorder=="ST4300"]=3
    # ziph.data = cbind(ziph.data.in, ziph.angles1, recorder.idx)         # new array: bind the input array with the two new variables above
    # # ziph.data2 = ziph.data[ziph.data$nClicks>=3,]                       # filter out detections based on fewer than 3 clicks (because detections with fewer clicks could not confidently be identified as Ziphius... this did not eliminate much data)
    # # trunc.angle = atan(W.trunc/(mean.depth-mean.arrayDepth)) *180/pi    # find truncation angle (in radians), given specified truncation distance above (W.trunc)
    # # ziph.data.trunc = ziph.data2[ziph.data2$ziph.angles1 <= trunc.angle & !is.na(ziph.data2$ziph.angles1),]     # filter out NAs detections beyond truncation angle
    # n.per.station = table(ziph.data.trunc$Station)                      # find the number of detections per station
    # detPerStation = cbind(stationNames, rep(0,n.sites))                 # new matrix:  column 1 = station number (1-21), column 2 = number detections per station
    # detPerStation[stationNames %in% names(n.per.station),2] = n.per.station
    # list(n.sites=n.sites, n=detPerStation[,2], A=A, k=k, stationNames=stationNames)
    list(nSites = nrow(nDets),
         n = nDets$nDets,
         A=A,
         k=nDets$nSnaps
    )
    list(n.sites = nrow(nDets),
         n = nDets$nDets,
         A=A,
         k=nDets$nSnaps,
         stationNames = as.numeric(nDets$station)
    )
}