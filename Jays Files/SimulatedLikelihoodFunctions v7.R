# filename= "Required Functions for Max SimulatedLikelihood v5.R"
#  This script included required functions for:

#  EstimateDetectionFunction from BWangles using maximum SimulatedLikelihood
#Barlow, Fregosi, Thomas, Harris and Griffiths
#"Acoustic detection range and population density of Cuvier's beaked whales estimated from
# near-surface hydrophones"
#Section IID2


###########################################################################
# compound half-normal detection function (is halfnormal when ShapeParam=0)
###########################################################################
CHN<- function(range,ScaleParam,ShapeParam) {
    prob<- 1- (1-exp(-range^2/(2*ScaleParam^2)))^(exp(ShapeParam))
}
##########################################################################


###################################################################################################
# Function estimates probability density for detection angles from simulation
###################################################################################################
SimDetAngleDensity<- function(nSim=1e4,param=NULL,Model=c('HN', 'C_HN'),DepthDistr=c('log-normal', 'normal'),
                              meanDepth=1217,sdDepth=354,SD_angleErr=0.59,
                              MainTitle='', TruncAngle=NULL, HP_depth=110, truncDist=4e3, fast=NULL, adOnly=FALSE) {
    Model <- match.arg(Model)
    if(is.null(fast)) {
        fast <- nSim > 3e4
    }
    if(is.null(param)) {
        if (Model == "HN") {
            param <- c(2000, 0)
        } else if (Model == "C_HN") {
            param <- c(1200, 40)
        }
    }
    MaxHRange<- 8000                          # maximum horizontal range (m)
    DepthDistr <- match.arg(DepthDistr)
    if(is.null(TruncAngle)) {
        TruncAngle<- atan(truncDist/(meanDepth-HP_depth))*180/pi
    }
    if(fast) {
        # set.seed(156)
        df <- data.frame(HRange = sqrt((MaxHRange^2)*runif(nSim)))
        df$Depth <- if(DepthDistr == "log-normal") {
            # set.seed(256)
            lognorm(nSim,mean=(meanDepth-HP_depth),sd=sdDepth)   # distribution depth (m)
        } else if (DepthDistr == "normal") {
            rnorm(nSim,mean=(meanDepth-HP_depth),sd=sdDepth)     # distribution depth (m)
        } else {
            cat("  ERROR, need to specify the depth distribution function as normal or log-normal")
        }
        df$Depth[df$Depth < 300] <- 300
        df <- mutate(df,
                     SRange = sqrt(HRange^2 + Depth^2),
                     # probDet = case_when(Model == 'C_HN' ~ CHN(SRange, param[1], param[2]),
                     #                  Model == 'HN' ~CHN(SRange, param[1], 0)),
                     Angle = abs(rnorm(nSim, mean = atan(HRange/Depth)*180/pi, sd=SD_angleErr))) %>%
            filter(Angle < TruncAngle) %>%
            mutate(HRindex = ceiling(100*HRange/MaxHRange))
        df$probDet <- if(Model == 'C_HN') {
            CHN(df$SRange, param[1], param[2])
        } else if(Model == 'HN') {
            CHN(df$SRange, param[1], 0)
        }
        # browser()
        # df$SRange <- sqrt(df$HRange^2 + df$Depth^2)
        # df$probDet <- if(Model == 'C_HN') {
        #     CHN(df$SRange, param[1], param[2])
        # } else if(Model == 'HN') {
        #     CHN(df$SRange, param[1], 0)
        # }
        # df$Angle <- atan(df$HRange/df$Depth) * 180 / pi
        # set.seed(356)
        # df$Angle <- abs(rnorm(nSim, mean = df$Angle, sd=SD_angleErr))
        # df <- df[df$Angle < TruncAngle,]
        # probDetHRange<- rep(1,100)
        sumDf <- data.frame(iAngle = c(ceiling(df$Angle), 1:91),
                            probDet = c(df$probDet, rep(0, 91)))
        sumprob <- group_by(sumDf, iAngle) %>% summarise(sum=sum(probDet)) %>% .$sum
        AngleDensity<- sumprob/sum(sumprob) #empirical prob. density function for angles from simulation
        if(adOnly) {
            return(list(AngleDensity=AngleDensity))
        }
        dx<- MaxHRange/100
        # df$HRindex<- ceiling(100*df$HRange/MaxHRange)
        HR= (1:100 - 0.5)*dx
        probDetDf <- group_by(df, HRindex) %>% summarise(pdhr=mean(probDet, na.rm=TRUE))
        probDetHRange <- probDetDf$pdhr
        if(length(probDetHRange) < 100) {
            new <- rep(NA, 100)
            new[1:100 %in% probDetDf$HRindex] <- probDetDf$pdhr
            probDetHRange <- new
        }
        # browser()
        # probDetHRange[probDetHRange$HRindex] <- probDetHRange$pdhr
        IntegralRangeTimesDetProb<- sum(HR[HR<=truncDist] * probDetHRange[HR<=truncDist],na.rm = TRUE)*dx
        EDR<- sqrt(2*IntegralRangeTimesDetProb)
        # calculate the EDR for 39 independent 2-minute detection opportunities for dive
        probDet39X= 1 - (1-probDetHRange)^39
        IntegralRangeTimesDetProb<- sum(((1:100 - 0.5)*dx) * probDet39X,na.rm = TRUE)*dx
        EDR_39X<- sqrt(2*IntegralRangeTimesDetProb)
        # normalized sum of det. prob. for empirical prob. density for integer angles from simulation
        # iAngle<- ceiling(df$Angle)            #angle index
        # probDet2<- c(df$probDet,rep(0,91))    #add zeros to make sure all angles are included in "by"
        # iAngle2<- c(iAngle,1:91)           #add corresponding angles
        # sumprob<- by(probDet2,iAngle2,sum) #add detection probabilities for each angle
        # sumprob<- as.numeric(sumprob)
        # sumDf <- data.frame(iAngle = c(ceiling(df$Angle), 1:91),
        #                     probDet = c(df$probDet, rep(0, 91)))
        # sumprob <- group_by(sumDf, iAngle) %>% summarise(sum=sum(probDet)) %>% .$sum
        # AngleDensity<- sumprob/sum(sumprob) #empirical prob. density function for angles from simulation
        # output values calculated in simulation
        # browser()
        OutList<- list(AngleDensity=AngleDensity,probDetHRange=probDetHRange,EDR=EDR,Angle=df$Angle,
                       probDet=df$probDet,EDR_39X=EDR_39X)
    } else {
        # set.seed(156)
        HRange<- sqrt((MaxHRange^2)*runif(nSim))    # expected triangular distribution for horiz range
        if (DepthDistr == "log-normal") {
            # set.seed(256)
            Depth<- lognorm(nSim,mean=(meanDepth-HP_depth),sd=sdDepth)   # distribution depth (m)
        } else if (DepthDistr == "normal") {
            Depth<- rnorm(nSim,mean=(meanDepth-HP_depth),sd=sdDepth)     # distribution depth (m)
        } else {
            cat("  ERROR, need to specify the depth distribution function as normal or log-normal")
        }

        Depth[Depth<300]<-300                          # truncate depth at biologically reasonable values
        SRange<- sqrt(HRange^2 + Depth^2)     # infered distribution for slant ranges
        #specify functional form of the detection function
        if (Model == "C_HN") {
            probDet<- CHN(SRange,param[1],param[2]) #  compound half normal probabilities of detection
        } else if (Model == "HN") {
            probDet<- CHN(SRange,param[1],0)        #  half normal probabilities of detection
        }
        # calculate distribution of angles
        Angle<- atan(HRange/Depth) * 180/pi      # simulated distribution of angles in degrees
        Angle<- abs(rnorm(nSim,mean=Angle,sd=SD_angleErr))  # add random angle error
        # set.seed(356)
        # truncate samples based on detection angle
        angleFilt <- Angle < TruncAngle

        probDet<- probDet[angleFilt]
        Angle<-   Angle[angleFilt]

        probDet2<- c(probDet,rep(0,91))    #add zeros to make sure all angles are included in "by"
        iAngle2<- factor(c(ceiling(Angle),1:91), levels=1:91)         #add corresponding angles
        # sumprob<- by(probDet2,iAngle2,sum) #add detection probabilities for each angle
        # sumprob<- as.numeric(sumprob)
        sumprob <- sapply(split(probDet2, iAngle2), sum)
        AngleDensity<- sumprob/sum(sumprob) #empirical prob. density function for angles from simulation
        if(adOnly) {
            return(list(AngleDensity=AngleDensity))
        }
        HRange<-  HRange[angleFilt]
        Depth<-   Depth[angleFilt]
        # calculate average prob of detection in 100 categories of horizontal range
        probDetHRange<- rep(NA,100)
        dx<- MaxHRange/100
        HRindex<- ceiling(100*HRange/MaxHRange)
        HR= (1:100 - 0.5)*dx
        # for (index in 1:100) {
        #     probDetHRange[index]<- mean(probDet[HRindex==index], na.rm=TRUE)
        # }
        HRindex <- factor(HRindex, levels=1:100)
        probDetHRange <- sapply(split(probDet, HRindex), function(p) mean(p, na.rm=TRUE))
        # if(length(pdr) < 100) {
        #     probDetHRange[as.numeric(names(pdr))] <- unname(pdr)
        # } else {
        #     probDetHRange <- unname(pdr)
        # }
        # browser()
        IntegralRangeTimesDetProb<- sum(HR[HR<=truncDist] * probDetHRange[HR<=truncDist],na.rm = TRUE)*dx
        EDR<- sqrt(2*IntegralRangeTimesDetProb)
        # calculate the EDR for 39 independent 2-minute detection opportunities for dive
        probDet39X= 1 - (1-probDetHRange)^39
        IntegralRangeTimesDetProb<- sum(((1:100 - 0.5)*dx) * probDet39X,na.rm = TRUE)*dx
        EDR_39X<- sqrt(2*IntegralRangeTimesDetProb)
        # normalized sum of det. prob. for empirical prob. density for integer angles from simulation
        # iAngle<- ceiling(Angle)            #angle index

        # output values calculated in simulation
        # browser()
        OutList<- list(AngleDensity=AngleDensity,probDetHRange=probDetHRange,EDR=EDR,Angle=Angle,
                       probDet=probDet,EDR_39X=EDR_39X)
    }
    # browser()
    # calculate the EDR from the integral (sum) of range time detection probability within trunc distance

    return(OutList)
}
################################################################################################


#######################################################################################
# Use simulation to estimate angle probability density at discrete values
#  of the scale (Param1) & shape (Param2) parameters and save as a dataframe.
#######################################################################################
GetProbDensityData<- function(nSim=1e4,Model='C_HN',Range=seq(1000, 3000, 10),MaxPar2index=NULL,
                              DepthDistr='log-normal',meanDepth=1217,sdDepth=354,SD_angleErr=0.59,
                              TruncAngle=NULL, HP_depth=110, truncDist=4e3, fast=NULL) {

    # MaxPar1index<- max(Range)
    if(is.null(fast)) {
        fast <- nSim > 3e4
    }
    if(is.null(MaxPar2index)) {
        MaxPar2index <- ifelse(Model == 'C_HN', 100, 1)
    }
    par2Range <- if(Model == 'C_HN') {
        (1:MaxPar2index) * 30 / MaxPar2index - 2
    } else{
        0
    }
    if(is.null(TruncAngle)) {
        TruncAngle<- atan(truncDist/(meanDepth-HP_depth))*180/pi
    }
    # ProbDensityMatrix<- array(NA,c(MaxPar1index,MaxPar2index,91))
    ProbDensityMatrix <- array(NA, c(length(Range), length(par2Range), 91))
    # pb <- txtProgressBar(min=0, max=length(Range) * length(par2Range), style=3)
    # ix <- 1
    for (ix1 in seq_along(Range)) {
        # cat(", Param1=",Param1)   #output first parameter to keep track of progress
        Param1 <- Range[ix1]
        # setTxtProgressBar(pb, value=ix)
        for (ix2 in seq_along(par2Range)) {
            # if (Model == "C_HN") {
            #     Param2<- Par2index*30/MaxPar2index - 2   #for chn, max of shape param is 30-2=28
            # }else if (Model == "HN") {
            #     Param2<- 0
            # }
            Param2 <- par2Range[ix2]
            #get simulated prob density for angles for given parameter estimates
            SimOut<- SimDetAngleDensity(nSim=nSim, param=c(Param1,Param2), Model=Model,
                                        DepthDistr=DepthDistr,meanDepth=meanDepth,sdDepth=sdDepth,SD_angleErr=SD_angleErr,
                                        TruncAngle=TruncAngle, HP_depth=HP_depth, truncDist=truncDist, fast=fast, adOnly=TRUE)
            ProbDensityAngle<- SimOut$AngleDensity        #get average prob density for integer angles
            ProbDensityAngle[ProbDensityAngle==0]<- 1/(10*nSim)    #replace zeros with very low value
            ProbDensityMatrix[ix1,ix2,]<- ProbDensityAngle
            # browser()
            # ix <- ix + 1
        }
    }
    return(ProbDensityMatrix)
}

#######################################################################################
# Calculate the sum of log-likelihood of observed angle distribution at discrete values
#  of the scale (Param1) & shape (Param2) parameters and save as a dataframe.
#######################################################################################
GetLikelihoodData<- function(ObsDetAngles=ObsDetAngles,Model=Model,
                             ProbDensityMatrix=ProbDensityMatrix,Range=Range,MaxPar2index=MaxPar2index) {
    par2Range <- if(Model == 'C_HN') {
        (1:MaxPar2index) * 30 / MaxPar2index - 2
    } else{
        0
    }
    LogLikeDF<- data.frame(LogLike=rep(NA,length(par2Range)*length(Range)))
    LogLikeDF$Param1<- NA
    LogLikeDF$Param2<- NA
    index<- 0
    for (ix1 in seq_along(Range)) {
        Param1 <- Range[ix1]
        for (ix2 in seq_along(par2Range)) {
            ProbDensityAngle<- ProbDensityMatrix[ix1,ix2,]
            # if (Model == "C_HN") {
            #     Param2<- Par2index*30/MaxPar2index - 2
            # }else if (Model == "HN") {
            #     Param2<- 0
            # }
            Param2 <- par2Range[ix2]
            index<- index +1
            LogLikeDF$LogLike[index]<- sum(log(ProbDensityAngle[ceiling(ObsDetAngles)]))
            LogLikeDF$Param1[index]<- Param1
            LogLikeDF$Param2[index]<- Param2
        }
    }
    LogLikeDF<- LogLikeDF[!is.na(LogLikeDF$LogLike),]
    return(LogLikeDF)
}



###########################################################################
# plot cumulative distributions of sample and expected from simulated likelihood expected
#    and associated KS statistics
###########################################################################
PlotCDF<- function(ObsDetAngles=ObsDetAngles,param=param,Model=Model,DepthDistr=DepthDistr,
                   meanDepth=meanDepth,sdDepth=sdDepth,SD_angleErr=SD_angleErr, MainTitle='',
                   TruncAngle, HP_depth, truncDist) {
    SimOut<- SimDetAngleDensity(nSim=1000000,param=param,Model=Model,DepthDistr=DepthDistr,
                                meanDepth=meanDepth,sdDepth=sdDepth,SD_angleErr=SD_angleErr,
                                MainTitle=MainTitle, TruncAngle=TruncAngle, HP_depth=HP_depth,
                                truncDist=truncDist)
    cdf_ObsDetAngless<- sapply(1:91,function(x) mean(ObsDetAngles<=x))
    inclAngle<- SimOut$Angle[runif(length(SimOut$probDet))<SimOut$probDet]
    cdf_Angles<- sapply(1:91,function(x) mean(inclAngle<=x))
    plot(1:91,cdf_Angles,type="l",lwd=3,main=MainTitle,xlab="Detection Angle (deg.)",
         ylab="Cumulative Distribution")
    lines(1:91,cdf_ObsDetAngless,col="red",lwd=3)
    KS_out<- ks.test(ObsDetAngles,inclAngle,exact=FALSE)  #test of observed and predicted angles
    OutList<- list(KS_out$statistic,KS_prob=KS_out$p.value)
    return(OutList)
}


###########################################################################
# Function for random log-normal distribution
###########################################################################
lognorm<- function(n,mean,sd) {
    meanlog <- log(mean^2 / sqrt(sd^2 + mean^2))
    sdlog<- sqrt(log(1 + (sd^2 / mean^2)))
    rlnorm(n,meanlog,sdlog)
}
###########################################################################


###################################################################################################
# evaluate the smoothed likelihood surface at given parameter values for the compound half-normal
#   called by function optim to find maximum likelihood
###################################################################################################
# SmoothLogLikelihood<- function(par) {
#     Param1<- par[1]
#     if (length(par) == 1) {Param2=0} else {Param2= par[2]}
#     if ((Param1 < Param1Limits[1])|(Param1 > Param1Limits[2])|
#         (Param2 < Param2Limits[1])|(Param2 > Param2Limits[2])) {
#         LogL<- -exp(100)         # penalize LogL if parameters are outside fitted range of gam_logLike
#     } else {
#         LogL<- predict.gam(gam_logLike,newdata=data.frame(Param1=Param1,Param2=Param2))
#     }
#     return(LogL)
# }
###################################################################################################
