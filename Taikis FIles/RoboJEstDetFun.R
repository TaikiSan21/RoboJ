#  EstimateDetectionFunction from BWangles using maximum SimulatedLikelihood
#PASCAL Analyses for Barlow, Moore et al.
#"Acoustic detection range and population density of Cuvier's beaked whales estimated from
# near-surface hydrophones"

# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
# setwd("C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/Files for Submission/Supplemental R-code & data")
# TestMode= TRUE  #Will run much faster in TestMode, but results will be less precise

# LOAD REQUIRED LIBRARIES
#library(MASS)
library(mgcv)

#LOAD REQUIRED FUNCTIONS
source("../Data/RoboJ/Jays Files/SimulatedLikelihoodFunctions v7.R")


####################################################################################################
#LOAD SETUP PARAMETERS
####################################################################################################
estimateDetFunction <- function(ea, SpeciesID = 'ZC',
                                recorders = c("All Recorders","SM3M Recorders","SM2 Recorders","ST4300 Recorders"),
                                meanDepth = 1217,
                                sdDepth = 354,
                                HP_depth = 110,
                                truncDist = 4000,
                                SD_angleErr = 0.59,
                                DepthDistr = 'log-normal',
                                Range = seq(1000, 3000, 10),
                                JK_nsamples = 20,
                                Model = 'C_HN',
                                UseSimulatedData = FALSE,
                                plot=FALSE,
                                TestMode=TRUE,
                                doSink=FALSE,
                                gamk=25) {
    if(plot) {
        op <- par(xaxs="i",yaxs="i",cex.axis=1.3,cex.lab=1.5,font=2,font.lab=2)
        on.exit(par(op))
    }

    TruncAngle<- atan(truncDist/(meanDepth-HP_depth))*180/pi

    if (!UseSimulatedData) {     #get click angle data and subset it
        # inputFile<- "AllEventAngles_drifts_1-22.csv"
        # Use angles observed from PASCAL DASBRs using one mean angle per 2-minute sample
        # ClickAngles<- read.csv(file=inputFile)
        # Limit sample to 2-min snapshots with at least 3 clicks
        ea<- ea[(ea$nClicks>2),]
        # BeakedWhales<- ea[ea$eventType %in% c("BB","BW","BW26-47","BW38",
        #                                                         "BW39V","BW43","BWunid","MS","ZC"),]
        # cat(" Tot number of files w/ all beaked whale species and 3 or more clicks = ",length(ea$EventId))
        # Limit sample to SpeciesID
        ea<- ea[ea$eventType == SpeciesID,]
    }

    if (Model=="HN") {MaxPar2index<-1}else{MaxPar2index<-100}    #number of discrete values of Param2
    if (TestMode) {
        nSim<- 1e6
        jackknifeTF= FALSE # number of simulations used to estimate detection angle density distribution
    } else {
        nSim<- 1000000     # number of simulations used to estimate detection angle density distribution
        jackknifeTF<- TRUE   # estimate CV with jackknife based on drifts with detections
    }

    #Create folder to save detection function outputs
    outDir<- paste(getwd(),"/DetFuncEst/",Sys.Date(), sep="")
    dir.create(outDir, showWarnings = FALSE,recursive = TRUE)
    # setwd(outDir)

    # estimate detection function for each subset, outputs are redirected to a
    # separate file for each data subset
    outs <- vector('list', length=length(recorders))
    names(outs) <- recorders
    ####################################################################################################
    # get angle probability densities from simulation for each value of detection function parameters
    #   and save as matrix.  Values are average probabilities for integer angles from 0 to 90.
    #   This takes a while but only needs to be done once bcs values don't depend on ObsDetAngles
    ####################################################################################################
    cat('Making prob density matrix...\n')
    ProbDensityMatrix<- GetProbDensityData(nSim=nSim,Model=Model,Range=Range,MaxPar2index=MaxPar2index,
                                           DepthDistr=DepthDistr,meanDepth=meanDepth,sdDepth=sdDepth,SD_angleErr=SD_angleErr,
                                           TruncAngle=TruncAngle, HP_depth=HP_depth, truncDist=truncDist)
    for (SubSetName in recorders) {
        MainTitle<- paste("PASCAL",SubSetName,Model,sep=" ")
        # Send output to a file
        outputFile<- paste(Model,"DetFunctFit","n",nSim,truncDist,SubSetName,".txt",sep="_")
        outputFile <- file.path(outDir, outputFile)
        if(doSink) {
            sink(outputFile)
            on.exit(sink(file=NULL))
        }

        ####################################################################################################
        # get sample of angles from simulation or actual data
        ####################################################################################################
        if (UseSimulatedData) {
            nSample<- 500         #number of samples to use from simulation
            # Use sample of angles from simulation with known parameters
            if (Model == "HN") {
                Param1<- 2000
                Param2<- 0
            } else if (Model == "C_HN") {
                Param1<- 1200
                Param2<- 40
            }
            SimOut<- SimDetAngleDensity(nSim=nSim,param=c(Param1,Param2),Model=Model,DepthDistr=DepthDistr,
                                        meanDepth=meanDepth,sdDepth=sdDepth,SD_angleErr=SD_angleErr, MainTitle=MainTitle,
                                        TruncAngle=TruncAngle, HP_depth=HP_depth, truncDist=truncDist)
            nOut<- length(SimOut$Angle)                            #number of angles output from simulation
            DetAngleOut<- SimOut$Angle[runif(nOut)<SimOut$probDet] #use rejection method to select angles
            if (length(DetAngleOut) >= nSample) {ObsDetAngles<- DetAngleOut[1:nSample]} else {stop}
            hist(ObsDetAngles,xlim=c(0,100),breaks=seq(0,100,5))
            EDRactual<- SimOut$EDR
            cat(" Actual Param1=",Param1," Param2=",Param2," and EDR= :",EDRactual,"\n")
        } else {
            # Limit sample to DASBRs subsets based on recorder type
            if (SubSetName == "SM3M Recorders") {
                thisEa<- ea[ea$Recorder == "SM3M",]
            } else if (SubSetName == "SM2 Recorders") {
                thisEa<- ea[ea$Recorder == "SM2Bat",]
            } else if (SubSetName == "ST4300 Recorders") {
                thisEa<- ea[ea$Recorder == "ST4300",]
            } else {thisEa<- ea}

            DriftsWithDetection<- unique(thisEa$Station)
            nDrifts<- length(DriftsWithDetection)
            # Observed detection angles are valid average angles
            ObsDetAngles<- thisEa$MeanAngle[!is.na(thisEa$MeanAngle)]
            cat(" Tot number of files w/ ",SpeciesID," and 3 or more clicks = ",length(ObsDetAngles))

            # plot distribution of detection angles before truncation
            if (SubSetName == "All Recorders") { upperFreq<-80 }else{upperFreq<-40}
            histout<-hist(thisEa$MeanAngle,breaks=50,col="black",xlab="Detection Angle (deg)",
                          xlim=c(0,100),ylim=c(0,upperFreq),main=MainTitle)
            lines(c(TruncAngle,TruncAngle),c(0,80),lty="dashed",col="white",lwd=2)
            maxcount<- histout$breaks[which.max(histout$counts)]
            HdistAtMaxCount<- (meanDepth-HP_depth) * tan(maxcount*pi/180)
            cat(" Maximum count at =",maxcount,"deg, corresponding to a Hdist of ",HdistAtMaxCount,"m \n")
            # truncate observed detection angles
            PctTrunc<- 100*sum(ObsDetAngles>TruncAngle,na.rm=T)/sum(!is.na(ObsDetAngles))
            cat(" Percentage of detection angles truncated =",PctTrunc,"% \n")
            ObsDetAngles<- ObsDetAngles[ObsDetAngles<TruncAngle & !is.na(ObsDetAngles)]
            nSample<- length(ObsDetAngles)
        }


        ####################################################################################################
        # get likelihood surface data for observed angles by simulated likelihood
        ####################################################################################################
        cat('Making likelihood surface...\n')
        LogLikeDF<- GetLikelihoodData(ObsDetAngles=ObsDetAngles,Model=Model,
                                      ProbDensityMatrix=ProbDensityMatrix,Range=Range,MaxPar2index=MaxPar2index)

        # find maximum likelihood value from discrete values of detection function parameters
        MaxLike<- which.max(LogLikeDF$LogLike)
        cat(" Maximum log likelihood from likelihood surface data is:",LogLikeDF$LogLike[MaxLike],
            " with Param1=",LogLikeDF$Param1[MaxLike],"and Param2=",LogLikeDF$Param2[MaxLike],"\n")

        ####################################################################################################
        # fit likelihood surface data with a GAM smooth and save GAM object as global variable
        ####################################################################################################
        Param1Limits<- c(min(LogLikeDF$Param1),max(LogLikeDF$Param1)) #limits to avoid (outside range)
        Param2Limits<- c(min(LogLikeDF$Param2),max(LogLikeDF$Param2)) #limits to avoid (outside range)
        cat('Fitting gam...\n')
        if (Model == "HN") {          #univariate fit
            gam_logLike<- gam(formula= LogLike ~ te(Param1,k=gamk), data=LogLikeDF)
            plot(Range,gam_logLike$fitted.values,type="l",ylab="Log-Likelihood")
        }else{                        #bivariate fit
            gam_logLike<- gam(formula= LogLike ~ te(Param1,Param2,k=gamk), data=LogLikeDF)
            plot(gam_logLike)
        }
        gam.check(gam_logLike)
        SmoothLogLikelihood<- function(par) {
            Param1<- par[1]
            if (length(par) == 1) {Param2=0} else {Param2= par[2]}
            if ((Param1 < Param1Limits[1])|(Param1 > Param1Limits[2])|
                (Param2 < Param2Limits[1])|(Param2 > Param2Limits[2])) {
                LogL<- -exp(100)         # penalize LogL if parameters are outside fitted range of gam_logLike
            } else {
                LogL<- predict.gam(gam_logLike,newdata=data.frame(Param1=Param1,Param2=Param2))
            }
            return(LogL)
        }


        ####################################################################################################
        # use optim to find maximum likelihood parameter estimates (and Standard Errors)
        #           using optim and smoothed likelihood surface from gam
        ####################################################################################################
        cat('Doing optim...\n')
        if (Model == "C_HN") {
            ParGuess<- c(LogLikeDF$Param1[MaxLike],LogLikeDF$Param2[MaxLike])
            FindMaxLik<- optim(par=ParGuess,fn=SmoothLogLikelihood,method="Nelder-Mead",hessian= TRUE,
                               control=list(trace=ifelse(doSink, 5, 1),fnscale=-1))
        } else if (Model == "HN") {
            ParGuess<- 1500           #  half normal probabilities of detection
            FindMaxLik<- optim(par=ParGuess,fn=SmoothLogLikelihood,method="Brent",upper=4000,lower=500,
                               hessian= TRUE,control=list(trace=ifelse(doSink, 5, 1),fnscale=-1))
        }
        ####################################################################################################
        # covariance matrix and uncertainties for parameter estimates (as been deleted due to singular problem)
        # fisher_info<- solve(-FindMaxLik$hessian)
        # par_sigma<- sqrt(diag(fisher_info))
        cat(" Parameter estimates (Param-1 & Param-2) are:",FindMaxLik$par,"\n")
        # cat(" Standard errors in parameter estimates (Param-1 & Param-2) are:",par_sigma,"\n")


        ####################################################################################################
        #  Output results from maximum likelihood fitting
        ####################################################################################################
        # use simulation to get predicted distribution of detection probability as function of range
        SimOutput<- SimDetAngleDensity(nSim=nSim,param=c(FindMaxLik$par[1],FindMaxLik$par[2]),
                                       Model=Model,DepthDistr=DepthDistr,meanDepth=meanDepth,sdDepth=sdDepth,SD_angleErr=SD_angleErr,
                                       TruncAngle=TruncAngle, HP_depth=HP_depth, truncDist=truncDist)
        maxL_EDR<- SimOutput$EDR
        maxL_EDR39X<- SimOutput$EDR_39X
        plot((1:100)*80,SimOutput$probDetHRange,type="l",lwd=5,main=MainTitle,xlab="Horizontal Range (m)",
             ylab="Probability of Detection",ylim=c(0,1.05),xlim=c(0,6000))
        plot((1:100)*80/1000,SimOutput$probDetHRange,type="l",lwd=5,main=NULL,xlab="Horizontal Range (km)",
             ylab="Probability of Detection",ylim=c(0,1.05),xlim=c(0,5))
        HRangeDetProbDF<- data.frame(HRange=(1:100)*80,DetProb=SimOutput$probDetHRange)
        write.csv(HRangeDetProbDF,paste("HRangeDetProb_",Model,"_",SubSetName,".csv",sep=""))
        # plot cumulative distribution of angles from observations and simulation
        CDFout<- PlotCDF(ObsDetAngles=ObsDetAngles,param=c(FindMaxLik$par[1],FindMaxLik$par[2]),Model=Model,
                         DepthDistr=DepthDistr,meanDepth=meanDepth,sdDepth=sdDepth,SD_angleErr=SD_angleErr,
                         TruncAngle=TruncAngle, HP_depth=HP_depth, truncDist=truncDist)

        maxL_ESA<- pi*maxL_EDR^2  #effective survey area
        cat(" Sample size of detection angles (#snapshots):",nSample,"\n")
        cat(" Detection function for: ",MainTitle,"\n")
        cat(" Maximum likelihood estimate of EDR:",maxL_EDR/1000,"km \n")
        cat(" Maximum likelihood estimate of ESA:",maxL_ESA/(1000^2),"km2 \n")
        cat(" ML goodness of fit from KS stat, p=",CDFout$KS_prob,"\n")
        cat(" The log-likelihood of ML estimates is",FindMaxLik$value,"\n")
        cat(" AIC =",-2*(FindMaxLik$value-length(FindMaxLik$par)),"\n")
        cat(" Prob detection at zero horiz range (gzero)=",HRangeDetProbDF$DetProb[1])
        if(doSink) {
            save.image(paste(Model,"DetFunctFit","n=",nSim,truncDist,SubSetName,".RData",sep="_"))
        }
        thisOut <- list(MaxLikeParam=FindMaxLik$par,nSample=nSample, MainTitle=MainTitle, maxL_EDR=maxL_EDR,
                        maxL_ESA=maxL_ESA, KS_prob=CDFout$KS_prob, MaxLik=FindMaxLik$value,
                        AIC = -2*(FindMaxLik$value-length(FindMaxLik$par)),
                        gzero = HRangeDetProbDF$DetProb[1],
                        HRangeDetProbDF=HRangeDetProbDF,
                        lld = LogLikeDF,
                        gam = gam_logLike)
        ###########################################################################
        # BEGIN JACKKNIFE SAMPLING
        ###########################################################################
        if (jackknifeTF) {
            JK_EDR<- array(1)
            JK_EDR_39X<- array(1)
            nAngles<- length(ObsDetAngles)
            JK_sample<- ceiling((1:nAngles) * JK_nsamples / nAngles)

            for (j in 1:JK_nsamples) {
                JK_ObsDetAngles<- ObsDetAngles[!(JK_sample == j)]

                ##########################################################################
                # get likelihood surface data for observed angles by simulated likelihood
                ##########################################################################
                LogLikeDF<- GetLikelihoodData(ObsDetAngles=JK_ObsDetAngles,Model=Model,
                                              ProbDensityMatrix=ProbDensityMatrix,Range=Range,MaxPar2index=MaxPar2index)
                # find maximum likelihood value from discrete values of detection function parameters
                MaxLike<- which.max(LogLikeDF$LogLike)
                cat(" Jackknife sample:",j,"\n")
                cat(" Sample size of detection angles (#snapshots):",length(JK_ObsDetAngles),"\n")
                cat(" Maximum log likelihood from likelihood surface data is:",LogLikeDF$LogLike[MaxLike],
                    " with Param1=",LogLikeDF$Param1[MaxLike],"and Param2=",LogLikeDF$Param2[MaxLike],"\n")


                ##########################################################################
                # fit likelihood surface data with a GAM smooth and save GAM object as global variable
                ##########################################################################
                if (Model == "HN") {          #univariate fit
                    gam_logLike<- gam(formula= LogLike ~ te(Param1,k=gamk), data=LogLikeDF)
                }else{                        #bivariate fit
                    gam_logLike<- gam(formula= LogLike ~ te(Param1,Param2,k=gamk), data=LogLikeDF)
                }

                ################################################################################################
                # use optim to find maximum likelihood parameter estimates (and Standard Errors)
                #           using optim and smoothed likelihood surface from gam
                ################################################################################################
                if (Model == "C_HN") {
                    ParGuess<- c(LogLikeDF$Param1[MaxLike],LogLikeDF$Param2[MaxLike])      #  compound half normal probabilities of detection
                    FindMaxLik<- optim(par=ParGuess,fn=SmoothLogLikelihood,method="Nelder-Mead",hessian= TRUE,
                                       control=list(trace=ifelse(doSink, 5, 1),fnscale=-1))
                } else if (Model == "HN") {
                    ParGuess<- 1500           #  half normal probabilities of detection
                    FindMaxLik<- optim(par=ParGuess,fn=SmoothLogLikelihood,method="Brent",upper=4000,lower=500,
                                       hessian= TRUE,control=list(trace=ifelse(doSink, 5, 1),fnscale=-1))
                }
                SimOutput<- SimDetAngleDensity(nSim=1000000,c(FindMaxLik$par[1],FindMaxLik$par[2]),Model=Model,
                                               DepthDistr=DepthDistr,meanDepth=meanDepth,sdDepth=sdDepth,
                                               SD_angleErr=SD_angleErr,
                                               TruncAngle=TruncAngle, HP_depth=HP_depth, truncDist=truncDist)
                JK_EDR[j]<- SimOutput$EDR
                JK_EDR_39X[j]<- SimOutput$EDR_39X

                cat(" Maximum likelihood estimate of EDR:",JK_EDR[j],"\n")

            }


            ################################################################################
            # Summary output for snapshot (2-min) detection probability
            JK_EDR
            cat(" Maximum likelihood EDR (1X)=",maxL_EDR,"m \n")
            nJK<- length(JK_EDR)
            mean_JK_EDR<- mean(JK_EDR)
            thisOut$mean_JK_EDR <- mean_JK_EDR
            JackKnife_var<- ((nJK-1)/nJK) * sum((JK_EDR-mean_JK_EDR)^2)
            thisOut$JackKnife_var <- JackKnife_var
            thisOut$JackKnife_cv <- sqrt(JackKnife_var)/mean_JK_EDR
            cat(" Jackknife mean EDR=",mean_JK_EDR,"  se =",sqrt(JackKnife_var),"  cv =",
                sqrt(JackKnife_var)/mean_JK_EDR,"\n")

            JK_ESA<- pi*JK_EDR^2  #effective survey area
            mean_JK_ESA<- mean(JK_ESA)
            thisOut$mean_JK_ESA <- mean_JK_ESA
            Jackknife_varESA<- ((nJK-1)/nJK) * sum((JK_ESA-mean_JK_ESA)^2)
            thisOut$Jackknife_varESA <- Jackknife_varESA
            thisOut$JackKnife_cvESA <- sqrt(Jackknife_varESA)/mean_JK_ESA
            cat(" Jackknife mean ESA=",mean_JK_ESA,"  se =",sqrt(Jackknife_varESA),"  cv =",
                sqrt(Jackknife_varESA)/mean_JK_ESA, "\n")

            ################################################################################
            # Summary output for dive detection probability (estimated a 39 independent snapshots)
            cat(" Maximum likelihood dive EDR (39X)=",maxL_EDR39X,"m \n")
            JK_EDR_39X
            nJK<- length(JK_EDR_39X)
            mean_JK_EDR_39X<- mean(JK_EDR_39X)
            thisOut$mean_JK_EDR_39X <- mean_JK_EDR_39X
            JackKnife_var<- ((nJK-1)/nJK) * sum((JK_EDR_39X-mean_JK_EDR_39X)^2)
            thisOut$JackKnife_var_39x <- JackKnife_var
            thisOut$JackKnife_cv_39x <- sqrt(JackKnife_var)/mean_JK_EDR_39X
            cat(" Jackknife mean dive EDR (39X)=",mean_JK_EDR_39X,"  se =",sqrt(JackKnife_var),"  cv =",
                sqrt(JackKnife_var)/mean_JK_EDR_39X,"\n")

            JK_ESA_39X<- pi*JK_EDR_39X^2  #effective survey area
            mean_JK_ESA_39X<- mean(JK_ESA_39X)
            thisOut$mean_JK_ESA_39X <- mean_JK_ESA_39X
            Jackknife_varESA_39X<- ((nJK-1)/nJK) * sum((JK_ESA_39X-mean_JK_ESA_39X)^2)
            thisOut$Jackknife_varESA_39X <- Jackknife_varESA_39X
            thisOut$JackKnife_cvESA_39x <- sqrt(Jackknife_varESA_39X)/mean_JK_ESA_39X
            cat(" Jackknife mean dive ESA (39X)=",mean_JK_ESA_39X,"  se =",sqrt(Jackknife_varESA_39X),"  cv =",
                sqrt(Jackknife_varESA_39X)/mean_JK_ESA_39X, "\n")
            if(doSink) {
                save.image(paste(Model,"DetFunctFit","n=",nSim,truncDist,SubSetName,"JK.RData",sep="_"))
            }
        }
        if(doSink) {
            sink()
            file.show(outputFile)
        }

        outs[[SubSetName]] <- thisOut
    }  # next subset
    outs
}

