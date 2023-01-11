# This program reads the click files with merged angle data from binaries and finds modal angles for each event
#   and averages all the angles using tangent-weighting
# NOTE: this version uses data from PamGuard v1 and (for some drifts) PamGuard v2
#   it uses the pooled Event Click file from PatchTogetherAllEventClickFiles.r

# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
setwd("C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/Files for Submission/Supplemental R-code & data")

  ClicksWithAngles= read.csv(file="AllEventClicks wRC.csv")
  ClicksWithAngles$ClickDateTime= as.character(ClicksWithAngles$ClickDateTime)
  ClicksWithAngles$BinaryFile= basename(as.character(ClicksWithAngles$BinaryFile))
  ClicksWithAngles$EventId= as.character(ClicksWithAngles$EventId)
# get file start date/time from BinaryFile name
  Binary_File_Date_Time= sub(pattern="Click_Detector_Click_Detector_Clicks_",replacement="",x=ClicksWithAngles$BinaryFile)   #get UTC from file start time
  Binary_File_Date_Time= sub(pattern=".pgdf",replacement="",x=Binary_File_Date_Time)                                                #get UTC from file start time
  BinaryDateTime= strptime(Binary_File_Date_Time, "%Y%m%d_%H%M%S", tz="gmt")
  #correct Soundtrap file time/dates to get UTC
  Station= ClicksWithAngles$Station
  BinaryDateTime[Station == 1]=  BinaryDateTime[Station == 1] + 10.8*60*60             #add 10.8 hrs to correct filename to UTC
  BinaryDateTime[(Station<=6 |Station==9 |Station==12 |Station==15 |Station==16 | Station>=26)]=  
    BinaryDateTime[(Station<=6 |Station==9 |Station==12 |Station==15 |Station==16 | Station>=26)] + 7.0*60*60              #add 7.0 hrs to correct filename to UTC
  ClicksWithAngles$FileDateTime= BinaryDateTime

# Identify 2-min files with EventClicks using their binary file name
    UniqueFiles= unique(ClicksWithAngles$BinaryFile)
    CharUniqueFiles= trimws(as.character(UniqueFiles))
    n=0
    
    for (j in 1:length(UniqueFiles)) {
      SubSet1= ClicksWithAngles[ClicksWithAngles$BinaryFile==UniqueFiles[j],]
      UniqueEvents= unique(SubSet1$EventId)
      # to Guard against two files having the same binary name, also require EventID to be the same
      for (k in 1:length(UniqueEvents)) {
        SubSet2= SubSet1[SubSet1$EventId==UniqueEvents[k],]
        FileDateTime= SubSet2$FileDateTime[1]
        Angles= SubSet2$AngleRefractCorrected
        lat= SubSet2$lat[1]
        long= SubSet2$long[1]
        Station= SubSet2$Station
        nCountLT60sec= SubSet2$nCountLT60     #need to add "sec" after 60 on next run using most recent Merge program
        nCountGT60sec= SubSet2$nCountGT60     #need to add "sec" after 60 on next run using most recent Merge program
        ClickTime= strptime(as.character(SubSet2$ClickDateTime),format="%Y-%m-%d %H:%M:%S",tz="GMT")
        eventType= trimws(as.character(SubSet2$eventType[1]))
        HistOut=hist(Angles,breaks=seq(0,181,1),plot=FALSE)
        FindMaxFreq= which.max(HistOut$counts[c(1:88,92:181)])      #find modal value, exlc values around 90 deg.
        if (FindMaxFreq > 88) FindMaxFreq= FindMaxFreq + 3
        ModalAngle= HistOut$mids[FindMaxFreq]
        #Calculate mean detection angle for angles less than 88 deg if Modal Angles is less than 88 det
        MeanAngle= NA
        if (ModalAngle < 88) {  
          MeanAngle= atan(mean(tan(Angles[Angles<88]*pi/180),na.rm=TRUE)) * 180 / pi  #average angles is within 2 deg of mode
        }
        nClicks= length(SubSet2$Angle)
        Comment= as.character(SubSet2$comment)
        if (n == 0) {
          if (length(lat)>0) {
            EventAngles= data.frame(FileDateTime=FileDateTime,BinaryFile=CharUniqueFiles[j],EventId=as.character(UniqueEvents[k]),
                        eventType,Angle=Angles,ModalAngle=as.numeric(ModalAngle),MeanAngle=as.numeric(MeanAngle),
                        Station=as.numeric(Station),ClickTime=ClickTime,nClicks=as.numeric(nClicks),
                        comment=Comment,lat=lat,long=long,nCountLT60sec=nCountLT60sec,nCountGT60sec=nCountGT60sec)
            AllEventAngles= EventAngles
            n= n+1
          }
        } else {
          if (length(lat)>0) {
            EventAngles= data.frame(FileDateTime=FileDateTime,BinaryFile=CharUniqueFiles[j],EventId=UniqueEvents[k],eventType,Angle=Angles,ModalAngle=as.numeric(ModalAngle),
                        MeanAngle=as.numeric(MeanAngle),Station=as.numeric(Station),ClickTime=ClickTime,nClicks=as.numeric(nClicks),
                        comment=Comment,lat=lat,long=long,nCountLT60sec=nCountLT60sec,nCountGT60sec=nCountGT60sec)
            AllEventAngles= rbind(AllEventAngles,EventAngles) 
            n= n+1
          }
        }
      }
    } 
    
# sort AllEventAngles by Drift and DateTime
  AllEventAngles= AllEventAngles[order(AllEventAngles$Station,AllEventAngles$FileDateTime),]
# remove click specific variables and find only unique events
  AllEventAngles= subset(AllEventAngles,select=c(-Angle,-ClickTime))
  AllEventAngles= unique(AllEventAngles)
# fix some event type errors and use current codes
  AllEventAngles$eventType= as.character(AllEventAngles$eventType)
  AllEventAngles$eventType[AllEventAngles$eventType=="BW34-50+"]= "BW39V"
  AllEventAngles$eventType[AllEventAngles$eventType=="BW34-50"]= "BW39V"
  AllEventAngles$eventType[AllEventAngles$eventType=="BW46"]= "MS"
  AllEventAngles$eventType[AllEventAngles$eventType=="BW50-75.1"]= "MS"
  AllEventAngles$eventType[AllEventAngles$eventType=="BW50-75"]= "MS"
  AllEventAngles$eventType[AllEventAngles$eventType=="?BW"]= "BWunid"

  summary(as.factor(AllEventAngles$eventType))
    
  HistOut=hist(AllEventAngles$ModalAngle,breaks=seq(0,181,1),plot=FALSE)
  HistOut$counts
  HistOut=hist(AllEventAngles$MeanAngle,breaks=seq(0,181,1),plot=FALSE)
  HistOut$counts
  
  # add recorder type to dataframe
  DriftNum= AllEventAngles$Station
  AllEventAngles$Recorder= as.character('ST4300')
  AllEventAngles$Recorder[DriftNum == 7 |DriftNum == 10 |DriftNum == 13 |DriftNum == 17]=  'SM3M'
  AllEventAngles$Recorder[DriftNum == 8 |DriftNum == 11 |DriftNum == 14 |DriftNum == 18]=  'SM2Bat'
  AllEventAngles$Recorder= factor(AllEventAngles$Recorder)

  # Add new columns for datetime, station number, analyst, ICI and multiple angles (multiple animals)
  AllEventAngles$Analyst= as.factor(substr(AllEventAngles$comment,1,1))
  MultipleAngles= regexpr("MA",AllEventAngles$comment)
  AllEventAngles$MultipleAngles= TRUE
  AllEventAngles$MultipleAngles[MultipleAngles==-1]= FALSE
  AllEventAngles$MultipleAngles[is.na(MultipleAngles)]= NA
  SplitComments= strsplit(as.character(AllEventAngles$comment),split=",")
  ICI= rep(NA,length(SplitComments))
  for (i in 1:length(SplitComments)) {ICI[i]= SplitComments[[i]][7]}
  AllEventAngles$ICI= as.numeric(ICI)
  DateTime= as.POSIXct(as.character(AllEventAngles$FileDateTime), format=("%Y-%m-%d %H:%M:%S"), tz="GMT")
  n= length(AllEventAngles$Station) 
  AllEventAngles$DeltaMin= NA
  AllEventAngles$DeltaMin[2:n]= as.numeric(difftime(DateTime[2:n],DateTime[1:(n-1)],units="mins"))
  AllEventAngles$DeltaMin[AllEventAngles$DeltaMin<0]= NA
  
  write.csv(AllEventAngles,file="AllEventAngles.csv")
  AllEventAngles_drifts_1_22= AllEventAngles[AllEventAngles$Station <= 22,]
  write.csv(AllEventAngles_drifts_1_22,file="AvgSnapshotAngles.csv")
  
