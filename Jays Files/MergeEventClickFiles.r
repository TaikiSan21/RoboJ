# This program reads the click files with merged angle data from binaries and merges them into a single file
#   some drifts are from PamGuard v1 and some drifts are re-analyzed with v2_00_16e
# V2 drifts include #7, 17 & 22

# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
RootWD= "C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/Files for Submission/Supplemental R-code & data"

DataFolder1= paste(RootWD,"/BW EventClicks wAngles&Loc from PamGuard v1",sep="")
setwd(DataFolder1)
  filelist= list.files(pattern='*BW EventClicks_wID&Angle.csv')
  CardName= gsub(pattern='_MASTER-BW EventClicks_wID&Angle.csv',replacement = "",filelist)
  CardName= gsub(pattern='\\d',replacement = "",CardName)
  CardName= gsub(pattern='Station-_',replacement = "",CardName)
  CardName[substr(CardName,1,4)!="Card"]= NA
  
# append event click files from PamGuard v1 
  for (i in 1:length(filelist)) {
    ClicksWithAngles= read.csv(file=filelist[i])
    #convert eventId to character and add card number for SM recordings to avoid duplicate events within a drift
    ClicksWithAngles$EventId= as.character(ClicksWithAngles$EventId)  
    ClickDateTime= strptime(ClicksWithAngles$ClickDateTime,format="%Y-%m-%d %H:%M:%S", tz="gmt")
    ClicksWithAngles= ClicksWithAngles[order(ClickDateTime),]
    if (!is.na(CardName[i])) ClicksWithAngles$EventId= paste(ClicksWithAngles$EventId,CardName[i],sep="_")
    if (i == 1) {
      AllClicksWithAngles= ClicksWithAngles
    } else {
      AllClicksWithAngles= rbind(AllClicksWithAngles,ClicksWithAngles)
    }
  }
# Keep only a common set of variables variables so rbind will work
  AllClicksWithAngles2= subset(AllClicksWithAngles,select= c(BinaryFile,ClickNo,Angle,ClickDateTime,EventId,Amplitude,eventType,comment,lat,long,Station,nClicks,nCountLT60sec,nCountGT60sec))

# append data analyzed with PamGuard v2
  DataFolder2= paste(RootWD,"/BW EventClicks wAngles&Loc from PamGuard v2",sep="")
  setwd(DataFolder2)
# remove data from original analyses of drifts analyzed with PamGuard v2 
  filelist2= list.files(pattern='*AllEventClicks wAngle&LatLon.csv')
  DriftNumbers= as.numeric(gsub(".*-([0-9]+).*$", "\\1", filelist2))
  AllClicksWithAngles2= AllClicksWithAngles2[!(AllClicksWithAngles2$Station %in% DriftNumbers),]
# append PamGuard v2 data
  for (iDrift in 1:length(DriftNumbers)) {
    ClicksWithAngles= read.csv(file=filelist2[iDrift])
# Keep only a common set of variables variables so rbind will work
    ClickDateTime= strptime(ClicksWithAngles$ClickDateTime,format="%Y-%m-%d %H:%M:%S", tz="gmt")
    ClicksWithAngles= ClicksWithAngles[order(ClickDateTime),]
    ClicksWithAngles2= subset(ClicksWithAngles,select= c(BinaryFile,ClickNo,Angle,ClickDateTime,EventId,Amplitude,eventType,comment,lat,long,Station,nClicks,nCountLT60sec,nCountGT60sec))
    AllClicksWithAngles2= rbind(AllClicksWithAngles2,ClicksWithAngles2)
  }
  
# write all event click file
  setwd(RootWD)
  write.csv(AllClicksWithAngles2,file="AllEventClicks.csv")
  summary(AllClicksWithAngles2)
  levels(as.factor(AllClicksWithAngles2$EventId))
  levels(as.factor(AllClicksWithAngles2$eventType))
  levels(as.factor(AllClicksWithAngles2$Station))

  summary(as.factor(AllClicksWithAngles2$eventType))
  
    

