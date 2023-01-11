#  Merge Pamguard DBs & Binaries & Locations v4.r  for PASCAL
#NOTE: this version is for use with Binary csv file from "ReadPamGuardBinaries v10 forPASCAL2016.R"
#  PamGuard files must be from version 2 (with unique identifiers: UID)

# Requires cleaned-up Spot location data from "CleanupPASCAL_SpotData.r"
#    to convert to UTC and eliminate outliers

# Truncate data based on deployment and retrieval times, and
# Add LatLon position to Events based on nearest times in Spot data, and
# merge with PamGuard DB event clicks (adding species ID), and 
# merge with PamGuard binary data to add click type and bearing angle
#setup
library(RSQLite)
library(PamBinaries)    #package 'PamBinaries' requires R >=  3.4.0
sqlite <- dbDriver("SQLite")
options(digits.sec=3)


##############################################################################################
# get deployment and retrieve times in Excel format and convert to R format
DeployRetrTimes= read.csv("C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/DateTimeDataReconciliation.csv")
DeployTime= as.POSIXct((DeployRetrTimes$UTC_DeployTime-2)*24*60*60,origin='1900-01-01',tz='gmt')
RetrTime= as.POSIXct((DeployRetrTimes$UTC_RetrTime-2)*24*60*60,origin='1900-01-01',tz='gmt')

##############################################################################################
# get cleaned up Spot location data
SpotData_wUTC= read.csv(file="C:\\Jay\\ACOUSTIC\\Buoy Recorder\\PASCAL\\AllSpotTracks wUTC.csv")
SpotData_wUTC$UTC= strptime(SpotData_wUTC$dateTime,format="%m/%d/%Y %H:%M",tz="gmt")

#####################################################################################################
# set working directory location of sqlite databases
#  setwd("E:\\PASCAL PAMGUARD Files\\Tests of Jaimes New BW Detector\\Drift-28")
# setwd("P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-17\\")
setwd("P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-22\\")
# setwd("P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-7\\")
# create list of all the SQLite database files within the folder
  files= list.files(pattern=glob2rx("Station*.sqlite3"),recursive=TRUE,full.names=FALSE)
  nfiles= length(files)
  DriftNum= array(NA)
  ifile= 1
  CardNames= dirname(files)
  for (ifile in 1:nfiles) {
  
    cat(ifile," of ",nfiles,"  ",files[ifile],"\n")
    conn <- dbConnect(sqlite,files[ifile])                       #connect to the database
    
    Events <- dbReadTable(conn, "Click_Detector_OfflineEvents")         #read offline events
    #convert UTC time info from Events into date/time and sort by date/time
    Events$UTC=  strptime(Events$UTC, "%Y-%m-%d %H:%M:%S", tz="gmt")  
    Events= Events[order(Events$UTC),]
    EventFileName= sub(pattern=".sqlite3",replacement="_Events.csv",x=files[ifile])
    write.csv(Events,file=EventFileName)                                #write events to csv file
    
    EventClicks <- dbReadTable(conn, "Click_Detector_OfflineClicks")         #read offline event clicks
    #convert UTC time info from Events into date/time and sort by date/time
    EventClicks$UTC=  strptime(EventClicks$UTC, "%Y-%m-%d %H:%M:%S", tz="gmt")   
    EventClicks= EventClicks[order(EventClicks$UTC),]
    ClickFileName= sub(pattern=".sqlite3",replacement="_EventClicks.csv",x=files[ifile])
    write.csv(EventClicks,file=ClickFileName)                                #write event clicks to csv file
    dbDisconnect(conn)

#determine Drift numbers from filenames
    DriftNum= as.character(as.numeric(gsub("[^0-9.]", "",  sub('\\.sqlite3$','',files[ifile]))))
    DriftNum[ifile]= as.numeric(DriftNum)
    NumberClicks= array(NA,length(files))
    HistSave= rep(0,24)
  
# for each file name in list, read event and click files from PamGuard database,
# adding eventType, nClicks, comments and Drift number to click files
    if (DriftNum[ifile]<=28) {    # note: drifts 29 and 30 are combined with other drifts
      Events= Events[,1:13]     #eliminate unused variables
#assign drift number for each event
      Events$Drift= DriftNum[ifile]
#correct Soundtrap file time/dates to get ClickDateTime
      if (DriftNum[ifile] == 1) {
        Events$ClickDateTime=  Events$UTC + 10.8*60*60             #add 10.8 hrs to correct filename to UTC
      } else if (DriftNum[ifile]<=6 |DriftNum[ifile]==9 |DriftNum[ifile]==12 |DriftNum[ifile]==15 |DriftNum[ifile]==16 | DriftNum[ifile]>=26) {
        Events$ClickDateTime=  Events$UTC + 7.0*60*60              #add 7.0 hrs to correct filename to UTC
      } else {
        Events$ClickDateTime=  Events$UTC
      }
#select spot data for this Drift
      StaSpotData= SpotData_wUTC[SpotData_wUTC$Drift==DriftNum[ifile],]
      if (DriftNum[ifile]==26) {
        StaSpotData= SpotData_wUTC[SpotData_wUTC$Drift==26 | SpotData_wUTC$Drift==29,]
      } else if (DriftNum[ifile]==27) {
        StaSpotData= SpotData_wUTC[SpotData_wUTC$Drift==27 | SpotData_wUTC$Drift==30,]    
      }
#find nearest match in Spot location
      n= length(Events$Drift)
      nearestTime= rep(NA,n)
      for (j in 1:n) {
        nearestTime[j]= which.min(abs(Events$ClickDateTime[j]-StaSpotData$UTC))
      }
      Events$lat= StaSpotData$lat[nearestTime]
      Events$long= StaSpotData$long[nearestTime]
      Events$SpotTime= StaSpotData$UTC[nearestTime]
      Events$deltaMin= as.numeric(abs(StaSpotData$UTC[nearestTime]-Events$ClickDateTime))/60

# make sure events are after deployment and before retrieval
      iDrift= DriftNum[ifile]
      Events= Events[(Events$ClickDateTime > DeployTime[iDrift]) & (Events$ClickDateTime < RetrTime[iDrift]),]  
    
    
# output Events with Lat/Lon info
      write.csv(Events,file=paste(sub('\\.sqlite3$','',files[ifile]),' Events_wLatLon.csv',sep=""))
# concatenate all events into one file
      if (ifile == 1) {
        AllEvents= Events
      }else{
        AllEvents= rbind(AllEvents,Events)
      }
    

# add event click variables    
      EventClicks$BinaryFile= trimws(as.character(EventClicks$BinaryFile))
      UniqueFiles= unique(EventClicks$BinaryFile)
      Binary_File_Date_Time= sub(pattern="Click_Detector_Click_Detector_Clicks_",replacement="",x=EventClicks$BinaryFile)   #get UTC from file start time
      Binary_File_Date_Time= sub(pattern=".pgdf",replacement="",x=Binary_File_Date_Time)                                                #get UTC from file start time
      EventClicks$FileDateTime= strptime(Binary_File_Date_Time, "%Y%m%d_%H%M%S", tz="gmt")
#correct Soundtrap file time/dates to get UTC
      if (DriftNum[ifile] == 1) {
        EventClicks$ClickDateTime=  EventClicks$UTC + 10.8*60*60             #add 10.8 hrs to correct filename to UTC
        EventClicks$FileDateTime=  EventClicks$FileDateTime + 10.8*60*60             #add 10.8 hrs to correct filename to UTC
      } else if (DriftNum[ifile]<=6 |DriftNum[ifile]==9 |DriftNum[ifile]==12 |DriftNum[ifile]==15 |DriftNum[ifile]==16 | DriftNum[ifile]>=26) {
        EventClicks$ClickDateTime=  EventClicks$UTC + 7.0*60*60              #add 7.0 hrs to correct filename to UTC
        EventClicks$FileDateTime=  EventClicks$FileDateTime + 7.0*60*60              #add 7.0 hrs to correct filename to UTC
      } else {
        EventClicks$ClickDateTime= EventClicks$UTC
      }
      EventClicks$eventType= NA
      EventClicks$nClicks= NA
      EventClicks$comment= NA
      EventClicks$lat= NA
      EventClicks$long= NA
      EventClicks$Drift= DriftNum[ifile]
#add new variables to Eventclicks files
      for (j in 1:length(Events$Id)) {
        EventTypeTF= (EventClicks$EventId == Events$Id[j])     #identify clicks from this event
        Events$nClicks[j]= sum(EventTypeTF)                      #get true count of number of clicks in event
        EventClicks$eventType[EventTypeTF]= trimws(as.character(Events$eventType[j]))
        EventClicks$comment[EventTypeTF]= trimws(as.character(Events$comment[j]))
        EventClicks$lat[EventTypeTF]= Events$lat[j]
        EventClicks$long[EventTypeTF]= Events$long[j]
        for (k in 1:length(UniqueFiles)) {
          FilesToCount= EventTypeTF & (EventClicks$BinaryFile==UniqueFiles[k])
          EventClicks$nClicks[FilesToCount]= sum(FilesToCount)    #count number of event clicks in each file for each event
        }
      }    
# eliminate events without any clicks      
     EventClicks= EventClicks[!is.na(EventClicks$nClicks),]
# make sure events are after deployment and before retrieval
     iDrift= EventClicks$Drift
     EventClicks= EventClicks[(EventClicks$FileDateTime > DeployTime[iDrift]) & (EventClicks$FileDateTime < RetrTime[iDrift]),]  


     summary(EventClicks)
    
#####################################################################################################
# read binary files for each DASBR/card and merge angle and amplitude information with EventClicks
    
# merge angle and amplitude information from binary files to event clicks
     BinaryClicks= read.csv(paste(dirname(files[ifile]),"\\","ClickBinaries.csv",sep=""))
     BinaryClicks$BinaryFile= trimws(as.character(BinaryClicks$BinaryFile))
     NumberClicks[ifile]= length(BinaryClicks$UID)
# get datetime of each click and correct for clock errors in seconds
    #BinaryClicks$datetime= strptime(as.character(BinaryClicks$datetime),format="%Y-%m-%d %H:%M:%OS",tz="gmt")
     BinaryClicks$datetime= convertPgDate(BinaryClicks$date)

     if (DriftNum[ifile] == 1) {
       BinaryClicks$datetime=  BinaryClicks$datetime + 10.8*60*60             #add 10.8 hrs to correct filename to UTC
     } else if (DriftNum[ifile]<=6 |DriftNum[ifile]==9 |DriftNum[ifile]==12 |DriftNum[ifile]==15 |DriftNum[ifile]==16 | DriftNum[ifile]>=26) {
       BinaryClicks$datetime=  BinaryClicks$datetime + 7.0*60*60              #add 7.0 hrs to correct filename to UTC
     }
# make time-of-day histograms for all clicks
    datezero= paste(substr(as.character(BinaryClicks$datetime),1,10),"00:00:00",sep=" ")          #beginning of the day
    datezero= strptime(datezero,format="%Y-%m-%d %H:%M:%S",tz="gmt")  
    histout= hist(as.numeric(difftime(BinaryClicks$datetime,datezero,units="hours"),breaks=0:24,plot=FALSE)) #time in hrs since beginning of the day
    HistSave= HistSave + histout$counts
# tally the total number of clicks of all types in each binary file
    UniqueBinaryFiles= unique(BinaryClicks$BinaryFile)
    TotClicksPerFile= rep(NA,length(UniqueBinaryFiles))
    for (j in 1:length(UniqueBinaryFiles)) {
      TotClicksPerFile[j]= length(BinaryClicks$BinaryFile[BinaryClicks$BinaryFile==UniqueBinaryFiles[j]])
    }
# create database with click counts per file 
    ClickCountPerFile= data.frame(BinaryFile= basename(UniqueBinaryFiles), TotClicks= TotClicksPerFile)
    ClickCountPerFile$BinaryFile= trimws(as.character(ClickCountPerFile$BinaryFile))
#create new smaller binary click dataframe with only needed data fields to speed up merge
    BinaryClicks2= data.frame(ClickNo=BinaryClicks$UID,Angle=BinaryClicks$angles*180/pi,type=BinaryClicks$type,
                              BinaryDateTime=BinaryClicks$datetime,BinaryFile=BinaryClicks$BinaryFile)
#create eventID name that is unique to a drift
    EventClicks$EventId= as.character(EventClicks$EventId)
    EventClicks$EventId= paste(EventClicks$EventId,CardNames[ifile],sep="_")
#create new smaller event click dataframe with only needed data
    EventClicks2= data.frame(ClickNo=EventClicks$UID,ClickDateTime=EventClicks$ClickDateTime,
                             eventType=EventClicks$eventType,EventId=EventClicks$EventId,                             
                             Drift=EventClicks$Drift,Amplitude=EventClicks$Amplitude,
                             comment=EventClicks$comment,FileDateTime=EventClicks$FileDateTime,
                             lat=EventClicks$lat,long=EventClicks$long,nClicks=EventClicks$nClicks)

# convert UIDs to characters to avoid too-large integers (in some cases)    
    BinaryClicks2$ClickNo= as.character(BinaryClicks2$ClickNo)
    EventClicks2$ClickNo= as.character(EventClicks2$ClickNo)
    
# merge clicks file with event file based on ClickNo to create new click file with event info
#   and eliminating clicks without events
    EventClicksMerged= merge(EventClicks2,BinaryClicks2,all.x=TRUE,all.y=FALSE) 
    MergedAllClicks= merge(BinaryClicks2,EventClicks2,all=TRUE,sort=TRUE) 
    
#    write.csv(EventClicksMerged,file="EventClicks Merged.csv")
#    write.csv(MergedAllClicks,file="AllClicks Merged.csv")

# calculate elapsed time of each click since beginning of file
# for continuously recording soundtraps, set file start time on the even minutes
    if (DriftNum[ifile]>=23) {
      timeA= as.character(EventClicksMerged$ClickDateTime)
      minutes= as.numeric(substr(timeA,15,16))
      seconds= as.numeric(substr(timeA,18,19))
      startminute= 2*floor(minutes/2)
      deltaTime= (minutes-startminute)*60 + seconds
# for other files time is taken from the binary file name
    } else {
      deltaTime= as.numeric(difftime(EventClicksMerged$ClickDateTime,EventClicksMerged$FileDateTime,units="secs"))  #time in sec since begging of file
    }
    EventClicksMerged$FileTime= deltaTime
# for each event, count number of clicks are present in either the first or the second minute of each file 
    EventClicksMerged$nCountLT60sec= NA
    EventClicksMerged$nCountGT60sec= NA
    allBinaries= unique(EventClicksMerged$BinaryFile)
    for (iBinary in 1:length(allBinaries)) {
      EventClicksMerged$nCountLT60sec[EventClicksMerged$BinaryFile==allBinaries[iBinary]]= 
        sum(deltaTime[EventClicksMerged$BinaryFile==allBinaries[iBinary]]<60)
      EventClicksMerged$nCountGT60sec[EventClicksMerged$BinaryFile==allBinaries[iBinary]]= 
        sum(deltaTime[EventClicksMerged$BinaryFile==allBinaries[iBinary]]>=60)
    }


#output files with event ID numbers and angles in EventClicks file
     if (ifile == 1) {
       AllClicksMerged= MergedAllClicks
       AllEventClicksMerged= EventClicksMerged
     } else {
       AllClicksMerged= rbind(AllClicksMerged,MergedAllClicks)
       AllEventClicksMerged= rbind(AllEventClicksMerged,EventClicksMerged)
     }
    cat(files[ifile],"\n")  
    
# calculate proportion number of 2-minute files with Ziphius click events 
    uniqueZCfiles= (unique(MergedAllClicks$BinaryFile[MergedAllClicks$eventType=="ZC"]))
    uniqueWAVfiles= (unique(MergedAllClicks$BinaryFile))
    cat(" Ziphius events in ",length(uniqueZCfiles)," of ",length(uniqueWAVfiles)," total files or ",
        100*length(uniqueZCfiles)/length(uniqueWAVfiles),"% \n")
  }
}
  
# Output merged files from all cards for all species 
Drift= AllEventClicksMerged$Drift[1]
write.csv(AllClicksMerged,file=paste("Drift-",Drift,"-AllClicks wAngle&LatLon.csv",sep=""))
write.csv(AllEventClicksMerged,file=paste("Drift-",Drift,"-AllEventClicks wAngle&LatLon.csv",sep=""))

# Output merged files from all drifts for beaked whales
BeakedWhales= c("ZC","MS","BWunid","BB","BW39V","BW43","BW","BW26-47","BW38")
BWClicksMerged= AllClicksMerged[AllClicksMerged$eventType %in% BeakedWhales,]
write.csv(BWEventClicksMerged,file=paste("Drift-",Drift,"-BWEventClicks wAngle&LatLon.csv",sep=""))


uniqueZCfiles= (unique(AllClicksMerged$BinaryFile[AllClicksMerged$eventType=="ZC"]))
uniqueWAVfiles= (unique(AllClicksMerged$BinaryFile))
cat(" Ziphius events in ",length(uniqueZCfiles)," of ",length(uniqueWAVfiles)," total files or ",
    100*length(uniqueZCfiles)/length(uniqueWAVfiles),"% \n")

  


