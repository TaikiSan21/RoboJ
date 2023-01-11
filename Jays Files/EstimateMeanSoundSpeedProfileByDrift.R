# Calculate sound speed profile at the mean location of each PASCAL drift
#  using mean data from September from the Global Ocean Sound Speed Profile Library 
#  (GOSSPL) Rdata structures (Barlow 2019)
#  https://swfsc-publications.fisheries.noaa.gov/publications/TM/SWFSC/NOAA-TM-NMFS-SWFSC-612.pdf

library("gam")    #needed to produce smoothed sound speed profile at 1-m scale

iMonth= 9  #this is September
SSlibraryLocation= "S:\\SoundSpeedWorldwideData\\WorldSoundSpeedData by Month"
aMonth= month.abb[iMonth]    #abbreviated month names for naming output files
load(paste(SSlibraryLocation,"\\WorldSoundSpeedData_",aMonth,".RData",sep=""))
OutputDirectory= "C:\\Jay\\ACOUSTIC\\Buoy Recorder\\PASCAL\\R_Analysis_Output\\"

MeanLocationsByDrift= read.csv(file= paste(OutputDirectory,"MeanLocationsByDrift.csv",sep=""))
Drifts= unique(MeanLocationsByDrift$Drift)
nDrifts= length(Drifts)
SoundSpeed= array(dim=c(length(DepthValues),nDrifts+1))
maxDepth= rep(NA,nDrifts)
SoundSpeed[,1]= DepthValues
iDrift= 1

# get mean monthly sound speeds at up to 78 depths for a given location from global database
for (i in 1:length(Drifts)) {
  iDrift= Drifts[i]
  Lat= MeanLocationsByDrift$meanLat[MeanLocationsByDrift$Drift==iDrift]
  Long= MeanLocationsByDrift$meanLon[MeanLocationsByDrift$Drift==iDrift]
  # plot sound speed profiles (south and west are negative lat/long)
  iLat= which.min(abs(LatValues-Lat))
  iLong= which.min(abs(LongValues-Long))
  SoundSpeed[,i+1]= SoundSpeedMatrix[iLong,iLat,]
  
  imaxDepth= which(is.na(SoundSpeedMatrix[iLong,iLat,1:78]))
  imaxDepth= imaxDepth[1]-1
  Depths= DepthValues[1:imaxDepth]
  maxDepth[i]= max(Depths)
  #plot(SoundSpeed[1:imaxDepth,i+1],-Depths,xlab="Sound Speed (m/s)",ylab="Depth",main=(paste("Sound Speed Profile for",aMonth,"Drift=",iDrift,sep=" ")))

}

SoundSpeedDF= data.frame(SoundSpeed)
names(SoundSpeedDF)= c('Depth',paste("Drift_",as.character(Drifts),sep=""))     
write.csv(SoundSpeedDF,file=paste(OutputDirectory,"SoundSpeedProfilesByDrift.csv",sep=""))

# use gam smooth to interpolate sound speed at 1-m intervals 
for (i in 1:length(Drifts)) {
  iDrift= Drifts[i]
  imaxDepth= which(maxDepth[i]==DepthValues)
  Depth= SoundSpeedDF$Depth[1:imaxDepth]
  SSpeed= SoundSpeedDF[1:imaxDepth,i+1]
  gamout= gam(SSpeed~s(Depth,30))
  newdata= data.frame(Depth=0:maxDepth[i])
  predicted= as.numeric(predict.gam(gamout,newdata=newdata))
  newdata$SoundSpeed= predicted
  write.csv(newdata,file=paste(OutputDirectory,"SoundSpeedProfileDrift_",iDrift,sep=""))
  plot(SSpeed,-Depth,xlab="Sound Speed (m/s)",ylab="Depth",main=(paste("Sound Speed Profile for",aMonth,"Drift=",iDrift,sep=" ")))
  lines(predicted,-(0:maxDepth[i]), col="red", lwd= 2)
}
