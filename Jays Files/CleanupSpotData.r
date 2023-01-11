# Cleans up Spot location data for individual drifts and merges it
#    convert to actual UTC, eliminate speed outliers and filter out locations before and after deployment
#    also estimates deployment duration and distances traveled.

# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
setwd("C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/Files for Submission/Supplemental R-code & data")

library(geosphere)

##############################################################################################
# get deployment and retrieve times in Excel format and convert to R format
DeployRetrTimes= read.csv("DeplRetrDateTimes.csv")
DeployTime= as.POSIXct((DeployRetrTimes$UTC_DeployTime-2)*24*60*60,origin='1900-01-01',tz='gmt')
RetrTime= as.POSIXct((DeployRetrTimes$UTC_RetrTime-2)*24*60*60,origin='1900-01-01',tz='gmt')


######################################################################################################
## read Spot geolocation info and save in a single file with UTC times added
# NOTE downloaded Spot data were in local time (UTC - 7)
sumKms= rep(NA,30); sumHrs= rep(NA,30)
for (iDrift in 1:30) {
  SpotData= read.csv(file=paste(getwd(),"/SpotLocationDataByDrift/spotTrack-station",as.character(iDrift),".csv",sep=""))
  
  # convert date/time to UTC and store as number  
  UTC=  strptime(SpotData$dateTime, "%m/%d/%Y %H:%M", tz="gmt") + (7*60*60)
  # eliminate positions that are < deploy and > retrieve time
  SpotData= SpotData[(UTC > DeployTime[iDrift]) & (UTC < RetrTime[iDrift]),]  
  
  #eliminate distance outliers
  n= length(SpotData$lat)
  dist= distGeo(cbind(SpotData$long[1:(n-1)],SpotData$lat[1:(n-1)]),cbind(SpotData$long[2:n],SpotData$lat[2:n])) / 1000   #distance in km
  dist[n]= 0
  #cat(iDrift,max(dist),which.max(dist),"\n")
  outliers= which(dist>2) + 1
  if (length(outliers) > 0) {
    SpotData$lat[outliers]= NA
    SpotData= SpotData[!is.na(SpotData$lat),]
  }
  # recalculate distances without outliers
  n= length(SpotData$lat)
  dist= distGeo(cbind(SpotData$long[1:(n-1)],SpotData$lat[1:(n-1)]),cbind(SpotData$long[2:n],SpotData$lat[2:n])) / 1000   #distance in km
  #cat(iDrift,max(dist),which.max(dist),"\n")
  
  # re-convert date/time to UTC and store as number  
  UTC=  strptime(SpotData$dateTime, "%m/%d/%Y %H:%M", tz="gmt") + (7*60*60)
  
  # compile all good positions into new dataframe
  if (iDrift == 1) {
    SpotData_wUTC= data.frame(Drift=iDrift,dateTime=SpotData$dateTime,spotID=SpotData$spotID,lat=SpotData$lat,long=SpotData$long,UTC)
  } else {
    SpotData_wUTC= rbind(SpotData_wUTC,data.frame(Drift=iDrift,dateTime=SpotData$dateTime,spotID=SpotData$spotID,lat=SpotData$lat,long=SpotData$long,UTC))
  }
  
  sumKms[iDrift]= sum(dist)
  sumHrs[iDrift]= (as.numeric(UTC[1])-as.numeric(UTC[n]))/(60*60)
  cat(iDrift,"   Tot Dist= ",sumKms[iDrift],"km   Tot Hrs=",sumHrs[iDrift],"hrs   AvgSpeed=",sumKms[iDrift]/sumHrs[iDrift]," km/hr \n")
}
TotKms= sum(sumKms[c(1:13,15:22)])
TotHrs= sum(sumHrs[c(1:13,15:22)])
write.csv(SpotData_wUTC,file="AllSpotData wUTC.csv")
cat("Overall","   Tot Dist= ",TotKms,"km   Tot Hrs=",TotHrs,"hrs   AvgSpeed=",TotKms/TotHrs," km/hr \n")



