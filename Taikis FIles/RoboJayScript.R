# RoboJay script
library(PAMpal)
# Can test on just one station, if that works then we can run all stations at once
db <- '../Data/RoboJ/Test/Station-7_Card-A_MASTER-BW - Copy.sqlite3'
bin <- '../Data/RoboJ/Test/Binaries_Station-7_Card-A_WithUID/'
#### These are commented out, we'll need to sort out the GPS adding for the real data
# gps <- '../Data/RoboJ/Jays Files/AllSpotData wUTC.csv'
# gpsDf <- read.csv(gps, stringsAsFactors = FALSE)
# gpsDf <- rename(gpsDf, Longitude=long, Latitude=lat)
# gpsDf <- filter(gpsDf, Drift == 7)
# PAMmisc::addPgGps(db, gpsDf, format='%Y-%m-%d %H:%M:%S')
####
pps <- PAMpalSettings(db, bin, sr_hz='auto', winLen_sec=.0025, filterfrom_khz=10, filterto_khz=NULL)
data <- processPgDetections(pps, mode='db')
data <- setSpecies(data, method='pamguard')
pascalRename <- data.frame(old = c('BW34-50', 'BW46', 'BW50-75', '?BW', 'Zc', 'Pm', 'SW', '?Pm', 'Oo'),
                           new = c('BW37V', 'MS', 'MS', 'BWunid', 'ZC', 'PM', 'PM', '?PM', 'OO'))
data <- setSpecies(data, method = 'reassign', value=pascalRename)
data <- filter(data, species != 'Ship')
data <- addGps(data)
# saveRDS(data, file='RoboJayStudy.rds')

eventClicks <- formatClickDf(data, pascal=TRUE)
# This should match AllEventClicks.csv but ours have UID instead of ClickNo
# write.csv(jayRenamer(eventClicks), 'AllEventClicks_RoboJ.csv', row.names = FALSE)
# eventClicks <- read.csv('../Data/RoboJ/Test/AllEventClicks_RoboJ_wDrift13K.csv', stringsAsFactors = FALSE)
# eventClicks$ClickDateTime <- as.POSIXct(eventClicks$ClickDateTime, format='%Y-%m-%d %H:%M:%S', tz='UTC')
# eventClicks$UTC <- eventClicks$ClickDateTime
# eventClicks$Latitude <- eventClicks$lat
# eventClicks$Longitude <- eventClicks$long
###################################################################################
################ Can stop here and compare AllEventClicks before continuing #######
# Download process in next step may be time consuming, so confirm that previous
# steps are working well first
###################################################################################

# Preps for the soundspeed profile download - marks detections into different groups
# first by the "splitBy" column, then new groups if separated by longer than "time"
# If download times seem to be long you can increase the "time" value to maybe 3 or 7 days

eventClicks <- splitSSPGroups(eventClicks, splitBy='eventId', time=3600*24*1)

# This creates a SSP for each row, has to download from HYCOM server so may be slow
sspList <- createSSPList(eventClicks, file='../Data/RoboJ/Test/sspListOneTest.rds')
# saveRDS(sspList, file='../Data/RoboJ/Test/sspList.rds')
sspList <- readRDS('../Data/RoboJ/Test/sspList.rds')
# Uses the SSP to correct angles using raytrace, then applies raytrace angle correction
# to existing angles. May be quite slow since raytracing can take some time
eventClicks <- doAngleCorrection(eventClicks, sspList)
# saveRDS(eventClicksRC, file='../Data/RoboJ/Test/eventClicksRC.rds')
# eventClicksRC <- readRDS('../Data/RoboJ/Test/eventClicksRC.rds')
# for(i in 1:nrow(pascalRename)) {
#     eventClicksRC$eventType[eventClicksRC$eventType == pascalRename$old[i]] <- pascalRename$new[i]
# }

# write.csv(jayRenamer(eventClicks), 'AllEventClickswRC_RoboJ.csv', row.names = FALSE)
# source('../Data/RoboJ/Taikis FIles/RoboJEstDetFun.R')
eventSummary <- calcModalAngle(eventClicks, pascal=TRUE)
# saveRDS(eventSummary, '../Data/RoboJ/Test/eventSummary.rds')
eventSummary <- readRDS('../Data/RoboJ/Test/eventSummary.rds')
eventSummary <- eventSummary[eventSummary$Station <= 22, ]

detHist <- createDetectionHistory(eventSummary)
testMyDetSmol <- newEstDetFunction(eventSummary, doJackknife = TRUE, JK_nsamples = 20, verbose=FALSE, nSim=1e7, progress=F)
# from robojesting functions - maybe good as a debugger / sanity checker if odd param values
ll <- makeLLArr(seq(from=1e3, to=3e3, length.out=50), seq(from=-2, to=10, length.out=25), FUN=testMyDetSmol$All$likeFun(testMyDetSmol$All$Angle))
ll <- readRDS('../Data/RoboJ/Test/parrFast.rds')
llfix <- ll
llfix[llfix < quantile(llfix, .95)] <- quantile(llfix, .95)
image(llfix)
saveRDS(testMyDet, '../Data/RoboJ/Test/testMyDet.rds')
detEst <- estimateDetFunction(eventSummary, recorders = 'All Recorders', gamk = 25)
detEstJay <- estimateDetFunction(jayData, recorders = 'All Recorders', gamk = 25) #1919 .609, 1661 1.34
detEstJay1e6 <- estimateDetFunction(jayData, recorders = 'All Recorders', gamk = 25) #1919 .609, 1661 1.34, 1427 2.22
detEst50 <- estimateDetFunction(eventSummary, recorders = 'All Recorders', gamk = 50)

detEst1e5_25 <- estimateDetFunction(eventSummary, recorders = 'All Recorders', gamk = 25) # 1875 .69 -2486, 1452 2.12 -2485, 1438, 2.17
detEst1e5_25_HN <- estimateDetFunction(eventSummary, Model = 'HN', recorders = 'All Recorders', gamk = 25)
detEst1e5_50_HN <- estimateDetFunction(eventSummary, Model = 'HN', recorders = 'All Recorders', gamk = 50)

esm3 <- rename(esm3, eventType = species)
ccesm3 <- estimateDetFunction(esm3, recorders='All Recorders', gamk=10)

plotDetFun(c(1490, 1.9), add=T, col='blue')

recordingFolder <- '../Data/RoboJ/Test/Recordings'
data <- addRecordings(data, folder=recordingFolder)