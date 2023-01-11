# RoboJay script
# Updated 4-20-2022 7PM - Adding CCES REGEX, Splitting Station 4, added more notes
library(PAMpal)
source('Taikis FIles/RoboJayFunctions.R')
# Place to save outputs
outPath <- 'Test'
# Can test on just one station, if that works then we can run all stations at once
db <- 'Test/Station-7_Card-A_MASTER-BW - Copy.sqlite3'
bin <- 'Test/Binaries_Station-7_Card-A_WithUID/'

pps <- PAMpalSettings(db, bin, sr_hz='auto', winLen_sec=.0025, filterfrom_khz=10, filterto_khz=NULL)
data <- processPgDetections(pps, mode='db')
data <- setSpecies(data, method='pamguard')
# Re-label species codes if necessary for consistency
speciesRename <- data.frame(old = c('BW34-50', 'BW46', 'BW50-75', '?BW', 'Zc', 'Pm', 'SW', '?Pm', 'Oo'),
                           new = c('BW37V', 'MS', 'MS', 'BWunid', 'ZC', 'PM', 'PM', '?PM', 'OO'))
data <- setSpecies(data, method = 'reassign', value=speciesRename)
# Make sure this doesnt have any codes that didnt get converted above
unique(species(data))

data <- addGps(data)
# Folder containing wavs. Easy if just one DB:
recordingFolder <- file.path(outPath, 'Recordings')
# If multiple DBs, first check order of DBs:
files(data)$db
# Then folder needs to be a vector in that order specifying which folder for each DB:
recordingFolder <- c(
    'Test/RecordingsDB1',
    'Test/RecordingsDB2'
)
# Alternatively set folder=NULL and it will do a pop-up asking you to choose folder for each DB if thats easier
data <- addRecordings(data, folder=recordingFolder)
# saveRDS(data, file=file.path(outPath, 'RoboJayStudy.rds'))

# some formatting is particular to pascal only
pascal <- FALSE
# stationPattern is a REGEXP to turn event names into a single "Station" name that is
# easier to write down.
# Default gets the name of the database from PAMpal's eventIds that are created when processing mode='db'
defaultPat <- '(.*)\\.OE.*|\\.DGL.*'
# This should convert CCES to single numbers
CCESPat <- '.*Drift-([0-9]{1,2})_.*[\\.OE|\\.DGL][0-9]+$'
eventClicks <- formatClickDf(data, pascal=pascal, stationPattern = defaultPat)

# Preps for the soundspeed profile download - marks detections into different groups
# first by the "splitBy" column, then new groups if separated by longer than "time"
# If download times seem to be long you can increase the "time" value to maybe 3 or 7 days
eventClicks <- splitSSPGroups(eventClicks, splitBy='station', time=3600*24*1)

# This creates a SSP for each row, has to download from HYCOM server so may be slow
# Will save the downloaded data to "file" so that it can be read back in and used in the future
sspList <- createSSPList(eventClicks, file=file.path(outPath, paste0(as.character(Sys.Date()), '_sspList.rds'))) #may need to add date

# Can read this from disk if running again in future on same data (change date for future)
# sspList <- readRDS(file.path(outPath, '2022-04-20_sspList.rds'))

# Uses the SSP to correct angles using raytrace, then applies raytrace angle correction
# to existing angles. May be quite slow since raytracing can take some time
eventClicks <- doAngleCorrection(eventClicks, sspList)

# "recorderInfo" has columns "station" "dutyCycle" and "recorder", "station" must match exact formatting
# deployDetails.csv is from our GDrive
recorderInfo <- formatDeployDetails('Test/deployDetails.csv')
recorderInfo <- filter(recorderInfo, Project == 'CCES')
# Split Station 4 - first 20days 2/2, then 2/18
st4new <- recorderInfo[recorderInfo$station == '4', ]
st4new$station <- '4_2'
st4new$dutyCycle <- '2/18'
st4new$deployTime <- st4new$deployTime + 20 * 24 * 3600
recorderInfo <- rbind(recorderInfo, st4new)
recorderInfo$station[recorderInfo$station == '4'] <- '4_1'
recorderInfo$dutyCycle[recorderInfo$station == '4_1'] <- '2/2'
recorderInfo$retrTime[recorderInfo$station == '4_1'] <- recorderInfo$deployTime[recorderInfo$station == '4_2']
# Make changes for eventClicks data. Set all to 4_1, then any after 4_2 deploy time to 4_2
eventClicks$station[eventClicks$station == '4'] <- '4_1'
eventClicks$station[eventClicks$station == '4_1' & eventClicks$UTC >= st4new$deployTime] <- '4_2'

eventSummary <- formatEventSummary(eventClicks, recorderInfo = recorderInfo, pascal=pascal, snapshot=2)

recordingDf <- files(data)$recordings
# Do any wav file filtering here - removing any files example from PASCAL:
# recordingDf$keepWavs <- (regexpr(pattern=glob2rx("*201883689*"),text=recordingDf$file) < 0) &   #remove Soundtrap ST300HF-1 files
#     (regexpr(pattern=glob2rx("*1678299174*"),text=recordingDf$file) < 0) &   #remove Soundtrap ST300HF-Jason files
#     (regexpr(pattern=glob2rx("*Stitched*"),text=recordingDf$file) < 0) &   #remove old stitched-together soundtrap files
#     (regexpr(pattern=glob2rx("*DASBR3*"),text=recordingDf$file) < 0)   #remove Station 14 with unuseable recordings

# This output shouldnt take any time to run
snapshots <- tallySnapshots(recordingDf, eventSummary, recorderInfo)

# Do any desired filtering after this point, examples below
eventSummary <- filter(eventSummary, species == 'ZC')
# eventSummary <- eventSummary[eventSummary$Station <= 22, ]

# Also takes no time to run
detHist <- createDetectionHistory(eventSummary)

# This will take quite a long time. To test and make sure things just run you can set doJackknife = FALSE and nSim = 1e6,
# but for full run (> 1hr ) we want doJackknife=TRUE and nSim=1e7
detFunction <- newEstDetFunction(eventSummary, subsetCol = 'recorder',
                                 doJackknife = TRUE, jk_nSamples = 5, nSim=1e7, progress=TRUE)

saveRDS(eventSummary, file.path(outPath, 'eventSummary.rds'))
saveRDS(detFunction, file.path(outPath, 'detFunction.rds'))
saveRDS(snapshots, file.path(outPath, 'snapshots.rds'))
saveRDS(detHist, file.path(outPath, 'detHist.rds'))