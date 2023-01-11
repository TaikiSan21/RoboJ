es <- readRDS('Test/eventSummaryCCES.rds')

ec <- readRDS('Test/eventClicksPascal.rds')
es <- formatEventSummary(ec, pascal=T, snapshot = 2)

# es <- filter(es, dutyCycle != '2/3')

detHist <- createDetectionHistory(filter(es, species == 'ZC'), snapshot = 2)
dh1 <- createDetectionHistory(es, snapshot=1)
View(dh1$dhMat)
ehData <- ehDataPrep(detHist, plot=F, plotDir = '../RoboJ/Test/EHPlot')

# change dataprep to output everything for GOF step


es <- formatEventSummary(ec, snapshot=1, pascal=TRUE)
detHist <- createDetectionHistory(filter(es, species == 'ZC'), snapshot = 1)
ehData <- ehDataPrep(detHist, plot=F, plotDir = '../RoboJ/Test/EHPlot')

es1 <- formatEventSummary(ec, snapshot=1, pascal=TRUE)
detHist1 <- createDetectionHistory(es1, snapshot = 1)
# have to fake this bc stdy didnt have end times for recordings
data <- readRDS('Test/RoboJayStudy.rds')
rec <- files(data)$recordings
rec <- rec[!is.na(rec$start),]
rec$end <- rec$start + 119 + runif(nrow(rec))
rec$length <- as.numeric(difftime(rec$end, rec$start, units='secs'))
snaps <- tallySnapshots(rec, es, recorderInfo = recorderInfo, snapshot=2)
snaps <- data.frame(station = as.character(1:22),
                    nSnaps = c(1629, 1630, 1762, 1841, 1918, 3078, 13669, 6807, 2812, 13996, 7073, 2801, 12133, 0, 3135, 2800, 14250, 6934, 2707, 3329, 3087, 3376)
)

erData <- erDataPrep(filter(es, species == 'ZC'), snaps)
detFun <- readRDS('Test/detFunction.rds')

# ok notes. I think snapshot change is okay. We just need to require that the
# "on" duration is a multiple of snapshot, but the "off" we can just divide
# as any value using updated code.

# Q about comment of snapshot minutse - I'm not finding where this proportion
# gets calculated, only seeing that total number of snapshots per station
# is input to density step (k). I think this is maybe inherent in the comparison
# to n below?

# estimates are made as GroupDensity * k (n snaps) * area. Then these are compared
# to total number of events (defined earlier with snapshot length) per station (n)
# when fitting.

# okay so result of above is i think we dont need to count minutes if we force
# the "on" durations to be mutliple of snapshot length?