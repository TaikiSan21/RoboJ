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

# FULL RUN EXAMPLE FROM HERE
ec <- readRDS('Test/eventClicksPascal.rds')
snap <- 1
es <- formatEventSummary(ec, pascal=T, snapshot = snap)
es <- filter(es, species=='ZC')
es <- filter(es, !(station %in% c(6,15)))
detFunction <- newEstDetFunction(es, subsetCol = 'recorder',
                                 doJackknife = TRUE, jk_nSamples = 25, nSim=1e7, progress=TRUE)
detHist <- createDetectionHistory(filter(es, !(station %in% c(6,15))), snapshot = snap, maxTime=62)
ehData <- ehDataPrep(detHist, plot=F, plotDir = '../RoboJ/Test/EHPlot')
library(R2OpenBUGS)

ehSpecs <- list(
    ni = 10000,  # chain length, including burn-in
    nt = 2,      # thinning rate
    nb = 2500,   # burn-in (e.g., 4000 * thin-rate 5 = 20K samples per chain)
    nc = 2       # number of chains
)

##### Set parameters to monitor #####
# original phi,p changed to *1
ehParams = c("phi1", "p1",
              "mean.availtime.min")

##### Source model file #####
ehModel <- makeEhModel()
##### Source initial-values file #####
# No initial-values file created.  Let BUGS set own random initial values by specifying inits=NULL in bugs() function.

##### Run BUGS analysis #####
ehOut <- bugs(data=ehData, inits=NULL,
            parameters.to.save=ehParams,
            model.file=ehModel, n.chains=ehSpecs$nc, n.iter=ehSpecs$ni,
            n.burnin=ehSpecs$nb, n.thin=ehSpecs$nt,
            saveExec=FALSE, debug=TRUE, OpenBUGS.pgm=NULL, working.directory=getwd(), clearWD=TRUE)
saveRDS(ehOut, file=here('availability_estimation', 'MCMC_availEst_2023-01-19_snap2.rds'))
ehOut <- readRDS(here('availability_estimation', '2023-01-17_MCMCsamples_EH.rds'))

# # rjags version
# library(rjags)
# ehJags <- jags.model(
#     file =ehModel,
#     data = ehData,
#     inits = list(p1=.3, phi1=.9),
#     n.chains = ehSpecs$nc,
#     n.adapt = 100 #100 from steph
# )
#
#
# # Burnin ------------------------------------------------------------------
#
# update(ehJags, n.iter = ehSpecs$nb)
#
#
# # Collect posterior samples -----------------------------------------------
#
#
# j <- jags.samples(
#     model = ehJags,
#     variable.names = ehParams,
#     n.iter = ehSpecs$ni,
#     thin = ehSpecs$nt
# )

# 15.39??/.5507 not sure where 14.0 was coming from before, original run was 15.23/.6623
# 2023-1-17 run 2 snap: .9444, .2582, 15.19, 889. Deviance fixed.
# 2023-1-18 run 1 snap: .9366, .379, 13.95, 1738.

# next input params:
# ESR and ESR_SD from Jay sim like model
# MAT + snap length, SD avail time from

erSpecs <- list(
    ni = 200000,  # chain length, including burn-in
    nt = 5,     # thinning rate
    nb = 50000,   # burn-in
    nc = 2       # number of chains
)
##### source the model file #####
# source(here('density_estimation', 'encounterRate-model_randomDensities_updated.R'))
erModel <- makeErModel()

##### Source initial-values file #####
# No initial-values file created.  Let BUGS set own random initial values by specifying inits=NULL in bugs() function.

##### Set parameters to monitor #####
erParams <- c("v.mean", "er.mean",                                 # effective search area and radius
               "logDgrp.hyper", "logDgrp.sig",                      # random effect parameters for log-density
               "logDgrp", "logDgrp.err", "Dper1000", "mu",          # parameters for individual DASBRs
               "mean.of.Errs", "SD.of.Errs", "mean.of.logDgrp",     # descriptive summary stats for individuals DASBR estimates (model checks)
               "Dgrp.mean", "D.mean.per1000", "N",                  # mean density and abundance for the study area (use)
               "Dgrp.mean.II","D.mean.per1000.II",  "N.II",         # mean density and abundance for the study area, alternate (do not use)
               "s", "g0")                                           # ancillary parameters
snaps <- data.frame(station = as.character(1:22),
                    nSnaps = c(1629, 1630, 1762, 1841, 1918, 3078, 13669, 6807, 2812, 13996, 7073, 2801, 12133, 0, 3135, 2800, 14250, 6934, 2707, 3329, 3087, 3376)
)
snaps$nSnaps <- snaps$nSnaps * 2 / snap
# Values from previous steps
# esrMean, esrSd, availMean (model+snaplen), availSd,
# Constants from other studies:
# deep dive period: ddpMean=191.4, ddpSd=28.2
# group size: gsizeMean=1.9, gsizeSd=.13
erData <- erDataPrep(es, snaps=snaps,
                     detFunction = detFunction$ST4300,
                     # esrMean = 3.13,
                     # esrSd = .313,
                     ehOut=ehOut,
                     ddpMean=191.4,
                     ddpSd=28.2,
                     gsizeMean=1.9,
                     gsizeSd = .13,
                     plot=FALSE, plotDir='Test/ERPlot')

# # erData$esrMean <- 3.13 # comp to 3.057, 3.000(snap1)
# erData$esrSd <- .313 # comp to .2686, .326(snap1)

# erData$availMean <- 15.23+2
# erData$availSd <- .6623

# 5,53 on all (3.0848, .3029)
# 5.023 on st4300(3.194, .2162)
# 5.336 on ols (3.13, .313)
##### Run BUGS analysis #####
# 5.633 on latest update, was 5.4 when using old esr/sd
erOut <- bugs(data=erData, inits=NULL, parameters.to.save=erParams, model.file=erModel,
            n.chains=erSpecs$nc, n.iter=erSpecs$ni, n.burnin=erSpecs$nb,
            n.thin=erSpecs$nt, saveExec=FALSE, debug=TRUE, OpenBUGS.pgm=NULL,
            working.directory=getwd(), clearWD=TRUE, bugs.seed=2)
saveRDS(erOut, file='density_estimation/MCMC_abundanceEst_2023-01-19_snap2.rds')
# 2023-1-19 2snap: d mean 5.53
# 2023-1-18 1snap: d mean 4.83
myOut <- outN
load('Jeffs Files/DensityEstimation/out/MCMCsamples_abundanceEstimates_4kTrunc.Rdata', verbose = TRUE)
out$mean$D.mean.per1000
