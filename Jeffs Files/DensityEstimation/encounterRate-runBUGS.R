# 2021 April 9, Jeff Moore
# This file:
# sources encounterRate-dataPrep.R to generate the BUGS data object
# sources BUGS model file (the density/abundance estimation model)
# specifies OpenBUGS run arguments (e.g., MCMC parameters, parameters to monitor)
# runs OpenBUGs to obtain estimates from the abundance model
# saves bugs() output object as an Rdata file
# OpenBugs .odc and .txt output files are manually saved within this process

##### Libraries #####
library(R2OpenBUGS)

##### Set working directory #####
setwd("C:/jeff/NOAA/PASCAL/analysis/finalAnalysisDataAndCode-March2021/DensityEstimation/")

##### Call encounterRate-dataPrep.R to read in encounter frequency data, create BUGS data object, and generate data summary plots #####
source("encounterRate-dataPrep.R")

##### Set BUGS MCMC parameters #####
# mcmc.specs <- list(
#   ni = 200000,  # chain length, including burn-in
#   nt = 5,     # thinning rate
#   nb = 50000,   # burn-in
#   nc = 2       # number of chains
# )

mcmc.specs <- list(
    ni = 200000,  # chain length, including burn-in
    nt = 5,     # thinning rate
    nb = 50000,   # burn-in
    nc = 2       # number of chains
)
##### source the model file #####
source("encounterRate-model randomDensities.R")

##### Source initial-values file #####
# No initial-values file created.  Let BUGS set own random initial values by specifying inits=NULL in bugs() function.

##### Set parameters to monitor #####
OB.params <- c("v.mean", "er.mean",                                 # effective search area and radius
               "logDgrp.hyper", "logDgrp.sig",                      # random effect parameters for log-density
               "logDgrp", "logDgrp.err", "Dper1000", "mu",          # parameters for individual DASBRs
               "mean.of.Errs", "SD.of.Errs", "mean.of.logDgrp",     # descriptive summary stats for individuals DASBR estimates (model checks)
               "Dgrp.mean", "D.mean.per1000", "N",                  # mean density and abundance for the study area (use)
               "Dgrp.mean.II","D.mean.per1000.II",  "N.II",         # mean density and abundance for the study area, alternate (do not use)
               "s", "g0")                                           # ancillary parameters

##### Run BUGS analysis #####
out <- bugs(data=erData, inits=NULL, parameters.to.save=OB.params, model.file=modelFilename,
            n.chains=mcmc.specs$nc, n.iter=mcmc.specs$ni, n.burnin=mcmc.specs$nb,
            n.thin=mcmc.specs$nt, saveExec=FALSE, debug=TRUE, OpenBUGS.pgm=NULL,
            working.directory=getwd(), clearWD=TRUE, bugs.seed=2)

# 2022-12-22 run snap2 d.mean 6.01, dgrp 2.844e-4

##### Save the BUGS output.  Use load() to read these back in at a later date #####
save(out, file = "out/MCMCsamples_abundanceEstimates_4kTrunc.Rdata")
load("out/MCMCsamples_abundanceEstimates_4kTrunc.Rdata")

