# 2021 April 8, Jeff Moore
# This file:
  # sources EH-dataPrep.R to generate the BUGS data object
  # sources BUGS model file (the CJS model)
  # specifies OpenBUGS run arguments (e.g., MCMC parameters, parameters to monitor)
  # runs OpenBUGs to obtain estimates from the availability model
  # saves bugs() output object as an Rdata file
  # OpenBugs .odc and .txt output files are manually saved within this process

##### Libraries #####
library(R2OpenBUGS)

##### Set working directory #####
setwd("C:/jeff/NOAA/PASCAL/analysis/finalAnalysisDataAndCode-March2021/AvailabilityEstimation/")

##### Call EH-dataPrep.R to read in encountery history data, create BUGS data object, and generate data summary plots #####
source("EH-dataPrep.R")

##### Set BUGS MCMC parameters #####
mcmc.specs <- list(
  ni = 5000,  # chain length, including burn-in
  nt = 2,      # thinning rate
  nb = 1000,   # burn-in (e.g., 4000 * thin-rate 5 = 20K samples per chain)
  nc = 2       # number of chains
)

##### Set parameters to monitor #####
# original phi,p changed to *1
OB.params = c("phi1", "p1",                                     # CJS parameters
              "mean.availtime.min" #, "mean.availtime.ints",     # Mean availability time (this is what we want)
              # "med.availtime.min", "med.availtime.ints"        # Median availability time
              )

##### Source model file #####
source("EH-model p(.) phi(.).R")
source(here('availability_estimation', 'EH-model_unequal.R'))

##### Source initial-values file #####
# No initial-values file created.  Let BUGS set own random initial values by specifying inits=NULL in bugs() function.

##### Run BUGS analysis #####
out <- bugs(data=ehData, inits=NULL,
            parameters.to.save=OB.params,
            model.file=modelFilename, n.chains=mcmc.specs$nc, n.iter=mcmc.specs$ni, n.burnin=mcmc.specs$nb, n.thin=mcmc.specs$nt,
            saveExec=FALSE, debug=TRUE, OpenBUGS.pgm=NULL, working.directory=getwd(), clearWD=TRUE)

##### Save the BUGS output.  Use load() to read these back in at a later date #####
save(out, file = "out/MCMCsamples_EH_p(.)phi(.).Rdata")  # phi(.) model
#load("out/MCMCsamples_EH_p(.)phi(.).Rdata")  # make sure higher-level working directory is set (above)



# phi1 .94 p1 .191 MAT 10.3  NOTE: .94^2 = .8836, 1-(1-.191)^2 = .3458 MAT CHANGED TO 15.05 WHEN I UPPED NOCC
# phi .883 p .346 MAT 15.09
# 2022-12-22 2snaprun phi1 .9371 p1 .2863 MAT 14.04. 10k, 2, 2.5k, 2. Chains settled well before end
# 2022-12-22 1snaprun phi1 .9349 p1 .3964 MAT 13.68 5k, 2, 1k, 2. STill pretty settled.