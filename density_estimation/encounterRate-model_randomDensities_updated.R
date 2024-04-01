# 2021 April 9, Jeff Moore
# called by encounterRate-runBUGS.R
# This is the final step of the distance-sampling model, fit to encounter rate, to estimate population size and density
# The detection function is not estimated here.  That is done separately in a maximum simulated liklihood framework.
# The effective search radius (ESR) estimates from the max sim lik method are passed to this code
# This code then just estimates population density and abundance given the ESR, encounter rates (data), and ancillary info such as group size
# Individual DASBR drifts are treated as random effects for density estimation

# Bibliography
# Barlow J. 2016...
# Barlow J et al. 2021...


##### give name to text file that will hold WinBUGS model code (should match name of R file to facilitate analyses) #####
# modelFilename="C:/jeff/NOAA/PASCAL/analysis/finalAnalysisDataAndCode-March2021/DensityEstimation/ziphiusDensityModel.txt"
modelFilename <- here('density_estimation', 'ziphiusDensityModel.txt')

################################################################################################
##### Define the model in BUGS code and write to text file, which will be called by bugs() #####
cat("

    model{

      ##### EFFECTIVE SEARCH RADIUS AND AREA #####
      v.mean <- er.mean * er.mean * 3.14159  # effective search area = pi * r^2
      er.mean  ~ dnorm(esrMean,tau.er)     # r, the effective search radius (3.13 km), from Jay's max sim lik model (for 4km truncation, Jan 2021)
      tau.er <- 1/(sd.er * sd.er)
      sd.er <- esrSd                  # Jay provided CV = 0.10 (so SD = 0.313)

      # This loop accommodates different ESRs for each DASBR, but we ultimately used a single ESR, so er[station] = er.mean
      for(z in 1:n.sites){
        er[stationNames[z]] <- er.mean
        v[stationNames[z]] <- 3.14159 * er[stationNames[z]] * er[stationNames[z]]
      } # z


      ##### MODEL FOR DASBR-SPECIFIC DENSITIES #####
      for(z in 1:n.sites){                              # loop through DASBR stations
        n[z] ~ dpois(mu[z])                             # liklihood.  Number of *groups* detected per station (n[z] are the data)
        mu[z] <- Dgrp[z] * v[stationNames[z]] * k[z]    # mu = D*v*k = expected number of detections, where k is number of snapshots (data), v is effective search area, and D is GROUP (not animal) density
        logDgrp[z] <- logDgrp.hyper + logDgrp.err[z]    # expected group density on log scale (uncorrected for g0)
        Dgrp[z] <- exp(logDgrp[z])                      # expected group density on real scale (uncorrected)
        logDgrp.err[z] ~ dnorm(0, logDgrp.tau)          # random DASBR effect for group density (uncorrected)
        Dper1000[z] <- Dgrp[z]*s/g0 * 1000              # animal density, g0 corrected, and rescaled to 'per 1000km2' (note: g0 is the availability corrxn, see below)
      } # Z


      ##### PRIORS FOR DENSITY MODEL #####
      logDgrp.hyper ~ dnorm(0,0.001)              # hyper-parameter for group density, log scale
      logDgrp.tau <- 1/(logDgrp.sig*logDgrp.sig)
      logDgrp.sig ~ dunif(0,4)                    # random effect variance for group density, log scale


      ##### SOME MODEL CHECKS #####
      mean.of.Errs <- mean(logDgrp.err[])         # mean of the estimated random effects (on log scale), should be close to 0
      SD.of.Errs <- sd(logDgrp.err[])             # sd of the estimated random effects, should be similar to logDgrp.sig
      mean.of.logDgrp <- mean(logDgrp[])          # mean of the log-densities across DASBRs, should be similar to logDgrp.hyper


      ##### DERIVED DENSITY ESTIMATE FOR STUDY AREA #####
      Dgrp.mean <- sum(Dgrp[])/n.sites  # uncorrected mean density of *groups* in study area (arithmetic mean across DASBR-specific densities)
      D.mean <- Dgrp.mean * s / g0      # animal density = uncorrected group density * mean group size * availability correction
      D.mean.per1000 <- D.mean*1000     # Animal density per 1000 km2
      N <- D.mean * A                   # Population size (A is study area size)
      # Note: In this code, g0 is the availability correction, although the term g0 is used differently in our paper (we call the availability corrxn something else in the paper)


      ##### AN ALTERNATE ESTIMATOR FOR DENSITY - DO NOT USE #####
      Dgrp.mean.II <- exp(logDgrp.hyper + logDgrp.sig*logDgrp.sig/2)                          # group density, uncorrected
      D.mean.per1000.II <- exp(logDgrp.hyper + logDgrp.sig*logDgrp.sig/2) * s / g0 * 1000     # animal density, per 1000km2
      N.II <- D.mean.per1000.II * A/1000                                                      # population size
      # Note: DO NOT USE THIS OUTPUT
      # This estimator seems theoretrically reasonable to me, but results are clearly biased (too high) (due to extreme skew of lognormal?, or due to correlation in the paramemters?)
      # The estimator is the simple formula for finding the mean of a lognormal distribution given lognormal parameters, i.e., mean = exp(mu + sig2/2)
      # A future question is to understand why Dgrp.mean and Dgrp.mean.II don't agree with each other


      ##### ANCILLARY INPUTS TO DENSITY MODELS ABOVE #####
      s ~ dnorm(gsizeMean, tau.s)           # group size.  Mean s comes from Barlow 2016 report
        tau.s <- 1/(gsizeSd * gsizeSd)
        # sd.s <- gsizeSd

      # from EH analysis... number of minutes that a beaker is available for detection...
      timeAvail ~ dnorm(availMean, tau.timeAvail)
        tau.timeAvail <- 1/(sd.timeAvail * sd.timeAvail)
        sd.timeAvail <- availSd
      # Note, the 17.23 = estimate of 15.23 mean avail time + 2 min (for snapshot duration).
      # See Jay's paper on snapshot length for explanation of adding the duration (Barlow 2021, in press)

      # effective g0, the proportion of time between start of successive dives that an animal is available to detection (clicking and oriented properly to hydrophone)
      g0 <- timeAvail / deepDivePeriod
      deepDivePeriod ~ dnorm(ddpMean,ddp.tau)
      ddp.tau <- 1/(ddpSd * ddpSd)
      # ddp.sd <- ddpSd
      # deep dive period (min) is from email from Jay, dated July 23, 2020

    }	# end model

", fill=TRUE, file=modelFilename)
################################################################################################