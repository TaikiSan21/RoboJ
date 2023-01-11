# 2021 April 8, Jeff Moore
# called by EH-runBUGS.R
# OpenBUGS CJS model fit to DASBR encounter history data
# Model phi(.)p(.)

# Survival parameter, phi, can be interpreted as probability of remaining 'available' following initial detection
# 1 - phi is probability of becoming permanently unavailable
# p is probability of detection given that it remains available
# Given phi, we estimate the mean (and median) "availability time" in terms of the number 2-min intervals or number of minutes (intervals x 2)


##### give name to text file that will hold WinBUGS model code #####
# modelFilename="C:/jeff/NOAA/PASCAL/analysis/finalAnalysisDataAndCode-March2021/AvailabilityEstimation/ziphius_EHmodel.txt"
library(here)
modelFilename <- here('availability_estimation', 'ziphius-EHmodel.txt')
################################################################################################
##### Define the model in BUGS code and write to text file, which will be called by bugs() #####
cat("

  model{

    for(i in 1:nind){
      for(t in 2:nocc){              # n.occ = number of occasions, including initial capture occasion
        phi.it[i,t] <- phi
        p.it[i,t] <- p * dcyc[i,t]     # p = p0 when dcyc[i,t]=1, and equals 0 when dcyc[i,t]=0 (detector is off).  Note, dcyc is a data object.
      }  # t
    }  # i

    ### PRIORS
    phi ~ dunif(0,1)
    p ~ dunif(0,1)

    ### LIKELIHOOD
    for(i in 1:nind){
      z[i,1] <- eh[i,1]                           # z indicates whether indiv i is available to sampling; this is 1 in the 1st occasion for all individuals since the encounter history is conditioned on the initial detection
      for(t in 2:nocc){
        mu.avail[i,t] <- z[i,t-1] * phi.it[i,t]   # mu.avail = prob i survived to t, given it was in sampling population at t-1 (z[i,t-1]=1)
        z[i,t] ~ dbern(mu.avail[i,t])             # index whether animal is still available to detection during capture occasion t
        mu.det[i,t] <- z[i,t] * p.it[i,t]            # mu.det =  prob of being detected (0 if animal no longer avaiable; p otherwise)
        eh[i,t] ~ dbern(mu.det[i,t])              # likelihood for the data (eh[i,t] is the data object)
      } # t
    } # i

    ### DERIVED PARAMETERS

      # Median availability time (in intervals and minutes)
      med.availtime.min <- (log(0.5)/log(phi)) * 2     # time in minutes for which 50% of animals are available to detection
      med.availtime.ints <- med.availtime.min / 2      # number of intervals (inluding 1st) for which 50% of animals are available to detection

      # Mean availablity time (in intervals and minutes)
      mean.availtime.min <- mean.availtime.ints * 2
      for(i in 1:nocc){
        y[i] <- d[i] * ((i-1) + 0.5)             # y(i) = d(i) * (i + 0.5), where d(i) is the proportion that die between i and i+1
        d[i] <-  pow(phi, (i-1)) * (1 - phi)     # d(i) = l(i) * q(i), where l(i) is survivorship to age i; q(i) = 1-phi(i)
      }
      mean.availtime.ints <- sum(y[])

  } # end model

", fill=TRUE, file=modelFilename)     # close cat() function
################################################################################################

