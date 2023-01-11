#2021 April 8, Jeff Moore
# posterior predictive check on sightings

##### Libraries #####
library(Hmisc)

##### Set working directory #####
setwd("C:/jeff/NOAA/PASCAL/analysis/finalAnalysisDataAndCode-March2021/DensityEstimation/")

load("out/MCMCsamples_abundanceEstimates_4kTrunc.Rdata")  # load the MCMC object for the abudance analysis
ppd.n <- matrix(NA,out$n.sims,ncol(out$sims.list$mu))     # posterior predictive distribution for n, the number of detections on each DASBR
for(j in 1:ncol(out$sims.list$mu)){                       # loop through DASBRs
    ppd.n[,j] <- rpois(out$n.sims, out$sims.list$mu[,j])  # generate Poisson values for n (the number of observed sightings), given mu (the expected number of sightings for each DASBR)
  }

# summarize simulations
n = detPerStation[,2]  # from bugs data object (the real data)
ppd.means = apply(out$sims.list$mu,2,mean)  # mean of the simulated Poisson values (at each DASBR)

# plot observed vs. mean of the predicted
plot(1:21,ppd.means, col="red", xlab="station", ylab="n per station",cex=1.5)  # posterior means for mu
points(1:21,n, col="blue",pch=19)  # observed n
legend(x=1,y=100,legend=c("posterior means for mu","observed n"),pch=c(21,19),col=c("red","blue"))

# save (add code for this)