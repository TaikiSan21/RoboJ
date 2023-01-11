# Optim gives param estimates
# Gam is used as function for optim
# Gam is fit on likelihood DF that has all comb of params in acceptable range
# Likelihood df comes from prob density matrix - plugging indetected angles for each param combo
# Prob density matrix comes from simulating over each param combo
# Simulation does uniform area Horiz Range + Depth Distribution to get
#     Simulated Slant Range which is fed into CHN to get det Prob
#     Theses slant ranges are correlated with different angles, so we sum
#     probability over all slant ranges in an angle
#     Params only go into CHN
eventSummary <- readRDS('../Data/RoboJ/Test/eventSummary.rds')
eventSummary <- eventSummary[eventSummary$Station <= 22, ]
testAngle <- eventSummary$MeanAngle[eventSummary$eventType == 'ZC' &
                                        eventSummary$nClicks > 2]
testAngle <- testAngle[!is.na(testAngle) & testAngle < 74.53]
testAngle <- read.csv('../Data/RoboJ/Test/ZcObsDetAngles.csv', stringsAsFactors = F)[[2]]

testCHN <- createSimLikeFun(DepthDistr='log-normal', Model='C_HN', nSim=1e6)(testAngle)
testless <- createSimLikeFun(DepthDistr='log-normal', Model='HN', nSim=1e7)(testAngle[runif(length(testAngle)) < .6])
testo <- createSimLikeFun(DepthDistr='log-normal', Model='C_HN', nSim=1e7)(testAngle)
oldtest <- oldSimLikeFun(testAngle, DepthDistr='log-normal', Model='C_HN', nSim=1e7)
topt <- optim(par=c(1e3,1),fn=testo, method='Nelder-Mead', hessian=FALSE, control=list(fnscale=-1, parscale=c(2e3, 15)))
topt$counts
library(microbenchmark)
library(optimx)
microbenchmark(
    nohess = optim(par=c(1e3,1),fn=testo, method='Nelder-Mead', hessian=FALSE, control=list(fnscale=-1, parscale=c(2e3, 15))),
    old = optim(par=c(1e3, 1), fn=testo, method='BFGS', hessian=FALSE, control=list(fnscale=-1, parscale=c(2e3, 15))),
    times=5)
topt
tf <- testo(topt$par, F)
df <- data.frame(Angle=tf$Angle, SRange=tf$SRange) %>%
    mutate(iAngle = factor(ceiling(Angle)))
plot(tf$AngleDensity, type='l')
library(ggplot2)
library(ggridges)
ggplot(df) +
    geom_density_ridges_gradient(aes(x=SRange, y=factor(iAngle), fill=..x..)) +
    scale_fill_viridis() +
    xlim(0, 4e3)

lld <- GetLikelihoodData(testAngle, Model='C_HN', ProbDensityMatrix = pd, Range=Range, MaxPar2index = 100)
mx <- which.max(lld$LogLike)
lld$Param1[mx]
lld$Param2[mx]
testCHN(c(lld$Param1[mx], lld$Param2[mx]))
lld$LogLike[mx]

hm7 <- optim(par=c(1000, 1), fn=testCHN, method='Nelder-Mead', hessian=TRUE, control=list(fnscale=-1))
hm6 <- optim(par=c(1000, 1), fn=testCHN, method='Nelder-Mead', hessian=TRUE, control=list(fnscale=-1))
hmHN <- optim(par=1000, fn=testHN, method='Brent', hessian=TRUE, control=list(fnscale=-1), upper=4e3, lower=500)


p1s=seq(from=1e3, to=3e3, by=10)
p2s=seq(from=-1.7, to=28, by=.3)
jayCHN <- createSimLikeFun(detAngle = jayAngle, DepthDistr='log-normal', Model='C_HN', nSim=1e7)
parrFast <- makeLLArr(p1s=p1s, p2s=p2s, FUN=jayCHN)
lld <- data.frame(Param1=rep(p1s, length(p2s)), Param2=rep(p2s, each=length(p1s)), LogLike=as.vector(parrFast))
gamtest <- gam(formula= LogLike ~ te(Param1, Param2, k=25), data=lld)
otest <- optim(c(1e3, 1), fn=jayCHN, method='Nelder-Mead', hessian=T, control=list(fnscale=-1))
otest
gamoptim <- optim(par=c(3, 0), fn=gamShitter, method='Nelder-Mead', hessian=TRUE, control=list(fnscale=-1))
gamoptim

library(profvis)
profvis(makeLLArr(seq(from=1e3, to=3e3, length.out=20), seq(from=-2, to=28, length.out=20)))


gamShitter(c(1e3, .2))
gamoptim <- optim(par=c(1000, 1), fn=gamShitter, method='Nelder-Mead', hessian=TRUE, control=list(fnscale=-1))


doAllOptim(testAngle, opStart=c(1e3, 1), nSim=1e5)



manySim2 <- doManySim(testAngle, 100, 100, 1e5, 1e5)
manySim <- doManySim(testAngle, 50, 50, 1e6, 1e7)
manySim2 <- bind_rows(manySim2, manySim7)
manySim2 <- arrange(manySim2, nSample)
manySim2$col <- 'red'
manySim2$col[manySim2$nSample == 1e6] <- 'blue'
manySim2$col[manySim2$nSample == 1e7] <- 'green'

library(ggplot2)
ggplot(manySim7) +
    geom_density(aes(x=par1, col=as.character(nSample)))
ggplot(manySim7) +
    geom_density(aes(x=par2, col=as.character(nSample)))

plot(x=1:2, y=c(manySim2$par1[1], manySim2$par2[1]*200), col='black', type='l')
for(i in 1:nrow(manySim2)) {
    lines(x=c(1, 2), y=c(manySim2$par1[i], manySim2$par2[i]*200), col=ifelse(manySim2$nSample[i] > 2e5, 'blue', 'black'))
}

for(i in 1:nrow(manySim2)) {
    plotDetFun(c(manySim2$par1[i], manySim2$par2[i]), model = 'C_HN', add = i>1, col = manySim2$col[i])
}
plotDetFun(c(2312, -.442), model='C_HN', add=TRUE, col='green')

jayData <- read.csv('../Data/RoboJ/Jeffs Files/DensityEstimation/ZcEventFilesWithAngles_drifts_1-22 - CorrectedAngles v2.csv', stringsAsFactors = FALSE)
jayData$MeanAngle <- 180 - jayData$MeanAngle
jayAngle <- jayData$MeanAngle

simJay <- doManySim(jayAngle, 1e2, 1e2)
plot(x=1:2, y=c(simJay$par1[1], simJay$par2[1]*200), col='black', type='l')
for(i in 1:nrow(simJay)) {
    lines(x=c(1, 2), y=c(simJay$par1[i], simJay$par2[i]*200), col=ifelse(simJay$nSample[i] > 2e5, 'blue', 'black'))
}

for(i in 1:nrow(simJay)) {
    plotDetFun(c(simJay$par1[i], simJay$par2[i]), model = 'C_HN', add = i>1, col = ifelse(simJay$nSample[i] > 2e5, 'blue', 'red'))
}
# 1919 .609, 1661 1.34, 1427 2.22
plotDetFun(c(1919, .61), add=T, col='darkgreen', lwd=3)
plotDetFun(c(1661, 1.34), col='darkgreen', lwd=3, add=T)
plotDetFun(c(1427, 2.22), add=T, col='green', lwd=3)


# 1427 2.22
predict.gam(detEstJay1e6$`All Recorders`$gam, newdata=data.frame(Param1=c(1427, 1509, 1577), Param2=c(2.22, 1.87, 1.61)))
str(detEstJay1e6$`All Recorders`$lld)
mat <- matrix(detEstJay1e6$`All Recorders`$lld$LogLike, ncol=201)
image(t(mat[5:40, 10:130]), col=hcl.colors(50, "YlOrRd", rev = TRUE))
