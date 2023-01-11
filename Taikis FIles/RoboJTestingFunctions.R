# roboj testing functions



estOptimDetFunction <- function(detAngle, nSim=1e6, Model='C_HN', opStart=c(1e3, 1), verbose=TRUE) {
    CHNfun <- createSimLikeFun(nSim=nSim, Model=Model)(detAngle)
    opResult <- if(Model == 'C_HN') {
        optim(opStart, fn=CHNfun, method='Nelder-Mead', hessian=TRUE, control=list(fnscale=-1))
    } else {
        optim(opStart[1], fn=CHNfun, method='Brent', upper=4e3, lower=500, hessian=TRUE, control=list(fnscale=-1))
    }
    if(verbose) {
        cat('Params :', opResult$par, ' Likelihood: ', opResult$value)
    }
    result <- list(par1=opResult$par[1], par2=opResult$par[2], value=opResult$value)
    vals <- CHNfun(opResult$par, like=FALSE)
    result$EDR <- vals$EDR
    result$EDR_39X <- vals$EDR_39X
    result
}

doManySim <- function(testAngle, n1=1e2, n2=1e2, samp1=1e5, samp2=1e7) {
    out1 <- vector('list', length=n1)
    pb <- txtProgressBar(min=0, max=n1, style=3)
    for(i in seq_along(out1)) {
        thisOut <- try(estOptimDetFunction(testAngle, opStart=c(1e3, 0), nSim=samp1, verbose=FALSE))
        setTxtProgressBar(pb, value=i)
        if(inherits(thisOut, 'try-error')) next
        out1[[i]] <- thisOut
    }
    out1 <- bind_rows(out1)
    out1$nSample <- samp1
    out2 <- vector('list', length=n2)
    cat('\n')
    pb <- txtProgressBar(min=0, max=n2, style=3)
    for(i in seq_along(out2)) {
        thisOut <- try(estOptimDetFunction(testAngle, opStart=c(1e3, 0), nSim=samp2, verbose=FALSE))
        setTxtProgressBar(pb, value=i)
        if(inherits(thisOut, 'try-error')) next
        out2[[i]] <- thisOut
    }
    out2 <- bind_rows(out2)
    out2$nSample <- samp2
    bind_rows(out1, out2)
}

gamShitter<- function(par, Param1Limits=c(1e3, 3e3), Param2Limits=c(-1.7, 28)) {
    Param1<- par[1]*1e3
    if (length(par) == 1) {Param2=0} else {Param2= par[2]}
    if ((Param1 < Param1Limits[1])|(Param1 > Param1Limits[2])|
        (Param2 < Param2Limits[1])|(Param2 > Param2Limits[2])) {
        LogL<- -exp(100)         # penalize LogL if parameters are outside fitted range of gam_logLike
    } else {
        LogL<- predict.gam(gamtest,newdata=data.frame(Param1=Param1,Param2=Param2))
    }
    return(LogL)
}

makeLLArr <- function(p1s=seq(from=1e3, to=3e3, by=10), p2s=seq(from=-2, to=28, by=.3), FUN=testCHN, progress=TRUE) {
    if(progress) {
        pb <- txtProgressBar(min=0, max=length(p1s)*length(p2s), style=3)
        ix <- 1
    }
    pArr <- matrix(NA, nrow=length(p1s), ncol=length(p2s))
    for(i in seq_along(p1s)) {
        for(j in seq_along(p2s)) {
            pArr[i, j] <- FUN(c(p1s[i], p2s[j]))
            if(progress) {
                setTxtProgressBar(pb, value=ix)
                ix <- ix + 1
            }
        }
    }
    pArr
}

denSimPlots <- function(angle, nSim=c(1e7, 2e7, 3e7)) {
    denOut <- vector('list', length=length(nSim))
    dfOut <- denOut
    for(i in seq_along(nSim)) {
        simFun <- createSimLikeFun(DepthDistr='log-normal', Model='C_HN', nSim=nSim[i])(angle)
        tf <- simFun(c(1e3, 1), F)
        dfOut[[i]] <- data.frame(Angle=factor(ceiling(tf$Angle)), SRange=tf$SRange, nSim=nSim[i])
        denOut[[i]] <- list(Angle=1:91, Density=tf$AngleDensity, nSim=nSim[i])
    }
    denOut <- bind_rows(denOut)
    dfOut <- bind_rows(dfOut)
    gDen <- ggplot(denOut, aes(x=Angle, y=Density, col=factor(nSim))) + geom_line(size=1)
    # gRidge <- ggplot(dfOut) +
    #     geom_density_ridges_gradient(aes(x=SRange, y=Angle, fill=..x..)) +
    #     scale_fill_viridis() +
    #     xlim(0, 4e3) +
    #     facet_wrap(~nSim)
    gRidge <- NA
    list(gDen=gDen, gRidge=gRidge, denOut=denOut, dfOut=dfOut)
}
