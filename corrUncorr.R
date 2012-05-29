library(animation)

addLines <- function(X, Y, numLines, lineColors){
    if(numLines){
        for(i in 1:numLines){
            lines(x= X[(2*i -1):(2 * i)], y= Y[(2*i -1):(2 * i)],
                col= lineColors, lwd= 2)
        }
    }
}

corrUncorr <- function(sdx, sdy, rho, numPoints= 100,
    pointColors= c("violet"), lineColors= c("red"),
    numSteps= 100, numLines= 2, aniInt= 0.075){
    # begun by Laurie Samuels, May 2012
    # in progress
    # no error handling!

    covXY <- rho * sdx * sdy
    sdMax <- max(sdx, sdy)

    X0 <- rnorm(numPoints)
    Y0 <- rnorm(numPoints)

    # let's get them centered exactly at 0, just because
    X0 <- X0 - mean(X0)
    Y0 <- Y0 - mean(Y0)
    M0 <- rbind(X0, Y0)

    # 1. Make M0 perfectly uncorrelated N(0,1)
    covMat0 <- matrix(c(sd(X0)^2, cov(X0, Y0), cov(X0, Y0), sd(Y0)^2), ncol= 2)
    # Cholesky transformation:
    # http://blogs.sas.com/content/iml/2012/02/08/use-the-cholesky-transformation-to-correlate-and-uncorrelate-variables/
    cholU0 <- chol(covMat0)
    cholL0 <- t(cholU0)
    cholLInv0 <- solve(cholL0)

    M1 <- cholLInv0 %*% M0

    # 2. Make data w/ specified sd's & cor from the uncorrelated N(0, 1) data
    covMat2 <- matrix(c(sdx^2, covXY, covXY, sdy^2), ncol= 2)
    cholU2 <- chol(covMat2)
    cholL2 <- t(cholU2)
    M2 <- cholL2 %*% M1

    X2 <- M2[1, ]
    Y2 <- M2[2, ]

    # 3. Now let's standardize them to have the same sd, but keep the cor
    # incrementally, for the animation
    incrStand <- matrix(nrow= 2 * numSteps, ncol= numPoints)
    sdIncr <- (sdMax - c(sd(X2), sd(Y2))) / numSteps
    for(i in 1:numSteps){
        # Xs:  odd-numbered rows
        incrStand[(2 * i) - 1, ] <- X2 * (1 + (i * sdIncr[1])/sd(X2))
        # Ys: even-numbered rows
        incrStand[2 * i, ] <- Y2 * (1 + (i * sdIncr[2])/sd(Y2))
    }

    # 4. Now let's get rid of the cor
    # incrementally, as above
    incrUncorr <- matrix(nrow= 2 * numSteps, ncol= numPoints)
    for(i in 1:numSteps){
        myCov <- (rho * sdMax^2) * (numSteps - i) / numSteps
        covMat <- matrix(c(sdMax^2, myCov, myCov, sdMax^2), ncol= 2)
        cholU <- chol(covMat)
        cholL <- t(cholU)
        newXY <- cholL %*% M1
        incrUncorr[((2 * i) - 1), ] <- newXY[1, ]
        incrUncorr[(2 * i), ] <- newXY[2, ]
    }

    # The animation

    # Find axis limits for graphs
    maxval <- abs(max(c(X2, Y2, 
        incrStand[nrow(incrStand) - 1, ], 
        incrStand[nrow(incrStand), ],
        incrUncorr[nrow(incrUncorr) - 1, ],
        incrUncorr[nrow(incrUncorr), ]))) 

    oopt = ani.options(interval = aniInt)

    # the data with the specified sds & cor
    plot(X2, Y2, col= pointColors, pch= 19, xlim= c(-maxval, maxval),
        ylim= c(-maxval, maxval), xlab= "", ylab= "")
    addLines(X2, Y2, numLines, lineColors)
    for(i in 1:10){
        ani.pause() # length of pause= ani.options(interval)
    }

    # incremental standardization of sds
    for(i in 1:numSteps){
        plot(incrStand[(2 * i) - 1, ], incrStand[(2 * i), ],
            col= pointColors, pch= 19, xlim= c(-maxval, maxval),
            ylim= c(-maxval, maxval), xlab= "", ylab= "")
        addLines(incrStand[(2 * i) - 1, ], incrStand[(2 * i), ],
            numLines, lineColors)
        ani.pause()  
    }
    for(i in 1:10){
        ani.pause()  
    }

    # incremental decorrelation
    for(i in 1:numSteps){
        plot(incrUncorr[(2 * i) - 1, ], incrUncorr[(2 * i), ],
            col= pointColors, pch= 19, xlim= c(-maxval, maxval),
            ylim= c(-maxval, maxval), xlab= "", ylab= "")
        addLines(incrUncorr[(2 * i) - 1, ], incrUncorr[(2 * i), ],
            numLines, lineColors)
        ani.pause()  
    }

    ## restore the options
    ani.options(oopt)
}


corrUncorr(3, 4, .7)


#######################################
# old code


#plot(X3, Y3, col= Colors, pch= 19, xlim= c(-maxval, maxval),
#    ylim= c(-maxval, maxval), xlab= "", ylab= "")
#ani.pause()  ## pause for a while (’interval’)


# Instead of doing it this way, do a series of
# in-between graphs-- it's too hard to see what's happening with the
# individual points
#Xs <- cbind(X2, X3, X4)
#Ys <- cbind(Y2, Y3, Y4)
#Colrs <- rep(Colors, each= numPoints)
#for(i in 1:(numPoints * 2 + 1)){
#    myIndices <- i:(i + numPoints - 1)
#    plot(Xs[myIndices], Ys[myIndices], col= Colrs[myIndices],
#        pch= 19, xlim= c(-maxval, maxval),
#        ylim= c(-maxval, maxval), xlab= "", ylab= "")
#    ani.pause()  ## pause for a while (’interval’)
#}

## with specified sds & cor
#plot(X2, Y2, col= Colors, pch= 19, xlim= c(-maxval, maxval),
#    ylim= c(-maxval, maxval), xlab= "", ylab= "")
#ani.pause()  ## pause for a while (’interval’)
#
## sd = same for both, cor still there
#plot(X3, Y3, col= Colors, pch= 19, xlim= c(-maxval, maxval),
#    ylim= c(-maxval, maxval), xlab= "", ylab= "")
#ani.pause()  ## pause for a while (’interval’)
#
## sd = same for both, uncorrelated
#plot(X4, Y4, col= Colors, pch= 19, xlim= c(-maxval, maxval),
#    ylim= c(-maxval, maxval), xlab= "", ylab= "")
#ani.pause()  ## pause for a while (’interval’)

## restore the options
ani.options(oopt)

