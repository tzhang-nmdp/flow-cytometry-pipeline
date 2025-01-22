# FlowQ test
require(pheatmap)
require(RColorBrewer)
require(scales)
require(lattice)
require(reshape2)
require(tabplot)
require(googleVis)
require(tidyverse)
require(CATALYST)
require(knitr)
require(ggcyto)
require(shinydashboard)

require(cytoqc)
require(flowAI)
require(PeacoQC)
require(flowCore)

require(FlowSOM)
require(spade)
setwd("C:/Users/tzhang/Desktop/Project/Graft composition")

## initialize the score data frame for a flow set
flowQA_firstScoreInit <- function(set, time){
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    if(missing(time))
        time <- findTimeChannel(set)
    
    fsApply(set, function(x) data.frame(cell_ID = 1:nrow(x), 
                                        time = x[,time], 
                                        score = rep(0, nrow(x))), 
            simplify=FALSE, use.exprs=TRUE)
}

flowQA_scoreInit <- function(set){
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    
    fsApply(set, function(x) data.frame(cell_ID = 1:nrow(x), score = rep(0, nrow(x))), 
            simplify=FALSE, use.exprs=TRUE)
}

## update the score frame
flowQA_scoreUpdate <- function(fsScore, nScore){
    if(length(fsScore) != length(nScore))
        stop("input must be of equal length!\n")
    for(i in names(fsScore)){
        fsScore[[i]] <- merge(fsScore[[i]], nScore[[i]])
    }
    return(fsScore)        
}

## summarize the scores, with sum score of scaling to 0-1
flowQA_scoreSummary <- function(fsScore){
    
    scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
    for (i in 1:length(fsScore)){
        s <- fsScore[[i]]
        standScores <- apply(subset(s, select = -c(cell_ID, score, time)), 2, scale01)
        s$score <- apply(standScores, 1, sum)
        fsScore[[i]] <- s
    } 
    return(fsScore)
}

## Cell number check
flowQA_cellnum <- function(set) { 
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    
    cellNumbers <- as.numeric(fsApply(set, nrow))
    cellNumbers[cellNumbers==1] <- NA
    cnframe <- data.frame(sampleName = sampleNames(set), cellNumber = cellNumbers)
}     


#' Margin events check
#' Check the percentage of margin events for selected markers and plot in heatmap
flowQA_marginevents <- function(set, channels=NULL, side="both", tol = NULL) {
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    if(is.null(channels)){
        parms <- setdiff(colnames(set[[1]]), c("time", "Time"))
    }else{
        if(!all(channels %in% colnames(set[[1]])))
            stop("Invalid channel(s)")
        parms <- channels
    }
    if(is.null(tol))
        tol <- -.Machine$double.eps
    
    nmarkers <- length(parms)
    nsamples <- length(set)
    perc <- matrix(nrow = nmarkers, ncol = nsamples)
    for(i in 1:nmarkers){
        marker <- parms[i]
        for(j in 1:nsamples){
            ff <- set[[j]]
            dat <- exprs(ff)[,marker]
            ranges <- range(ff)[,marker]
            if(length(dat) ==0){
                res <- NULL
            }else{
                res <- switch(side,
                              both={(dat > (ranges[1]-tol)) & (dat<(ranges[2]+tol))},
                              upper = dat < (ranges[2]+tol), 
                              lower = dat > (ranges[1]-tol))
                res[is.na(res)] <- TRUE        
            }   
            res[is.na(res)] <- TRUE
            trueCount <- sum(res==T)
            count <- length(dat)
            perc[i,j] <- 1- trueCount/count 
        }
    }
    
    colnames(perc) <- sampleNames(set)
    rownames(perc) <- parms
    
    return(perc)   
}

## score = percentage of margin channels of each cell
marginEventsScore <- function(set, parms, side, tol){
    meScore <- flowQA_scoreInit(set)
    nmarkers <- length(parms)
    for(i in 1:nmarkers){
        marker <- parms[i]
        for(j in 1:length(set)){
            ff <- set[[j]]
            dat <- exprs(ff)[,marker]
            ranges <- range(ff)[,marker]
            s <- meScore[[j]]
            res <- switch(side,
                          both={(dat > (ranges[1]-tol)) & (dat<(ranges[2]+tol))},
                          upper = dat < (ranges[2]+tol), 
                          lower = dat > (ranges[1]-tol))
            res[is.na(res)] <- TRUE
            s$score <- s$score + as.numeric(!res)
            meScore[[j]] <- s   
        }
    }
    
    meScore <- lapply(meScore, 
                      function(x) {
                          names(x) <- c("cell_ID", "meScore")
                          x$meScore <- x$meScore / nmarkers # percentage of markers
                          x  })     
}
    
## Detect disturbances in the flow of cells over time
flowQA_timeflow <- function(set, time, binSize) {
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    if(missing(time))
        time <- findTimeChannel(set)
    if(missing(binSize) || is.null(binSize) || is.na(binSize)) 
        binSize <- min(max(1, floor(median(fsApply(set, nrow)/100))), 500)
    
    timeFlowData <- fsApply(set[fsApply(set,nrow) !=1], timeFlowCheck, time=time,
                            binSize=binSize, use.exprs=FALSE, simplify=FALSE) 
    
    timeflow_res <- list(timeFlowData = timeFlowData, binSize = binSize)
    return(timeflow_res)
}    

timeFlowCheck <- function(x, time, binSize) {
    xx <- sort(exprs(x)[, time])   # time channel
    lenx <- length(xx)             # num of time ticks
    ux <- unique(xx)               # num of unique time tick
    nrBins <- floor(lenx/binSize)  # num of bins
    
    if(length(ux) < nrBins || lenx < binSize){
        warning("Improper bin size, auto change bin size", call.=FALSE)
        binSize <- min(max(1, lenx/100), 500)
        nrBins <- floor(lenx/binSize)
    }
    
    tbins <- seq(min(xx), max(xx), len=nrBins + 1)    # time bins
    tbCounts <- hist(xx, tbins, plot = FALSE)$counts  # number of events per time bin
    expEv <- lenx/(nrBins)   ##median(tbCounts)       # expected number of events per bin
    return(list(frequencies=cbind(tbins[-1], tbCounts), expFrequency=expEv, bins=nrBins))
}

# ## score
# fsScore <- timeFlowScore(fsScore, tfScore, set, time, binSize, varCut)
timeFlowScore <- function(set, binSize, varCut) {
    tfScore <- flowQA_scoreInit(set)
    time <- findTimeChannel(set)
    
    for (i in 1:length(set)){
        x <- set[[i]]
        s <- tfScore[[i]]
        xod <- order(exprs(x)[, time])
        s <- s[xod, ]
        xx <- exprs(x)[, time][xod]    # time channel
        lenx <- length(xx)             # num of time ticks
        ux <- unique(xx)               # num of unique time tick
        nrBins <- floor(lenx/binSize)  # num of bins
        
        if(length(ux) < nrBins || lenx < binSize){
            warning("Improper bin size, auto change bin size", call.=FALSE)
            binSize <- min(max(1, lenx/100), 500)
            nrBins <- floor(lenx/binSize)
        }
        
        tbins <- seq(min(xx), max(xx), len=nrBins + 1)    # time bins
        tbCounts <- hist(xx, tbins, plot = FALSE)$counts  # number of events per time bin
        expEv <- lenx/(nrBins)  ## median(tbCounts)        # expected number of events per bin
        
        binCheck <- (tbCounts - expEv) / expEv / varCut # bin score based on variation to mean
        #binCheck[binCheck < 0 & binCheck > -0.8] <- 0.1  ### speed decreas is OK in certain range
        binCheck <- abs(binCheck)
        cellCheck <- unlist(sapply(1:nrBins, function(x) rep(binCheck[x], tbCounts[x])))
        s$score <- cellCheck  # cell score
        tfScore[[i]] <- s
    }
    
    tfScore <- lapply(tfScore, function(x) {
        names(x) <- c("cell_ID", "tfScore")
        x  })
}


## Plot frequency values for each flowFrame
## y, list of outputs from timeFlowCheck with name of each flow frame
## varCut, The cutoff in the adjusted variance(standY) represented by the width of rectangule
timeFlowPlot <- function(y, timeThres, varCut = 1, main="Time Flow Frequency", 
                         col = NULL, ylab = NULL, lwd = 2, ...) {
    
    if(missing(ylab) | is.null(ylab)) 
        ylab <- truncNames(names(y))
    if(missing(col) | is.null(col)){
        colp <- brewer.pal(8, "Dark2")
        col <- colorRampPalette(colp)(length(y))
    }else{
        if(length(col)!=1 || length(col)!=length(y))
            stop("'col' must be color vector of length 1 or same length as the flowSet")
    }
    
    for(i in names(y)){
        freq <- y[[i]]$frequencies
        range <- timeThres[[i]]
        time <- freq[,1]
        badCellid <- time < range[1] | time > range[2]
        freq[badCellid,2] <- 0
        y[[i]]$frequencies <- freq
    }
    
    stX <- lapply(y, function(x) x$frequencies[,1])
    stY <- lapply(y, function(x) x$frequencies[,2] / (median(x$frequencies[,2][x$frequencies[,2] > 0])) - 1)
    
    stacks <- ((length(y):1)-1) * max(diff(range(stY)), varCut*2)*1.02
    stYY <- mapply(function(x,s) x+s, stY, stacks, SIMPLIFY=FALSE)
    if(!is.list(stYY)) stYY <- list(stYY)
    
    ## draw the plot frame
    opar <- par(c("mar", "mgp", "mfcol", "mfrow", "las"))
    on.exit(par(opar))
    layout(matrix(1:2), heights=c(3.5, 1))
    par(mar=c(1,5,3,3), mgp=c(2,0.5,0), las=1)
    xlim <- c(0, max(unlist(stX)))
    ylim <- range(c(stYY), stacks+varCut+0.2, stacks-varCut-0.2)
    plot(stX[[1]], stYY[[1]], xlab="", ylab="", type="n", xaxt="n", 
         lwd=lwd, xlim=xlim, ylim=ylim, main=main, yaxt="n", ...)
    xl <- par("usr")[1:2]
    xl <- xl + c(1,-1)*(diff(xl)*0.01)
    axis(2, stacks, ylab, cex.axis=0.6)
    for(j in 1:length(y)){
        rect(xl[1], stacks[j]-varCut, xl[2], stacks[j]+varCut, 
             col=desat("lightgray", by=50), border=NA)
        lines(stX[[j]], stYY[[j]], col=col[j], lwd=lwd, ...) 
    }
    
    ## compute QA score and spearman rank correlation coefficient
    qaScore <- sapply(stY, function(z) sum(abs(z) > varCut)/length(z))*100
    #srho <- mapply(cor, stX, stY, method="spearman", use="pairwise") 
    
    # a legend indicating the problematic samples
    par(mar=c(5,3,0,3), las=2)
    top <- 2  # top value cut of the barplot
    barplot(qaScore, axes=FALSE, col=col, cex.names=0.6, names.arg = truncNames(names(qaScore)),
            ylim=c(0, min(c(top, max(qaScore, na.rm=TRUE)))), border=col, space=0.4)
    wh <- which(qaScore >= top)
    points((wh+wh*0.4)-0.5, rep(top-(top/12), length(wh)), pch=17, col="white", cex=0.5)
}


## function qaProcess.time
flowQA_timeline <- function(set, channels=NULL, time, varCut=1, binSize = NULL) {
    ## some sanity checking
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    if(is.null(channels)){
        parms <- setdiff(colnames(set[[1]]), c("time", "Time"))
    }else{
        if(!all(channels %in% colnames(set[[1]])))
            stop("Invalid channel(s)")
        parms <- channels
    }
    
    if(missing(time) || is.null(time))
        time <- findTimeChannel(set)
    
    if(!is.numeric(varCut) || length(varCut)!=1)
        stop("'varCut' must be numeric scalar")
    
    if(missing(binSize) || is.null(binSize) || is.na(binSize)) 
        binSize <- min(max(1, floor(median(fsApply(set, nrow)/100))), 500)
    
    tlScore <- flowQA_scoreInit(set)
    
    ## score
    tlScore <- timeLineScore(tlScore, set, parms, time, binSize, varCut)
    
    ## create timeline plot for each marker
    timeLineData <- timeLineCheck(set, parms, time, binSize, varCut)
    
    timeline_res <- list(timeLineData = timeLineData, varCut = varCut, 
                         binSize = binSize, tlScore = tlScore)
    return(timeline_res)
} 


timeLineScore <- function(tlScore, set, parms, time, binSize, varCut){
    nmarkers <- length(parms)
    for (i in 1:length(set)){
        x <- set[[i]]
        s <- tlScore[[i]]
        exp <- exprs(x)
        xx <- exp[, time]
        ord <- order(xx)
        xx <- xx[ord]                   # time data
        yy <- exp[, parms][ord, ]       # channels data
        s <- s[ord, ]              
        lenx <- length(xx)              # num of time ticks
        ux <- unique(xx)                # num of unique time tick
        nrBins <- floor(lenx/binSize)   # num of bins
        
        ## channel expression median (locM & varM) value agains time bin
        if(length(ux) < nrBins || lenx < binSize){
            warning("Improper bin size, auto change bin size to unique time tick", call.=FALSE)
            tmpy <- split(as.data.frame(yy), xx)
            splitScore <- split(s, xx)
        }else{
            cf <- c(rep(1:nrBins, each=binSize), rep(nrBins+1, lenx-nrBins*binSize))
            stopifnot(length(cf) == lenx)
            tmpy <- split(as.data.frame(yy),cf)
            splitScore <- split(s, cf)
        }
        
        yy <- t(sapply(tmpy, function(y) apply(y, 2, mean)))
        standy <- apply(yy, 2, function(z) abs(z-median(z)) / median(z)  / varCut)
        
        for(k in 1:length(splitScore)){
            sl <- nrow(splitScore[[k]])
            splitScore[[k]]$score <- rep(sum(standy[k,]), sl)
        }
        
        s <- do.call(rbind, splitScore)
        tlScore[[i]] <- s
    }
    
    tlScore <- lapply(tlScore, function(x) {
        names(x) <- c("cell_ID", "tlScore")
        # x$tlScore <- x$tlScore / nmarkers  # average
        x  })
    
    return(tlScore)         
}

## Modified from flowCore:::prepareSet
## binned <- timeCheck(x, param, time, bs, locM=median, varM=mad)
## median absolut varaition, mad = median(abs(x - median(x)))
#' @param x a flowFrame object 
#' @param parm channel names for the selected markers
#' @param time channel name of Time
#' @param binSize the size of bin
#' @param locM local median value (median value)
#' @param varM median variance value (median absolute deviation)
timeLineCheck <- function(set, parms, time, binSize, varCut=varCut) {        
    tlResults <- fsApply(set, function(x) {
        exp <- exprs(x)
        xx <- exp[, time]
        ord <- order(xx)
        xx <- xx[ord]                   # time data
        yy <- exp[, parms][ord, ]       # channels data
        lenx <- length(xx)              # num of time ticks
        ux <- unique(xx)                # num of unique time tick
        nrBins <- floor(lenx/binSize)   # num of bins
        
        ## channel expression median (locM & varM) value agains time bin
        if(length(ux) < nrBins || lenx < binSize){
            warning("Improper bin size, auto change bin size to unique time tick", call.=FALSE)
            tmpy <- split(as.data.frame(yy), xx)
            xx <- unique(xx)
            binSize <- 1
        }else{
            cf <- c(rep(1:nrBins, each=binSize), rep(nrBins+1, lenx-nrBins*binSize))
            stopifnot(length(cf) == lenx)
            tmpx <- split(xx,cf)
            tmpy <- split(as.data.frame(yy),cf)
            xx <- sapply(tmpx, mean)               
        }
        ## method for standardize the binned y value
        ## old: yy <- median(y); vy <- mad(y); stdy <- (yy - median(yy)) / mean(vy)
        yy <- t(sapply(tmpy, function(y) apply(y, 2, mean)))
        standy <- apply(yy, 2, function(z) (z-median(z)) / median(z) / varCut)
        return(list(res = cbind(time = xx,standy), bins = nrBins)) 
    }, simplify=FALSE)
    
    return(tlResults)
}

timeLinePlot <- function(x, timeThres = NULL, main = "Time Line Plot", channels = NULL, 
                         col = NULL, ylab=names(x), lwd =1, ...) {
    
    if(missing(col) | is.null(col)){
        colp <- brewer.pal(8, "Dark2")
        col <- colorRampPalette(colp)(max(sapply(x, ncol)))  #col <- rainbow(ncol(x[[1]])) 
    }else{
        if(length(col)!=1 || length(col)!=ncol(x[[1]]))
            stop("'col' must be color vector of same length as the markers")
    }
    if(is.null(channels)){
        parms <- setdiff(colnames(x[[1]]), c("time", "Time"))
    }else{
        if(!all(channels %in% colnames(x[[1]])))
            stop("Invalid channel(s)")
        parms <- channels
    }
    if(missing(ylab) | is.null(ylab))
        ylab <- truncNames(names(x))
    
    if(!(is.null(timeThres))){
        for(i in names(x)){
            res <- x[[i]]
            range <- timeThres[[i]]
            time <- res[,1]
            badCellid <- time < range[1] | time > range[2]
            x[[i]] <- res[!(badCellid), c("time",parms)]
        }
    }
    
    actualRange <- max(c(sapply(x, function(y) diff(range(y[,-1], na.rm=TRUE))), 2)) * 1.01
    stacks <- ((length(x):1)-1) * actualRange
    for (i in seq_len(length(x))){
        x[[i]][,-1] <- x[[i]][,-1] + stacks[i]    
    }
    xlim <- c(0, max(sapply(x, function(z) max(z[,1], na.rm=TRUE))))
    ylim <- range(c(sapply(x, function(z) range(z[,-1], na.rm=TRUE))))
    
    ## Create the plot, either one of the 4 possible types.
    opar <- par(c("mar", "mgp", "mfcol", "mfrow", "las"))
    on.exit(par(opar))
    layout(cbind(1,2), widths=c(5,1))
    par(mar=c(3,5,3,1), mgp=c(2,0.5,0), las=1)
    plot(x[[1]], xlab="", ylab="", type="n", xaxt="s", 
         lwd=lwd, xlim=xlim, ylim=ylim, main=main, yaxt="n", ...)
    axis(2, stacks, ylab, cex.axis=0.6)
    xl <- par("usr")[1:2]
    xl <- xl + c(1,-1)*(diff(xl)*0.01)
    for(j in 1:length(x)){
        rect(xl[1], stacks[j] - 1, xl[2], stacks[j] + 1, 
             col=desat("lightgray", by=20), border=NA)
        for (k in 2:ncol(x[[j]])){
            lines(x[[j]][, 1], x[[j]][, k], col=col[k], lwd=lwd) 
        }            
    }
    
    par(mar=c(2, 0, 3, 1))
    plot.new()
    legend('topleft', legend = colnames(x[[1]])[-1], col=col[-1], lty = 1, cex = 0.4)
}


#' tabplot of the expression value for selected sample
flowQA_tabplot <- function(set, sample = 1, channels = NULL, binSize = 500, 
                           from = 0, to = 100, scales = "auto"){
    snames <- sampleNames(set)
    if(is.numeric(sample)){
        if (sample > length(snames) || sample < 0)
            stop("wrong sample ID")
    } else if (is.character(sample)){
        if (!(sample %in% snames))
            stop("sample not in the x")
    }
    
    frame <- set[[sample]]
    time <- findTimeChannel(frame)
    
    if(is.null(channels)){
        parms <- colnames(frame)
    }else{
        if(!all(channels %in% colnames(frame)))
            stop("Invalid channel(s)")
        parms <- unique(c(time, channels))
    }
    
    dat <- as.data.frame(frame@exprs[ ,parms])
    nBins <- round(nrow(dat) / binSize)
    
    tableplot(dat, sortCol = 1, decreasing = FALSE, nBins = nBins, 
              scales = scales, from = 0, to = 100)
}


#' tabplot of the quality scores for selected sample
flowQA_scoreplot <- function(fsScore, sample = 1, timeThres = NULL, scoreThres = 3, nBins = 100, 
                             from = 0, to = 100, sortID = 1, scales = "lin"){
    snames <- names(fsScore)
    if(is.numeric(sample)){
        if (sample > length(snames) || sample < 0)
            stop("wrong sample ID")
    } else if (is.character(sample)){
        if (!(sample %in% snames))
            stop("sample not in the x")
    }
    dat <- fsScore[[sample]]
    if(!(is.null(timeThres))){
        range <- timeThres[[sample]]
        time <- dat[ ,"time"]
        badCellid <- time < range[1] | time > range[2]
        dat <- dat[!(badCellid), ]
        dat <- subset(dat, score <= scoreThres)
    }
    colOrder <- c(sortID, setdiff(1:ncol(dat), sortID))
    dat <- as.data.frame(dat[,colOrder])
    tableplot(dat, sortCol = 1, decreasing = FALSE, nBins = nBins, 
              from = 0, to = 100, scales = scales)
}

flowQA_scorecheck <- function(flowqa, mescoreCuts = NULL, tfscoreCuts = NULL, tlscoreCuts = NULL,
                              timeCuts = NULL, idCuts = NULL){
    iterms <- c("cell_ID", "time", "score", "meScore", "tfScore", "tlScore")
    
    scoreS <- flowQA_scoreSummary(flowqa$fsScore)
    fset <- flowqa$flowSet
    for(i in sampleNames(fset)){
        scorei <- scoreS[[i]]
        thresi <- thres[[i]]
        fcs <- fset[[i]]
    }        
}


## filter the data from QA
## -----------------------------------------------------------------------------------------------------
qcfilter <- function(x) {
    UseMethod("qcfilter", x)
}

qcfilter.cellnum <- function(x) {     
    outsample <- x$thresCount[1] < x$cellNumbers & x$cellNumbers < x$thresCount[2] 
    cellCut_fset <- x$flowSet[outsample]       
}

qcfilter.marginevents <- function(x) {
    tol <- x$toleranceValue
    mffset <- fsApply(x$flowSet, function(y){
        params <- parameters(y)
        keyval <- keyword(y)
        sub_exprs <- exprs(y)
        ranges <- range(y)
        for(m in seq_len(length(x$channels))){
            mexp <- sub_exprs[,m]
            summary(mexp)
            mrng <- ranges[,m]
            mrng
            okCells <- switch(x$side,
                              both= (mexp > (mrng[1]-tol)) & (mexp<(mrng[2]+tol)),
                              upper = mexp < (mrng[2]+tol), 
                              lower = mexp > (mrng[1]-tol))
            table(okCells)
            sub_exprs <- sub_exprs[okCells, ]
        }
        flowFrame(exprs = sub_exprs, parameters = params, description=keyval)
    }, simplify=TRUE)
}

qcfilter.qcScore <- function(x) {
    
}



# write.flowSet(fset, outDir)

## -----------------------------------------------------------------------------------------------------
## Guess which channel captures time in a exprs, flowFrame or flowset
findTimeChannel <- function(xx, strict=FALSE) {
    time <- grep("^Time$", colnames(xx), value=TRUE,
                 ignore.case=TRUE)[1]
    if(is.na(time)){
        if(is(xx, "flowSet")||is(xx, "ncdfFlowList"))
            xx <- exprs(xx[[1]])
        else if (is(xx, "flowFrame"))
            xx <- exprs(xx)
        cont <- apply(xx, 2, function(y) all(sign(diff(y)) >= 0))
        time <- names(which(cont))
    }
    if(!length(time) && strict)
        stop("Unable to identify time domain recording for this data.\n",
             "Please define manually.", call.=FALSE)
    if(length(time)>1)
        time <- character(0)
    return(time)
}

## Truncate the names to make them fit on the plot
truncNames <- function(names) {
    nc <- nchar(names)
    names[nc>11] <- paste(substr(names[nc>11], 1, 8), "...", sep="")
    if(any(duplicated(names))){
        ns <- split(names, names)
        is <- split(seq_along(names), names)
        ns<- lapply(ns, function(x)
            if(length(x)>1) paste(x, seq_along(x), sep="_") else x)
        names <- unlist(ns)[unlist(is)]
    }
    return(names)
}

## Lighten up or darken colors
desat <- function(col, by=50) {
    rgbcol <- col2rgb(col)
    up <- max(rgbcol-255)
    down <- min(rgbcol)
    if(by>0)
        rgbcol <- rgbcol + min(by, up)
    if(by<0)
        rgbcol <- rgbcol - min(abs(by), down)
    return(rgb(rgbcol[1,], rgbcol[2,], rgbcol[3,], maxColorValue=255))
}

SPADE.strip.sep <- function(name) {
  ifelse(substr(name,nchar(name),nchar(name))==.Platform$file,substr(name,1,nchar(name)-1),name)
}

SPADE.normalize.out_dir <- function(out_dir) {
  out_dir <- SPADE.strip.sep(out_dir)
  out_dir_info <- file.info(out_dir)
  if (is.na(out_dir_info$isdir)) {
    dir.create(out_dir)
  }
  if (!file.info(out_dir)$isdir) {
    stop(paste("out_dir:",out_dir,"is not a directory"))
  }
  out_dir <- paste(SPADE.strip.sep(out_dir),.Platform$file,sep="")
}

spade_plot_test<- function (graph, files, file_pattern = "*anno.Rsave", out_dir = ".", 
                       layout = SPADE.layout.arch, attr_pattern = "percent|medians|fold|cvs", 
                       scale = NULL, pctile_color = c(0.02, 0.98), normalize = "global", 
                       size_scale_factor = 1, edge.color = "grey", bare = FALSE, 
                       palette = "bluered") 
{
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.null(scale) && (!is.vector(scale) || length(scale) != 
                          2)) {
    stop("scale must be a two element vector")
  }
  if (!is.vector(pctile_color) || length(pctile_color) != 2) {
    stop("pctile_color must be a two element vector with values in [0,1]")
  }
  if (length(files) == 1 && file.info(SPADE.strip.sep(files))$isdir) {
    files <- dir(SPADE.strip.sep(files), full.names = TRUE, 
                 pattern = glob2rx(file_pattern))
  }
  out_dir <- SPADE.normalize.out_dir(out_dir)
  load_attr <- function(save_file) {
    anno <- NULL
    l <- load(save_file)
    stopifnot(l == "anno")
    return(anno)
  }
  boundaries <- NULL
  if (normalize == "global") {
    boundaries <- c()
    all_attrs <- c()
    for (f in files) {
      attrs <- load_attr(f)
      for (i in grep(attr_pattern, colnames(attrs))) {
        n <- colnames(attrs)[i]
        all_attrs[[n]] <- c(all_attrs[[n]], attrs[, i])
      }
    }
    for (i in seq_along(all_attrs)) {
      boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]], 
                                                    probs = pctile_color, na.rm = TRUE)
    }
  }
  if (is.function(layout)) 
    graph_l <- layout(graph)
  else graph_l <- layout
  if (palette == "jet") 
    palette <- colorRampPalette(c("#00007F", "blue", 
                                  "#007FFF", "cyan", "#7FFF7F", "yellow", 
                                  "#FF7F00", "red", "#7F0000"))
  else if (palette == "bluered") 
    palette <- colorRampPalette(c("blue", "#007FFF", 
                                  "cyan", "#7FFF7F", "yellow", "#FF7F00", 
                                  "red"))
  else stop("Please use a supported color palette.  Options are 'bluered' or 'jet'")
  colorscale <- palette(100)
  for (f in files) {
    attrs <- load_attr(f)
    vsize <- attrs$percenttotal
    vsize <- vsize/(max(vsize, na.rm = TRUE)^(1/size_scale_factor)) * 
      3 + 2
    vsize[is.na(vsize) | (attrs$count == 0)] <- 1
    for (i in 16){ # as.numeric(channel_index)) { #grep(attr_pattern, colnames(attrs))) {
      print(i)
      name <- colnames(attrs)[i]
      attr <- attrs[, i]
      if (!is.null(scale)) 
        boundary <- scale
      else if (normalize == "global") {
        boundary <- boundaries[[name]]
      }
      else boundary <- quantile(attr, probs = pctile_color, 
                                na.rm = TRUE)
      if (length(grep("^medians|percent|cvs", name))) 
        boundary <- c(min(boundary), max(boundary))
      else boundary <- c(-max(abs(boundary)), max(abs(boundary)))
      boundary <- round(boundary, 2)
      if (boundary[1] == boundary[2]) {
        boundary <- c(boundary[1] - 1, boundary[2] + 
                        1)
      }
      grad <- seq(boundary[1], boundary[2], length.out = length(colorscale))
      color <- colorscale[findInterval(attr, grad, all.inside = TRUE)]
      color[is.na(attr) | (attrs$count == 0)] <- "grey"
      if (grepl("^percenttotalratiolog$", name)) {
        color[is.na(attr) & attrs$count > 0] <- tail(colorscale, 
                                                     1)
      }
      fill_color <- color
      is.na(fill_color) <- is.na(attr)
      frame_color <- color
      #           pdf(paste(out_dir, basename(f), ".", name,  ".pdf", sep = ""))
      graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[, 
                                                        2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
      plot(graph, layout = graph_l, vertex.shape = "circle", 
           vertex.color = fill_color, vertex.frame.color = frame_color, 
           edge.color = edge.color, vertex.size = vsize, 
           vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1, 
           asp = graph_aspect)
      if (!bare) {
        if (length(grep("^medians", name))) 
          name <- sub("medians", "Median of ", 
                      name)
        else if (length(grep("^fold", name))) 
          name <- sub("fold", "Arcsinh diff. of ", 
                      name)
        else if (grepl("^percenttotal$", name)) 
          name <- sub("percent", "Percent freq. of ", 
                      name)
        else if (grepl("^percenttotalratiolog$", 
                       name)) 
          name <- "Log10 of Ratio of Percent Total of Cells in Each Cluster"
        else if (grepl("^cvs", name)) 
          name <- sub("cvs", "Coeff. of Variation of ", 
                      name)
        if (grepl("_clust$", name)) 
          name <- sub("_clust", "\n(Used for tree-building)", 
                      name)
        title(main = paste(strsplit(basename(f), ".fcs")[[1]][1], 
                           sub = name, sep = "\n"))
        # subplot(image(grad, c(1), matrix(1:length(colorscale), 
        #                                  ncol = 1), col = colorscale, xlab = ifelse(is.null(scale), 
        #                                                                             paste("Range:", pctile_color[1], "to", 
        #                                                                                   pctile_color[2], "pctile"), ""), 
        #               ylab = "", yaxt = "n", xaxp = c(boundary, 
        #                                               1)), x = "right,bottom", size = c(1, 
        #                                                                                 0.2))
      }
      #            dev.off()
    }
  }
}


