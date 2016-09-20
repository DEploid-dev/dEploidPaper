rm(list=ls())
source("common.r")

getIndex <- function (mixedSample, ref1, ref2){
    sum1 = sum(mixedSample[,1] != ref1) + sum(mixedSample[,2] != ref2)
    sum2 = sum(mixedSample[,2] != ref1) + sum(mixedSample[,1] != ref2)
#    print(sum1)
#    print(sum2)
    if ( sum1 < sum2 ){
        return (c(1,2))
    } else {
        return (c(2,1))
    }
}


fun.divide.to.seg <- function(hapLength, by = 50){
    myseq = seq(1, hapLength, by = by)
    if ((hapLength-1)%%by != 0){
        myseq = c(myseq, hapLength)
    }
    return(myseq)
}


#dataDir = "./"
panel = read.table(paste("labStrains.14.panel.txt",sep=""),header=T)
endAt = cumsum(table(panel[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

Ref1 = panel[,5] # HB3
Ref1Name = "HB3"
Ref2 = panel[,6] # 7G8
Ref2Name = "7G8"

sample = "PG0402-C"
suffix = "lab"
subSamples = c(20, 50, 80, 100)
#subSamples = c(50)
for ( subSample in subSamples ) {
    vcfPrefix = paste(sample, ".subSample", subSample, sep="")
    prefix = vcfPrefix
    vcfName = paste(vcfPrefix,".vcf.gz", sep="")
    coverage = fun.extract.vcf(vcfName)

    totalCov = coverage$refCount + coverage$altCount
    spacing = ceiling(max(totalCov)/15)
    mybins = seq(0, ceiling(max(totalCov)/spacing)*spacing, by = spacing)
    nBins = length(mybins)

    #mybins = quantile(c(totalCov, max(totalCov)+1), prob=(1:nBins)/nBins)
    binIndex = findInterval(totalCov, mybins)
    eventsType = c("0/0", "0/1", "1/0", "1/1")


    eventCountMat = c()#data.frame ()
    eventSumMat = c()
    eventArray = c()
    #prefix = paste(vcfPrefix,".asia.out", sep="")
    for ( seed in 1:15 ){
        outprefix = paste("repeats/", sample, ".seed", seed, ".subSample", subSample,".", suffix, ".out", sep="")

    #    tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
    #    prop = as.numeric(tmpProp[dim(tmpProp)[1],])
        hap.corrected = as.matrix(read.table(paste(outprefix,".hap",sep=""), header=T)[,c(-1,-2)])

    #    colIndex = which(prop>0.01)
    #    prop.corrected = prop[colIndex]
    #    hap.corrected = hap[,colIndex,drop=FALSE]

        chrom = 1

        tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
#        tmpProp = prop.corrected
        tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
        tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]

#        rearranged.Index = getIndex(tmpHap, tmpRef1, tmpRef2)
#        tmpHap = tmpHap[,rearranged.Index,drop=FALSE]
#        tmpProp = tmpProp[rearranged.Index]


        haplength = dim(tmpHap)[1]
        index.of.seg = fun.divide.to.seg(haplength)

        truth = c()
        infered = c()
        for ( i in 1:(length(index.of.seg)-1) ){
            tmpIndex = c(index.of.seg[i]:index.of.seg[i+1])

            tmptmpHap = tmpHap[tmpIndex,]
            tmptmpRef1 = tmpRef1[tmpIndex]
            tmptmpRef2 = tmpRef2[tmpIndex]
            rearranged.Index = getIndex(tmptmpHap, tmptmpRef1, tmptmpRef2)
            tmptmpHap = tmptmpHap[,rearranged.Index,drop=FALSE]

            truth = c(truth, paste(tmptmpRef1, "/", tmptmpRef2, sep=""))
            infered = c(infered, paste(tmptmpHap[,1], "/", tmptmpHap[,2], sep=""))
        }

        for ( event in eventsType ){
            eventCount = rep(0, nBins)
            eventSum = rep(0, nBins)
            eventIndex = which(truth==event)
            tmpevent = (truth[eventIndex]!=infered[eventIndex])
#            cat("length(eventIndex) ",length(eventIndex), " length(tmpevent)", length(tmpevent), "\n")
            for ( i in 1:nBins ){
                tmpIndex = which(binIndex[eventIndex]==(i))
        #        cat(length(tmpIndex),"\n")
                eventSum[i] = sum(tmpevent[tmpIndex]*1)
                eventCount[i] = length(tmpIndex)
            }
#            cat(sum(eventCount),"\n")
#            eventSum[eventCount<5] = 0
            eventCountMat = rbind(eventCountMat, eventCount)
            eventSumMat = rbind(eventSumMat, eventSum)
            eventArray = c(eventArray, event)
    #        lines(mybins,eventSum/(eventCount+0.00000001), col=color)
        #    lines(mybins,eventSum/(sum(truth!=infered)), col=color)
        #lines(mybins,eventSum, col=color)
        #    print(eventSum)
        #    print(eventCount)
    #        color = color+1
    #        print(sum(eventSum))
        }
    }
#    colnames(eventCountMat) = as.character(1:15)
#    rownames(eventCountMat) = as.character(1:60)
#    eventCount.frame = data.frame(eventCountMat)
#    eventCount.frame[[event]] = eventArray

#    colnames(eventSumMat) = as.character(1:15)
#    rownames(eventSumMat) = as.character(1:60)
#    eventSum.frame = data.frame(eventSumMat)
#    eventSum.frame[[event]] = eventArray
#    mytitle = paste(prefix, Ref1Name, round(tmpProp[1], digits=3), "/", Ref2Name, round(tmpProp[2],digits=3))
    png(paste(prefix, ".errorVsCoverage.png",sep=""), width=800, height=800)
    #par(mfrow=c(1,2))

    layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = TRUE))
    plot(c(min(mybins),max(mybins)),c(0, 0.3), type="n", ylab="# of sites was wrongly inferred", xlab="Total coverage")
    color = 1
    for ( event in eventsType ){
        tmpCount = colSums(eventCountMat[eventArray == event,])
        tmpSum = colSums(eventSumMat[eventArray == event,])
#        tmpCount = colMeans(eventCountMat[eventArray == event,])
#        tmpSum = colMeans(eventSumMat[eventArray == event,])
        tmpSum[tmpCount<100] = 0
print(tmpSum)
print(tmpCount)
        lines(mybins,tmpSum/(tmpCount+0.00000001), col=color)
        color = color+1
    }

    legend("topright", legend=c("HB3/7G8",eventsType), col=c(1,1:4), lty=c(0,1,1,1,1))

    hist(totalCov, breaks=mybins, ylim = c(0, 400), col = rgb(1,0,0,0.5), xlab = "Coverage / Alternative allele count" )
    #obj = mpileAlt[,3][mpileAlt[,3]>0]
    hist(coverage$altCount, breaks=seq(0, ceiling(max(coverage$altCount)/spacing)*spacing, by = spacing), add = T, col = rgb(0,0,1,0.5))
    legend("topright", c("Total coverage", "Alt count"), fill=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

    dev.off()
}
