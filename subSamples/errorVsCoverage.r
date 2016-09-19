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


dataDir = "./"
panel = read.table(paste(dataDir,"labStrains.14.panel.txt",sep=""),header=T)
endAt = cumsum(table(panel[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

Ref1 = panel[,5] # HB3
Ref1Name = "HB3"
Ref2 = panel[,6] # 7G8
Ref2Name = "7G8"

sample = "PG0417-C"
suffix = "asia"
subSamples = c(20, 50, 80, 100)
for ( subSample in subSamples ) {
vcfPrefix = paste(sample, ".subSample", subSample, sep="")
#prefix = paste(vcfPrefix,".asia.out", sep="")
prefix = paste(vcfPrefix,".", suffix, ".out", sep="")

tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
prop = as.numeric(tmpProp[dim(tmpProp)[1],])

hap = as.matrix(read.table(paste(prefix,".hap",sep=""), header=T)[,c(-1,-2)])

colIndex = which(prop>0.01)

prop.corrected = prop[colIndex]
hap.corrected = hap[,colIndex,drop=FALSE]

chrom = 1


tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
tmpProp = prop.corrected
tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]

rearranged.Index = getIndex(tmpHap, tmpRef1, tmpRef2)
tmpHap = tmpHap[,rearranged.Index,drop=FALSE]
tmpProp = tmpProp[rearranged.Index]

fun.divide.to.seg <- function(hapLength, numSeg = 501){
    return(floor(seq(1, hapLength, length.out =numSeg)))
}

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


vcfName = paste(vcfPrefix,".vcf", sep="")
coverage = fun.extract.vcf(vcfName)

totalCov = coverage$refCount + coverage$altCount
#nBins = 6

png(paste(prefix, ".errorVsCoverage.png",sep=""), width=800, height=800)
#par(mfrow=c(1,2))

layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = TRUE))


spacing = ceiling(max(totalCov)/15)

mybins = seq(0, ceiling(max(totalCov)/spacing)*spacing, by = spacing)
nBins = length(mybins)
#mybins = quantile(c(totalCov, max(totalCov)+1), prob=(1:nBins)/nBins)
binIndex = findInterval(totalCov, mybins)
eventsType = c("0/0", "0/1", "1/0", "1/1")
color = 1
mytitle = paste(prefix, Ref1Name, round(tmpProp[1], digits=3), "/", Ref2Name, round(tmpProp[2],digits=3))
plot(c(min(mybins),max(mybins)),c(0, 0.6), type="n", ylab="# of sites was wrongly inferred", xlab="Total coverage", main = mytitle)

for ( event in eventsType ){
    eventCount = rep(0, nBins)
    eventSum = rep(0, nBins)
    eventIndex = which(truth==event)
    tmpevent = (truth[eventIndex]!=infered[eventIndex])
    cat("length(eventIndex) ",length(eventIndex), " length(tmpevent)", length(tmpevent), "\n")
    for ( i in 1:nBins ){
        tmpIndex = which(binIndex[eventIndex]==(i))
#        cat(length(tmpIndex),"\n")
        eventSum[i] = sum(tmpevent[tmpIndex]*1)
        eventCount[i] = length(tmpIndex)
    }
    cat(sum(eventCount),"\n")
    eventSum[eventCount<5] = 0

    lines(mybins,eventSum/(eventCount+0.00000001), col=color)
#    lines(mybins,eventSum/(sum(truth!=infered)), col=color)
#lines(mybins,eventSum, col=color)
#    print(eventSum)
#    print(eventCount)
    color = color+1
    print(sum(eventSum))
}
legend("topright", legend=c("HB3/7G8",eventsType), col=c(1,1:4), lty=c(0,1,1,1,1))

print(sum(truth!=infered))
print("")
#nBins = 5

#gatkCol = rgb(1,0,0,0.5)
#mpileCol = rgb(0,0,1,0.5)
hist(totalCov, breaks=mybins, ylim = c(0, 400))#, col = gatkCol, xlab = "Alternative allele count" )
#obj = mpileAlt[,3][mpileAlt[,3]>0]
hist(coverage$altCount, breaks=seq(0, ceiling(max(coverage$altCount)/spacing)*spacing, by = spacing), add = T, col = rgb(0,0,1,0.5))
#legend("topright", c("GATK", "MpileUp"), fill=c(gatkCol, mpileCol))

dev.off()


}

#mybins = quantile(c(coverage$altCount, max(coverage$altCount)+1), prob=(1:nBins)/nBins)
#binIndex = findInterval(coverage$altCount, mybins)
#eventsType = c("0/0", "0/1", "1/0", "1/1")
#color = 1
#mytitle = paste(prefix, Ref1Name, round(tmpProp[1], digits=3), "/", Ref2Name, round(tmpProp[2],digits=3))
#plot(c(min(mybins),max(mybins)),c(0,0.6), type="n", ylab="Error rate", xlab="ALT count", main = mytitle)

#for ( event in eventsType ){
#    eventCount = rep(0, nBins)
#    eventSum = rep(0, nBins)
#    eventIndex = which(truth==event)
#    tmpevent = (truth[eventIndex]!=infered[eventIndex])
#    cat("length(eventIndex) ",length(eventIndex), " length(tmpevent)", length(tmpevent), "\n")
#    for ( i in 1:(nBins) ){
#        tmpIndex = which(binIndex[eventIndex]==(i-1))
##        cat(length(tmpIndex),"\n")
#        eventSum[i] = sum(tmpevent[tmpIndex]*1)
#        eventCount[i] = length(tmpIndex)
#    }
#        cat(sum(eventCount),"\n")
#    lines(mybins,eventSum/eventCount, col=color)
#    color = color+1
#    print(sum(eventSum))
#}
#legend("topright", legend=c("HB3/7G8",eventsType), col=c(1,1:4), lty=c(0,1,1,1,1))
#print(sum(truth!=infered))


#dev.off()
