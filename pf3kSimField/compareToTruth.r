rm(list=ls())

fun.divide.to.seg <- function(hapLength, by = 3001){
#    return(floor(seq(1, hapLength, length.out =numSeg)))
    myseq = seq(1, hapLength, by = by)
    if ((hapLength-1)%%by != 0){
        myseq = c(myseq, hapLength)
    }
    return(myseq)
}

#getIndex3 <- function (mixedSample, ref1, ref2, ref3, prop){
##1	1	2	2	3	3
##2	3	1	3	2	1
##3	2	3	1	1	2
##    prop = c(.25,.25,.5)
#    prop = c(1,1,1)
#    mySums = c ( sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
#                 sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] ,
#                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
#                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
#                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
#                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] )
##print(mySums)
##print("")
#    case = which.min(mySums)
#    if ( case == 1 ){
#        return (c(1,2,3))
#    } else if ( case == 2 ){
#        return (c(1,3,2))
#    } else if ( case == 3 ){
#        return (c(2,1,3))
#    } else if ( case == 4 ){
#        return (c(2,3,1))
#    } else if ( case == 5 ){
#        return (c(3,2,1))
#    } else if ( case == 6 ){
#        return (c(3,1,2))
#    }
#}


getIndex2 <- function (mixedSample, ref1, ref2){
    sum1 = sum(mixedSample[,1] != ref1) + sum(mixedSample[,2] != ref2)
    sum2 = sum(mixedSample[,2] != ref1) + sum(mixedSample[,1] != ref2)
    if ( sum1 < sum2 ){
        return (c(1,2))
    } else {
        return (c(2,1))
    }
}


#exportSites <- function ( CHROM, POS, tmpIndex, logFileName){
#    for ( index in tmpIndex ){
#        cat ( as.character(CHROM[index]), ",", POS[index], "\n", file = logFileName, append = T)
#    }
#}


fun.computeErrors2 <- function(hap, ref1, ref2){
    switchError = c(0, 0)
    mutError = c(0, 0)

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    if ( nhap == 2 ){

        index.of.seg = fun.divide.to.seg(haplength)
        countSwitch = matrix(NA, ncol = length(index.of.seg)-1, nrow = 2)

        for ( i in 1:(length(index.of.seg)-1) ){
            tmpIndex = c(index.of.seg[i]:index.of.seg[i+1])
            # k = 1
            if ( sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref2[tmpIndex]) ){
                hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
                countSwitch[1,i] = 1
            } else {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref2[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 7
                countSwitch[1,i] = 2
            }

            # k = 2
            if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) ){
                hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                countSwitch[2,i] = 2
            } else {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
                countSwitch[2,i] = 1
            }
        }
        strain = cbind(c(NA,NA), countSwitch)
        strainNext = cbind(countSwitch, c(NA,NA))
        switchError = rowSums(strain != strainNext, na.rm=T)
        mutError = colSums(hap==2)
    }

    return ( list ( hap = hap,
                    switchError = switchError,
                    mutError = mutError ))

}


fun.plotHapWithProp <- function( hap, prop, fig.title, max.at ){

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, max.at)
    yrange = c(0, 1)
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 2, cex.axis=2)

    xleft = 0:(haplength-1)
    xright = xleft+1
    ycum = c(0, cumsum(prop))

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }
}


asia1 = read.table("asia1.14.panel.txt", header=T, check.names=F)
endAt = cumsum(table(asia1[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

Ref1 = asia1[["PH0064-C"]]
Ref1Name = "PH0064-C"
Ref2 = asia1[["PH0193-C"]]
Ref2Name = "PH0193-C"
RefNames=c(Ref1Name, Ref2Name)

#Ref2 = asia1[["PH0064-C"]]
#Ref2Name = "PH0064-C"
#Ref1 = asia1[["PH0193-C"]]
#Ref1Name = "PH0193-C"


#if ( currentDir %in% c("PG0389-C", "PG0390-C", "PG0391-C", "PG0392-C", "PG0393-C", "PG0394-C") ){
#}

#cases = c("PG0412-C.14.noPanel",
#          "PG0412-C.14.asia",
#          "PG0412-C.14.asiaPlus",
#          "PG0412-C.14.labPanel")

#suffix = ".noError"
#suffix = ".withError"
suffix = ".noErrormiss0001"
#suffix = ".withErrormiss0001"
#sampleName = "25v75"
sampleName = "75v25"
#sampleName = "55v45"
#cases = paste(sampleName, "panel", 1, suffix, sep="")
cases = c(paste(sampleName, "panel", c(1,2,3), suffix, sep=""),
          paste(sampleName, "noPanel", suffix, sep = ""))

png(paste("differentPanelForSample.", sampleName, suffix, ".png", sep=""), width = 1920, height = 1080)
par ( mfrow = c(length(cases),1))

for ( prefix in cases ){

    tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])

    hap = as.matrix(read.table(paste(prefix,".hap",sep=""), header=T)[,c(-1,-2)])

    colIndex = which(prop>0.01)

    prop.corrected = prop[colIndex]
    hap.corrected = hap[,colIndex,drop=FALSE]

    switchError = 0
    mutError = 0
    for ( chrom in 1:length(beginAt)){
        tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
        tmpProp = prop.corrected
        tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
        tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]

        if ( length(prop.corrected) == 2 ){
            rearranged.Index = getIndex2(hap.corrected, Ref1, Ref2)
            tmpHap = tmpHap[,rearranged.Index,drop=FALSE]
            tmpProp = tmpProp = prop.corrected[rearranged.Index]
        }

        hapAndError = fun.computeErrors2( tmpHap, tmpRef1, tmpRef2)

        tmpTitle = paste(prefix, rownames(table(asia1[,1]))[chrom], sum(hapAndError$switchError), "switch errors", RefNames, hapAndError$mutError, "miss copy errors")

        fun.plotHapWithProp (hapAndError$hap, tmpProp,
             tmpTitle,
             max(chromLength))
        switchError = switchError + hapAndError$switchError
        mutError = mutError + hapAndError$mutError
    }
    write.table(c(switchError, mutError), file = paste(prefix,".errorCount", sep=""), quote = F, row.names=F, col.names=F)
}

dev.off()
