rm(list= ls())
source("common.r")

getIndex <- function (mixedSample, ref1, ref2, ref3, prop){
#1	1	2	2	3	3
#2	3	1	3	2	1
#3	2	3	1	1	2
#    prop = c(.25,.25,.5)
    prop = c(1,1,1)
    mySums = c ( sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
                 sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] ,
                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] )
#print(mySums)
#print("")
    case = which.min(mySums)
    if ( case == 1 ){
        return (c(1,2,3))
    } else if ( case == 2 ){
        return (c(1,3,2))
    } else if ( case == 3 ){
        return (c(2,1,3))
    } else if ( case == 4 ){
        return (c(2,3,1))
    } else if ( case == 5 ){
        return (c(3,2,1))
    } else if ( case == 6 ){
        return (c(3,1,2))
    }
}


fun.plotHapWithProp <- function( hap, prop, ref1, ref2, ref3, startIndex, endIndex, fig.title, max.at ){

    switchError = 0
    mutError = 0

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]
    xrange = c(0, max.at)
    yrange = c(0, 1)
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="Hap proportion", main=fig.title, xlab = "SNP")

    xleft = 0:(haplength-1)
    xright = xleft+1
    ycum = c(0, cumsum(prop))

    if ( nhap == 3 ){
#        hap[hap[,1] != ref1, 1 ] = 2
#        hap[hap[,2] != ref2, 2 ] = 2
#        hap[hap[,3] != ref3, 3 ] = 2

        index.of.seg = fun.divide.to.seg(haplength)
        countSwitch = matrix(NA, ncol = length(index.of.seg)-1, nrow = 3)

        for ( i in 1:(length(index.of.seg)-1) ){
            tmpIndex = c(index.of.seg[i]:index.of.seg[i+1])
            # k = 1
            if ( sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref2[tmpIndex]) & sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref3[tmpIndex]) ){
                hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
                countSwitch[1,i] = 1
            } else if ( sum(hap[tmpIndex,1] != ref2[tmpIndex]) < sum(hap[tmpIndex,1] != ref3[tmpIndex]) ) {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref2[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 7
                countSwitch[1,i] = 2
            } else if ( sum(hap[tmpIndex,1] != ref3[tmpIndex]) <= sum(hap[tmpIndex,1] != ref2[tmpIndex]) ) {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref3[tmpIndex]] = 2
#                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 8
                countSwitch[1,i] = 3
            }

            # k = 2
            if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) & sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref3[tmpIndex]) ){
                hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                countSwitch[2,i] = 2
            } else if ( sum(hap[tmpIndex,2] != ref3[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) ) {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref3[tmpIndex]] = 2
#                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 8
                countSwitch[2,i] = 3
            } else if ( sum(hap[tmpIndex,2] != ref1[tmpIndex]) <= sum(hap[tmpIndex,2] != ref3[tmpIndex]) ) {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
                countSwitch[2,i] = 1
            }

            # k = 3
            if ( sum(hap[tmpIndex,3] != ref3[tmpIndex]) < sum(hap[tmpIndex,3] != ref2[tmpIndex]) & sum(hap[tmpIndex,3] != ref3[tmpIndex]) < sum(hap[tmpIndex,3] != ref1[tmpIndex]) ){
                hap[tmpIndex,3][hap[tmpIndex,3] != ref3[tmpIndex]] = 2
#                hap[tmpIndex,3][hap[tmpIndex,3] == 0] = 8
                countSwitch[3,i] = 3
            } else if ( sum(hap[tmpIndex,3] != ref2[tmpIndex]) < sum(hap[tmpIndex,3] != ref1[tmpIndex]) ) {
                hap[tmpIndex,3][hap[tmpIndex,3] != ref2[tmpIndex]] = 2
                hap[tmpIndex,3][hap[tmpIndex,3] == 0] = 7
                countSwitch[3,i] = 2
            } else if ( sum(hap[tmpIndex,3] != ref1[tmpIndex]) <= sum(hap[tmpIndex,3] != ref2[tmpIndex]) ) {
                hap[tmpIndex,3][hap[tmpIndex,3] != ref1[tmpIndex]] = 2
                hap[tmpIndex,3][hap[tmpIndex,3] == 0] = 5
                countSwitch[3,i] = 1
            }
        }
        strain = cbind(c(NA,NA,NA), countSwitch)
        strainNext = cbind(countSwitch, c(NA,NA,NA))
        switchError = rowSums(strain != strainNext, na.rm=T)
        mutError = sum(hap==2)
    }


    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }

    return ( list ( switchError = switchError,
                    mutError = mutError ))
}


exportSites <- function ( CHROM, POS, tmpIndex, logFileName){
    for ( index in tmpIndex ){
        cat ( as.character(CHROM[index]), ",", POS[index], "\n", file = logFileName, append = T)
    }
}


dataDir = "./"
panel = read.table(paste(dataDir,"labStrains.14.panel.txt",sep=""),header=T)

endAt = cumsum(table(panel[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

Ref1 = panel[,4] # Dd2
Ref1Name = "Dd2"
Ref2 = panel[,5] # HB3
Ref2Name = "HB3"
Ref3 = panel[,6] # 7G8
Ref3Name = "7G8"

sampleName = "PG0396-C"
cases = paste(sampleName, c(".14.noPanel",
                            ".14.asia",
                            ".14.asiaPlus",
#                            ".14.asiaPlus2",
                            ".14.labPanel"), sep="")

png(paste("differentPanelForSample.", sampleName, ".png", sep=""), width = 1920, height = 1080)
par ( mfrow = c(length(cases),1))

for ( prefix in cases ){
print(prefix)
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
    tmpRef3 = Ref3[beginAt[chrom]:endAt[chrom]]

    #    cat("chrom", chrom, which(tmpHap[,1] != tmpRef1),", sum is ",sum(tmpHap[,1] != tmpRef1), "\n")
    tmpTitle = ""
    if ( length(prop.corrected) == 3 ){
        rearranged.Index = getIndex(hap.corrected, tmpRef1, tmpRef2, tmpRef3, tmpProp)
        tmpHap = tmpHap[,rearranged.Index]
        tmpProp = tmpProp = prop.corrected[rearranged.Index]
        tmpTitle = paste(prefix, "Chrom", chrom, ",", sum(tmpHap[,1] != tmpRef1), "sites differ from", Ref1Name, ",", sum(tmpHap[,2] != tmpRef2),"site differ from", Ref2Name, ",", sum(tmpHap[,3] != tmpRef3),"site differ from", Ref3Name)
    }

    fun.plotHapWithProp (tmpHap, tmpProp,
         tmpRef1, tmpRef2, tmpRef3,
         beginAt[chrom], endAt[chrom],
         tmpTitle,
         max(chromLength))
}

dev.off()

