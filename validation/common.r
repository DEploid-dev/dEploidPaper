fun.divide.to.seg <- function(hapLength, numSeg = 501){
    return(floor(seq(1, hapLength, length.out =numSeg)))
}

getIndex3 <- function (mixedSample, ref1, ref2, ref3, prop){
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


getIndex2 <- function (mixedSample, ref1, ref2){
    sum1 = sum(mixedSample[,1] != ref1) + sum(mixedSample[,2] != ref2)
    sum2 = sum(mixedSample[,2] != ref1) + sum(mixedSample[,1] != ref2)
    if ( sum1 < sum2 ){
        return (c(1,2))
    } else {
        return (c(2,1))
    }
}


exportSites <- function ( CHROM, POS, tmpIndex, logFileName){
    for ( index in tmpIndex ){
        cat ( as.character(CHROM[index]), ",", POS[index], "\n", file = logFileName, append = T)
    }
}


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
#            } else if ( sum(hap[tmpIndex,1] != ref3[tmpIndex]) <= sum(hap[tmpIndex,1] != ref2[tmpIndex]) ) {
#                hap[tmpIndex,1][hap[tmpIndex,1] != ref3[tmpIndex]] = 2
##                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 8
#                countSwitch[1,i] = 3
#            }

            # k = 2
            if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) ){
                hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                countSwitch[2,i] = 2
            } else {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
#                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 8
                countSwitch[2,i] = 1
            }
#            } else if ( sum(hap[tmpIndex,2] != ref1[tmpIndex]) <= sum(hap[tmpIndex,2] != ref3[tmpIndex]) ) {
#                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
#                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
#                countSwitch[2,i] = 1
#            }
        }
        strain = cbind(c(NA,NA), countSwitch)
        strainNext = cbind(countSwitch, c(NA,NA))
        switchError = rowSums(strain != strainNext, na.rm=T)
        mutError = rowSum(hap==2)
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
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="Hap proportion", main=fig.title, xlab = "SNP", cex.main = 4)

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

fun.computeErrors3 <- function(hap, ref1, ref2, ref3){
    switchError = c(0, 0, 0)
    mutError = c(0, 0, 0)

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

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
        mutError = rowSum(hap==2)
    }

    return ( list ( hap = hap,
                    switchError = switchError,
                    mutError = mutError ))

}
