rm(list= ls())
# This script is used to compare Jason's sample with the panel
#    R --slave "--args panel" < ../compareToPanel.r
#    R --slave "--args panelFree" < ../compareToPanel.r

#args=(commandArgs(TRUE))

#addPrefix = args[1] # "_panelFreeExclude"

#if ( addPrefix == "NOTHING" ){
#    addPrefix = ""
#}

getIndex <- function (mixedSample, ref1, ref2){
    sum1 = sum(mixedSample[,1] != ref1) + sum(mixedSample[,2] != ref2)
    sum2 = sum(mixedSample[,2] != ref1) + sum(mixedSample[,1] != ref2)
    if ( sum1 < sum2 ){
        return (c(1,2))
    } else {
        return (c(2,1))
    }
}


fun.plotHapWithProp <- function( hap, prop, ref1, ref2, startIndex, endIndex, fig.title, max.at ){

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]
    xrange = c(0, max.at)
    yrange = c(0, 1)
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="Hap proportion", main=fig.title, xlab = "SNP")
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, xlab = "", ylab="")

    xleft = 0:(haplength-1)
    xright = xleft+1
    ycum = c(0, cumsum(prop))

    if ( nhap == 2 ){
        hap[hap[,1] != ref1, 1 ] = 2
        hap[hap[,2] != ref2, 2 ] = 2
    }

    if ( nhap == 1 && currentDir == "PG0398-C" ){
        hap[hap[,1] != ref1, 1 ] = 2
    }

    if ( nhap == 1 && currentDir == "PG0415-C" ){
        hap[hap[,1] != ref2, 1 ] = 2
    }


    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }

#    if ( nhap == 2 ){
#        flip.at = which((ref1 != hap[,1] & ref2 != hap[,2]) & ref1 != ref2 )

#        hap[flip.at,1] = 3
#        hap[flip.at,2] = 3
#    #    cat( "flip.at = ", flip.at , "\n")

#        for ( k in c(1:nhap) ){
#          ybottom = ycum[k]
#          ytop = ycum[k+1]
#          rect(xleft, ybottom, xright, ytop, col =  hap[,k], border = "transparent")
#        }
#    }
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

Ref1 = panel[,5] # HB3
Ref1Name = "HB3"
Ref2 = panel[,6] # 7G8
Ref2Name = "7G8"

#Ref1 = panel[,3] # 3d7
#Ref1Name = "3d7"
#Ref2 = panel[,4] # Dd2
#Ref2Name = "Dd2"

#if ( currentDir %in% c("PG0389-C", "PG0390-C", "PG0391-C", "PG0392-C", "PG0393-C", "PG0394-C") ){
#}

#cases = c("PG0412-C.14.noPanel",
#          "PG0412-C.14.asia",
#          "PG0412-C.14.asiaPlus",
#          "PG0412-C.14.labPanel")

sampleName = "PG0406-C"
cases = paste(sampleName, c(".14.noPanel",
                            ".14.asia",
                            ".14.asiaPlus",
                            ".14.asiaPlus2",
                            ".14.labPanel"), sep="")

png(paste("differentPanelForSample.", sampleName, ".png", sep=""), width = 1920, height = 1080)
par ( mfrow = c(length(cases),1))

for ( prefix in cases ){

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
    #    cat("chrom", chrom, which(tmpHap[,1] != tmpRef1),", sum is ",sum(tmpHap[,1] != tmpRef1), "\n")
    tmpTitle = ""
    if ( length(prop.corrected) == 2 ){
        rearranged.Index = getIndex(hap.corrected, Ref1, Ref2)
        tmpHap = tmpHap[,rearranged.Index,drop=FALSE]
        tmpProp = tmpProp = prop.corrected[rearranged.Index]
        tmpTitle = paste(prefix, "Chrom", chrom, ",", sum(tmpHap[,1] != tmpRef1), "sites differ from", Ref1Name, ",", sum(tmpHap[,2] != tmpRef2),"site differ from", Ref2Name)
    }

#    tmpTitle=""
    fun.plotHapWithProp (tmpHap, tmpProp,
         tmpRef1, tmpRef2,
         beginAt[chrom], endAt[chrom],
         tmpTitle,
         max(chromLength))
}

dev.off()
