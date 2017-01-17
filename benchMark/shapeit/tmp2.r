rm(list=ls())

args = (commandArgs(TRUE))
#R --slave "--args sampleName" < tmp2.r

prefix = args[1]

source("../../validation/common.r")


fun.computeErrors2 <- function(hap, ref1, ref2, prop){
#    cat("prop parsed in : ", prop, "\n")
    switchError = c(0, 0)
    mutError = c(0, 0)

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]
threshold = 0.16
    if ( nhap == 2 ){
#        if (min(prop)<threshold){
#            index.of.seg = fun.divide.to.seg(haplength, by = 100)
#        } else {
#            index.of.seg = fun.divide.to.seg(haplength, by = 100)
#        }
index.of.seg = fun.divide.to.seg(haplength)
        countSwitch = matrix(NA, ncol = length(index.of.seg)-1, nrow = 2)

        for ( i in 1:(length(index.of.seg)-1) ){
            tmpIndex = c(index.of.seg[i]:index.of.seg[i+1])

            if (min(prop)<threshold){
                # k = 1
                if ( sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref2[tmpIndex]) & prop[1] >= (1-threshold)){
#                    print("here1")
                    hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
                    hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
                    countSwitch[1,i] = 1
                } else if ( sum(hap[tmpIndex,1] != ref1[tmpIndex]) > sum(hap[tmpIndex,1] != ref2[tmpIndex]) & prop[1] < threshold){
#                    print("here2")
                    hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
                    hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
                    countSwitch[1,i] = 1
                } else {
#                    cat(prop[1]," ")
#                    cat ("sum(hap[tmpIndex,1] != ref1[tmpIndex]) = ", sum(hap[tmpIndex,1] != ref1[tmpIndex])," ")
#                    cat ("sum(hap[tmpIndex,1] != ref2[tmpIndex]) = ", sum(hap[tmpIndex,1] != ref2[tmpIndex])," ")
#                    print("here3")
                    hap[tmpIndex,1][hap[tmpIndex,1] != ref2[tmpIndex]] = 2
                    hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 7
                    countSwitch[1,i] = 2
                }

                # k = 2
                if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) & prop[2] >= (1-threshold) ){
#                    print("here1")
                    hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                    hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                    countSwitch[2,i] = 2
                } else if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) > sum(hap[tmpIndex,2] != ref1[tmpIndex]) & prop[2] < threshold){
#                    print("here2")
                    hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                    hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                    countSwitch[2,i] = 2
                } else {
#                    print("here3")
                    hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
                    hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
    #                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 8
                    countSwitch[2,i] = 1
                }
            } else {
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
    #                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 8
                    countSwitch[2,i] = 1
                }
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



panel = read.table("labStrains.eg.panel.txt", header=T)
#panel = read.table("labStrains.eg.panel.txt", header=T)

endAt = cumsum(table(panel[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

Ref1 = panel[,5] # HB3
Ref1Name = "HB3"
Ref2 = panel[,6] # 7G8
Ref2Name = "7G8"

print(paste(prefix, ".prop", sep=""))
#sampleName = "PG0396-C"
#for ( seed in 1:15 ){

#    prefix = paste("PG0390-C/PG0390-C_seed", seed, "k2", sep="")
    if ( length(grep ( paste(c("PG0389-C", "PG0390-C", "PG0391-C", "PG0392-C", "PG0393-C", "PG0394-C") , collapse="|"), prefix)) > 0){
        Ref1 = panel[,3] # 3d7
        Ref1Name = "3d7"
        Ref2 = panel[,4] # Dd2
        Ref2Name = "Dd2"
    }

    png(paste(prefix, "compareHap.png", sep=""), width = 1920, height = 1080)
    ncol = ceiling(length(endAt)/2)
    par(mfrow = c(ncol,length(endAt)/ncol))
#    tmpProp = c(.5,.5)
    tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])
cat("read in prop: ", prop, "\n")
#    hap = as.matrix(read.table(paste(prefix,".gt.txt",sep=""), header=T)[,c(-1,-2)])
    hap = as.matrix(read.table(paste(prefix,".haps",sep=""), header=F)[,c(-1,-2)])

#    colIndex = which(prop>=0.00)
    colIndex= c(1,2)
    prop.corrected = prop[colIndex]
    hap.corrected = hap[,colIndex,drop=FALSE]

    correctedIndex = c(1,2)
#    correctedIndex = getIndex2(hap.corrected, Ref1, Ref2)
    printed.prop = prop.corrected[correctedIndex]
    switchError = c(0, 0)
    mutError = c(0, 0)
    for ( chrom in 1:length(beginAt)){
        tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
#        tmpProp = prop.corrected
        tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
        tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]

        if ( length(prop.corrected) == 2 ){
            rearranged.Index = getIndex2(tmpHap, tmpRef1, tmpRef2)
            tmpHap = tmpHap[,rearranged.Index,drop=FALSE]
#            tmpProp = prop.corrected[rearranged.Index]
        }
        hapAndError = fun.computeErrors2( tmpHap, tmpRef1, tmpRef2, prop)

        tmpTitle = paste(rownames(table(panel[,1]))[chrom], sum(hapAndError$switchError), "switch errors", sum(hapAndError$mutError), "miss copy errors")

        fun.plotHapWithProp (hapAndError$hap, tmpProp,
             tmpTitle,
             max(chromLength))
#	cat("switchError ", switchError, "\n")
#	cat("hapAndError$switchError ", hapAndError$switchError,"\n")
        switchError = switchError + hapAndError$switchError
#	cat("mutError ", mutError,"\n")
#	cat("hapAndError$mutError ", hapAndError$mutError,"\n")
        mutError = mutError + hapAndError$mutError
    }
    if ( length(prop.corrected) == 2 ){
        tmp = cbind(printed.prop, switchError, mutError)
        tmpIndex = c(1,2)
        if ( length(grep ( paste(c("PG0393-C", "PG0394-C", "PG0399-C", "PG0401-C", "PG0402-C", "PG0404-C", "PG0405-C") , collapse="|"), prefix)) > 0){
            print("here flipping")
            tmpIndex = c(2,1)
        }
        write.table(tmp[tmpIndex,], file = paste(prefix,".errorCount", sep=""), quote = F, row.names=F, col.names=F)
    }
    dev.off()

#}
warnings()
