rm(list=ls())

args = (commandArgs(TRUE))

prefix = args[1]

source("common.r")

panel = read.table("/well/mcvean/joezhu/pf3k/pf3k_5_1_final/labStrains.eg.panel.txt", header=T)
#panel = read.table("labStrains.eg.panel.txt", header=T)

endAt = cumsum(table(panel[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

Ref1 = panel[,5] # HB3
Ref1Name = "HB3"
Ref2 = panel[,6] # 7G8
Ref2Name = "7G8"



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
    tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])

    hap = as.matrix(read.table(paste(prefix,".hap",sep=""), header=T)[,c(-1,-2)])

    colIndex = which(prop>0.01)

    prop.corrected = prop[colIndex]
    hap.corrected = hap[,colIndex,drop=FALSE]

    printed.prop = prop.corrected[getIndex2(hap.corrected, Ref1, Ref2)]

    switchError = c(0, 0)
    mutError = c(0, 0)
    for ( chrom in 1:length(beginAt)){
        tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
        tmpProp = prop.corrected
        tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
        tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]

        if ( length(prop.corrected) == 2 ){
            rearranged.Index = getIndex2(tmpHap, tmpRef1, tmpRef2)
            tmpHap = tmpHap[,rearranged.Index,drop=FALSE]
            tmpProp = prop.corrected[rearranged.Index]
        }

        hapAndError = fun.computeErrors2( tmpHap, tmpRef1, tmpRef2)

        tmpTitle = paste(prefix, rownames(table(panel[,1]))[chrom], sum(hapAndError$switchError), "switch errors", sum(hapAndError$mutError), "miss copy errors")

        fun.plotHapWithProp (hapAndError$hap, tmpProp,
             tmpTitle,
             max(chromLength))
	cat("switchError ", switchError, "\n")
	cat("hapAndError$switchError ", hapAndError$switchError,"\n")
        switchError = switchError + hapAndError$switchError
	cat("mutError ", mutError,"\n")
	cat("hapAndError$mutError ", hapAndError$mutError,"\n")
        mutError = mutError + hapAndError$mutError
    }
    if ( length(prop.corrected) == 2 ){
        write.table(cbind(printed.prop, switchError, mutError), file = paste(prefix,".errorCount", sep=""), quote = F, row.names=F, col.names=F)
    }
    dev.off()

#}
warnings()
