rm(list=ls())
#png("switchVsMisCopyErrlog.png", width=800, height = 800)
#png("switchVsMisCopyErrlogHi.png", width=4000, height = 4000, res = 350)
tiff("switchVsMisCopyErrlogHi.tif", width=4000, height = 4000, res = 350)


simpleMix = paste("PG03", 89:94, "-C", sep="")
#mixedSampleNames = c(paste("PG03", 89:99, "-C", sep=""), paste("PG040", 0:9, "-C", sep=""), paste("PG04", 10:15, "-C", sep=""))

#plot(c(0,1), c(0, 2000), type="n")
mytable = data.frame ( prop = numeric(0), coverage = numeric(0), switchError = numeric(0), missCopy = numeric(0), sampleName = character(), strainName = character(), labStrain = character(), stringsAsFactors = F)
mytablei = 1


metaData = read.table("pf3k_release_5_metadata.txt", header = T, comment.char="", sep="\t", stringsAsFactors=F)

sample2 = read.table("../labSampleNames2Strains", header=F, stringsAsFactors=F)$V1
for (samplei in sample2){
    for ( seed in 1:15 ){
        prefix = paste(samplei, "_seed", seed, "k2", sep="")
        fileName = paste(prefix, ".errorCount", sep ="")
        if ( file.exists(fileName) ){
            tmpInfo = read.table(fileName, header=F)
            for ( i in 1:dim(tmpInfo)[1] ){
                mytable[mytablei,]$prop = tmpInfo$V1[i]
                mytable[mytablei,]$coverage = metaData$mean_coverage[metaData$sample == samplei]
                mytable[mytablei,]$switchError = tmpInfo$V2[i]
                mytable[mytablei,]$missCopy = tmpInfo$V3[i]
                mytable[mytablei,]$sampleName = samplei
                mytable[mytablei,]$strainName = paste(samplei,".", i, sep="")
                labStrain = ""
                if ( i == 1 ){
                    labStrain = "hb3"
                    if ( samplei %in% simpleMix) { labStrain = "3d7" }
                } else {
                    labStrain = "7g8"
                    if ( samplei %in% simpleMix) { labStrain = "dd2" }
                }
                mytable[mytablei,]$labStrain = labStrain
                mytablei = mytablei+1
            }
        }
    }
}

sample3 = read.table("../labSampleNames3Strains", header=F, stringsAsFactors=F)$V1
for (samplei in sample3){
    for ( seed in 1:15 ){
        prefix = paste(samplei, "_seed", seed, "k3", sep="")
        fileName = paste(prefix, ".errorCount", sep ="")
        if ( file.exists(fileName) ){
            tmpInfo = read.table(fileName, header=F)
            for ( i in 1:dim(tmpInfo)[1] ){
                mytable[mytablei,]$prop = tmpInfo$V1[i]
                mytable[mytablei,]$coverage = metaData$mean_coverage[metaData$sample == samplei]
                mytable[mytablei,]$switchError = tmpInfo$V2[i]
                mytable[mytablei,]$missCopy = tmpInfo$V3[i]
                mytable[mytablei,]$sampleName = samplei
                mytable[mytablei,]$strainName = paste(samplei,".", i, sep="")
                labStrain = ""
                if ( i == 1 ){
                    labStrain = "hb3"
                } else if ( i == 2) {
                    labStrain = "7g8"
                } else {
                    labStrain = "dd2"
                }
                mytable[mytablei,]$labStrain = labStrain
                mytablei = mytablei+1
            }
        }
    }

}

mytable[["color"]] = 2
mytable$color[mytable$sampleName %in% c("PG0395-C","PG0396-C","PG0397-C")] = 4
mytable[["marker"]] = factor(mytable$labStrain)

par(mfrow = c(2,1))

plot(mytable$prop, mytable$switchError, ylab= "Number of switches when copying", xlab="Strain proportion", col=mytable$color, log="y", xlim = c(0,1), type="n", cex.lab=1.5)
strains = unique(mytable$strainName)
for ( strain in strains ){
    tmpIndex = which(mytable$strainName == strain)
    tmpMarker = as.numeric(unique(mytable$marker[tmpIndex]))
    tmpColor = unique(mytable$color[tmpIndex])
    points( mean(mytable$prop[tmpIndex]), mean(mytable$switchError[tmpIndex]), pch = tmpMarker, col=tmpColor, cex=2, lwd=2)
}
legend("topright", legend = levels(mytable$marker), pch = c(1,2,3,4), cex=2, pt.lwd=2)
legend("top", legend = c("Mixture of 2", "Mixture of 3"), text.col = c(2,4), cex=2)


plot(mytable$prop, mytable$missCopy/18571, ylab= "Genotype error rate", xlab="Strain proportion", col=mytable$color, log="y", xlim = c(0,1), type="n", cex.lab=1.5)
strains = unique(mytable$strainName)
for ( strain in strains ){
    tmpIndex = which(mytable$strainName == strain)
    tmpMarker = as.numeric(unique(mytable$marker[tmpIndex]))
    tmpColor = unique(mytable$color[tmpIndex])
    print(strain)
#    if (as.character(strain) == "PG0402-C.2"){
#        text( mean(mytable$prop[tmpIndex]), mean(mytable$missCopy[tmpIndex])/18571, label = "PG0402-C.7G8", col=tmpColor, cex=1.5)
#    }else{
        points( mean(mytable$prop[tmpIndex]), mean(mytable$missCopy[tmpIndex])/18571, pch = tmpMarker, col=tmpColor, cex=2, lwd=2)
#    }
}
legend("topright", legend = levels(mytable$marker), pch = c(1,2,3,4), cex=2, pt.lwd=2)
legend("top", legend = c("Mixture of 2", "Mixture of 3"), text.col = c(2,4), cex=2)

dev.off()


