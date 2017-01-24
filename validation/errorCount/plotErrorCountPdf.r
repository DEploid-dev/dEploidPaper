rm(list=ls())
#pdf("switchVsMisCopyErrlogHi.pdf", width=8, height = 8)
pdf("HapErrors.pdf", width=8, height = 8)


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
deploid2Color= 2
deploid3Color = 6
shapeitColor = 3
beagleColor = 4
mytable[["color"]] = deploid2Color
mytable$color[mytable$sampleName %in% c("PG0395-C","PG0396-C","PG0397-C")] = deploid3Color
mytable[["marker"]] = factor(mytable$labStrain)


mytable2 = data.frame ( prop = numeric(0), coverage = numeric(0), switchError = numeric(0), missCopy = numeric(0), sampleName = character(), strainName = character(), labStrain = character(), stringsAsFactors = F)
mytablei = 1
print(simpleMix)
for (samplei in sample2){
    fileName = paste("../../benchMark/beagle/", samplei, ".errorCount", sep ="")
    if ( file.exists(fileName) ){
        tmpInfo = read.table(fileName, header=F)
        for ( i in 1:dim(tmpInfo)[1] ){
            mytable2[mytablei,]$prop = tmpInfo$V1[i]
            mytable2[mytablei,]$coverage = metaData$mean_coverage[metaData$sample == samplei]
            mytable2[mytablei,]$switchError = tmpInfo$V2[i]
            mytable2[mytablei,]$missCopy = tmpInfo$V3[i]
            mytable2[mytablei,]$sampleName = samplei
            mytable2[mytablei,]$strainName = paste(samplei,".", i, sep="")
            labStrain = ""
            if ( i == 1 ){
                labStrain = "hb3"
                if ( samplei %in% simpleMix) { labStrain = "3d7" }
            } else {
                labStrain = "7g8"
                if ( samplei %in% simpleMix) { labStrain = "dd2" }
            }
            mytable2[mytablei,]$labStrain = labStrain
            mytable2[mytablei,]$color = 5
            mytablei = mytablei+1
        }
    }
}

mytable2[["marker"]] = factor(mytable2$labStrain)


mytable3 = data.frame ( prop = numeric(0), coverage = numeric(0), switchError = numeric(0), missCopy = numeric(0), sampleName = character(), strainName = character(), labStrain = character(), stringsAsFactors = F)
mytablei = 1
print(simpleMix)
for (samplei in sample2){
    fileName = paste("../../benchMark/shapeit/", samplei, ".errorCount", sep ="")
    if ( file.exists(fileName) ){
        tmpInfo = read.table(fileName, header=F)
        for ( i in 1:dim(tmpInfo)[1] ){
            mytable3[mytablei,]$prop = tmpInfo$V1[i]
            mytable3[mytablei,]$coverage = metaData$mean_coverage[metaData$sample == samplei]
            mytable3[mytablei,]$switchError = tmpInfo$V2[i]
            mytable3[mytablei,]$missCopy = tmpInfo$V3[i]
            mytable3[mytablei,]$sampleName = samplei
            mytable3[mytablei,]$strainName = paste(samplei,".", i, sep="")
            labStrain = ""
            if ( i == 1 ){
                labStrain = "hb3"
                if ( samplei %in% simpleMix) { labStrain = "3d7" }
            } else {
                labStrain = "7g8"
                if ( samplei %in% simpleMix) { labStrain = "dd2" }
            }
            mytable3[mytablei,]$labStrain = labStrain
            mytable3[mytablei,]$color = 5
            mytablei = mytablei+1
        }
    }
}

mytable3[["marker"]] = factor(mytable3$labStrain)


layout(rbind(1,2,3), heights=c(10,10,1))  # put legend on bottom 1/8th of the chart

#par(mfrow = c(2,1))

plot(mytable$prop, mytable$switchError, ylab= "Number of switches when copying", xlab="Strain proportions", col=mytable$color, log="y", ylim= c(1, 700), xlim = c(0,1.05), type="n", cex.lab=1.6)

strains = unique(mytable2$strainName)
for ( strain in strains ){
    tmpIndex = which(mytable2$strainName == strain)

    tmpMarker = as.numeric(unique(mytable2$marker[tmpIndex]))
    tmpColor = unique(mytable2$color[tmpIndex])
    if (mytable2$prop[tmpIndex] > 0){
        points( mean(mytable2$prop[tmpIndex]), mean(mytable2$switchError[tmpIndex]), pch = tmpMarker, col=beagleColor, cex=1.1, lwd=2)
    } else {
        print(mytable2[tmpIndex,])
    }
}

strains = unique(mytable3$strainName)
for ( strain in strains ){
    tmpIndex = which(mytable3$strainName == strain)

    tmpMarker = as.numeric(unique(mytable3$marker[tmpIndex]))
    tmpColor = unique(mytable3$color[tmpIndex])
    if (mytable3$prop[tmpIndex] > 0){
        points( mean(mytable3$prop[tmpIndex]), mean(mytable3$switchError[tmpIndex]), pch = tmpMarker, col=shapeitColor, cex=1.1, lwd=2)
    } else {
        print(mytable3[tmpIndex,])
    }
}

strains = unique(mytable$strainName)
for ( strain in strains ){
    tmpIndex = which(mytable$strainName == strain)
    tmpMarker = as.numeric(unique(mytable$marker[tmpIndex]))
    tmpColor = unique(mytable$color[tmpIndex])
    points( mean(mytable$prop[tmpIndex]), mean(mytable$switchError[tmpIndex]), pch = tmpMarker, col=tmpColor, cex=1.1, lwd=2)
}


legend("topright", legend = levels(mytable$marker), pch = c(1,2,3,4), cex=1.4, pt.lwd=2)
#legend("top", legend = c("Beagle", "Shapeit", "DEploid (Mixture of 2)", "DEploid (Mixture of 3)"), text.col = c(3,8,2,4), cex=1.4, ncol=2)
legend("top", legend = c("DEploid (Mixture of 2)", "DEploid (Mixture of 3)", "BEAGLE (Mixture of 2)", "SHAPEIT (Mixture of 2)"), text.col = c(deploid2Color,deploid3Color,beagleColor,shapeitColor), cex=1.4, ncol=2)


plot(mytable$prop, mytable$missCopy/18571, ylab= "Genotype error rate", xlab="Strain proportions", col=mytable$color, log="y", ylim= c(0.02, 0.5), xlim = c(0,1.05), type="n", cex.lab=1.6)

strains = unique(mytable3$strainName)
for ( strain in strains ){
    tmpIndex = which(mytable3$strainName == strain)

    tmpMarker = as.numeric(unique(mytable3$marker[tmpIndex]))
    tmpColor = unique(mytable3$color[tmpIndex])
    if (mytable3$prop[tmpIndex] > 0){
        points( mean(mytable3$prop[tmpIndex]), mean(mytable3$missCopy[tmpIndex])/18571, pch = tmpMarker, col=shapeitColor, cex=1.1, lwd=2)
    } else {
        print(mytable3[tmpIndex,])
    }
}

strains = unique(mytable2$strainName)
for ( strain in strains ){
    tmpIndex = which(mytable2$strainName == strain)

    tmpMarker = as.numeric(unique(mytable2$marker[tmpIndex]))
    tmpColor = unique(mytable2$color[tmpIndex])
    if (mytable2$prop[tmpIndex] > 0){
        points( mean(mytable2$prop[tmpIndex]), mean(mytable2$missCopy[tmpIndex])/18571, pch = tmpMarker, col=beagleColor, cex=1.1, lwd=2)
    } else {
        print(mytable2[tmpIndex,])
    }
}

strains = unique(mytable$strainName)
for ( strain in strains ){
    tmpIndex = which(mytable$strainName == strain)
    tmpMarker = as.numeric(unique(mytable$marker[tmpIndex]))
    tmpColor = unique(mytable$color[tmpIndex])
#    print(strain)
#    if (as.character(strain) == "PG0402-C.2"){
#        text( mean(mytable$prop[tmpIndex]), mean(mytable$missCopy[tmpIndex])/18571, label = "PG0402-C.7G8", col=tmpColor, cex=1.5)
#    }else{
        points( mean(mytable$prop[tmpIndex]), mean(mytable$missCopy[tmpIndex])/18571, pch = tmpMarker, col=tmpColor, cex=1.1, lwd=2)
#    }
}
legend("topright", legend = levels(mytable$marker), pch = c(1,2,3,4), cex=1.4, pt.lwd=2)
#legend("top", legend = c("Mixture of 2", "Mixture of 3", "Beagle", "Shapeit"), text.col = c(2,4,3,8), cex=1.4, ncol=4)
legend("top", legend = c("DEploid (Mixture of 2)", "DEploid (Mixture of 3)", "BEAGLE (Mixture of 2)", "SHAPEIT (Mixture of 2)"), text.col = c(deploid2Color,deploid3Color,beagleColor,shapeitColor), cex=1.4, ncol=2)

par(mar=c(0, 0, 0, 0))
plot.new()

#par(mar=c(0, 0, 0, 0))
#plot.new()
#legend("center", legend = levels(mytable$marker), pch = c(1,2,3,4), cex=1.4, pt.lwd=2, ncol=4)

#par(mar=c(0, 0, 0, 0))
#plot.new()
#legend("center", legend = c("Mixture of 2", "Mixture of 3", "Beagle", "Shapeit"), text.col = c(2,4,3,8), cex=1.4, ncol=4)

#par(mar=c(1, 2, 0, 2))
#plot.new()
#legend("left", legend = levels(mytable$marker), pch = c(1,2,3,4), cex=1.4, pt.lwd=2, ncol=4)

##par(mar=c(0, 0, 0, 0))
##plot.new()
#legend("right", legend = c("Mixture of 2", "Mixture of 3", "Beagle", "Shapeit"), text.col = c(2,4,3,8), cex=1.4, ncol=4)

dev.off()


