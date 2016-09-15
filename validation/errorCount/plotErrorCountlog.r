rm(list=ls())
png("switchVsMisCopyErrlog.png", width=1500, height = 1500)
plot(c(10,300), c(1150,3200), type="n", log="xy", xlab="Number of Switch Error", ylab ="Number of Miss Copying Error")
#plot(c(0,110), c(50,1150), type="n")

sample3 = read.table("../labSampleNames3Strains", header=F, stringsAsFactors=F)$V1
for (samplei in sample3){
    tmpErrors = c()
    tmpSeed=c()
    for ( seed in 1:15 ){
        prefix = paste(samplei, "_seed", seed, "k3", sep="")
        fileName = paste(prefix, ".errorCount", sep ="")
        if ( file.exists(fileName) ){
            a = system(paste("tail -1 ", prefix, ".log",sep=""), intern=T)
            prop = as.numeric(unlist(strsplit(a, split="\t")))
            nStrain = sum(which(prop>0.01))
#            cat (seed, " ", prop, "\n")
            tmpErrors = cbind(tmpErrors, read.table(fileName, header=F)$V1)
            tmpSeed=c(tmpSeed, seed)
        }
    }
    centre = rowMeans(tmpErrors)
#    points(centre[1], centre[2])
    errorBar = apply(tmpErrors,1,sd)
#    lines(c(centre[1]-errorBar[1], centre[1]+errorBar[1]), c(centre[2],centre[2]))
#    lines(c(centre[1], centre[1]), c(centre[2]-errorBar[2],centre[2]+errorBar[2]))
    text(centre[1], centre[2], label=samplei, col="blue")
}


sample2 = read.table("../labSampleNames2Strains", header=F, stringsAsFactors=F)$V1
for (samplei in sample2){
    tmpErrors = c()
    tmpSeed=c()
    for ( seed in 1:15 ){
        prefix = paste(samplei, "_seed", seed, "k2", sep="")
        fileName = paste(prefix, ".errorCount", sep ="")
        if ( file.exists(fileName) ){
            a = system(paste("tail -1 ", prefix, ".log",sep=""), intern=T)
            prop = as.numeric(unlist(strsplit(a, split="\t")))
            nStrain = sum(which(prop>0.01))
#            cat (seed, " ", prop, "\n")
            tmpErrors = cbind(tmpErrors, read.table(fileName, header=F)$V1)
            tmpSeed=c(tmpSeed, seed)
        }
    }
    centre = rowMeans(tmpErrors)
#    points(centre[1], centre[2])
#    points(tmpErrors[1,], tmpErrors[2,])
#    text(tmpErrors[1,], tmpErrors[2,], labels=tmpSeed)
    errorBar = apply(tmpErrors,1,sd)
#    lines(c(centre[1]-errorBar[1], centre[1]+errorBar[1]), c(centre[2],centre[2]))
#    lines(c(centre[1], centre[1]), c(centre[2]-errorBar[2],centre[2]+errorBar[2]))
    if ( samplei %in% c("PG0389-C", "PG0394-C", "PG0398-C", "PG0399-C", "PG0400-C", "PG0401-C", "PG0413-C", "PG0414-C", "PG0415-C") ){
        text(centre[1], centre[2], label=samplei, col = "red")
    } else {
        text(centre[1], centre[2], label=samplei)
    }
}
dev.off()


