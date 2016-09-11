rm(list = ls())

#groupName = "asiaGroup1"

sampleNames = as.character(read.table("labSampleNames", header = F)$V1)

root = "/well/mcvean/joezhu/pf3k/pf3k_5_1_final/dEploidOut/"

#PG0389-C_seed2k2

for (k in c(2,5)){
    suffix = paste("k", k, sep="")
    chosenSeed = c()
    printSampleNames = c()

    for ( sampleName in sampleNames ){
        dic = c()
        for ( seed in 1:15){
            tmpFileName = paste(root, sampleName, "/", sampleName, "_seed", seed, suffix, "dic.log", sep = "")
            if ( !file.exists(tmpFileName) ){
                print(paste("file:", tmpFileName, " does not exist"))
                next
            }
            dic = c (dic, read.table ( tmpFileName, header = F)$V2[3] )
        }
        if ( length(dic) > 0 ){
            chosenSeed = c( chosenSeed, which.min(dic) )
            printSampleNames = c(printSampleNames, sampleName)
        }
    }

    write.table ( cbind(printSampleNames, chosenSeed), file = paste("labSample", suffix, "seed", sep=""), row.names = F, col.names = F, sep="\t", quote = F)
}
