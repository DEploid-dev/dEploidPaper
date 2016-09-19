rm(list=ls())
source("common.r")

fun.writeSubVcf <- function( vcfName, newVcfName, coverage ){
    catCmd = "cat"
    if ( grepl("gzip", system(paste("file --mime-type", vcfName), T) ) == TRUE ){
        catCmd = "zcat"
    }

    skipNum = as.numeric(system(paste(catCmd, vcfName, " | head -5000 | grep \"##\" | wc -l"), T))
    vcf  = read.table( gzfile(vcfName), skip=skipNum, header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)

    system ( paste("grep \"##\"", vcfName, ">", newVcfName) )
    vcf[["FORMAT"]] = rep("GT:AD", length(vcf[["FORMAT"]]))
    vcf[,10] = paste( ".:", coverage$refCount, ",", coverage$altCount, sep = "")
    write.table(vcf, file = newVcfName, append = T, sep = "\t", quote = F, row.names = F)

}

set.seed(1)
sampleName = "PG0413-C"
originalVcf = paste("../validation/", sampleName, ".14.vcf", sep="")

coverage = fun.extract.vcf(originalVcf)
for ( i in 1:10 ){
    currentP = i/10
    tmpCoverage = coverage
    tmpCoverage$refCount = rbinom(length(coverage$refCount), coverage$refCount, currentP)
    tmpCoverage$altCount = rbinom(length(coverage$altCount), coverage$altCount, currentP)
    tmpNewVcfName = paste(sampleName, ".subSample", i, "0.vcf", sep="")
    fun.writeSubVcf (originalVcf, tmpNewVcfName, tmpCoverage)
}
