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
    system(paste("bgzip", newVcfName))
}

set.seed(1)

sampleName = "PG0406-C"
originalVcf = paste("../validation/", sampleName, ".eg.vcf", sep="")

coverage = fun.extract.vcf(originalVcf)

expectedTotalCov = c(150, 80, 30)
for ( i in expectedTotalCov ){
    WSAF = coverage$altCount / (coverage$altCount+coverage$refCount+0.000001)
    tmpCoverage = coverage
    tmpTotalCoverage = rpois(length(coverage$altCount), i);
    tmpCoverage$altCount = rbinom(length(coverage$altCount), tmpTotalCoverage, WSAF)
    tmpCoverage$refCount = tmpTotalCoverage - tmpCoverage$altCount
    tmpNewVcfName = paste(sampleName, ".subSample.expectedCov", i, ".vcf", sep="")
    fun.writeSubVcf (originalVcf, tmpNewVcfName, tmpCoverage)
}
