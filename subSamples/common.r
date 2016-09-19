fun.extract.vcf <- function ( vcfName, ADFieldIndex = 2 ){
    # Assume that AD is the second field
    catCmd = "cat"
    if ( grepl("gzip", system(paste("file --mime-type", vcfName), T) ) == TRUE ){
        catCmd = "zcat"
    }

    skipNum = as.numeric(system(paste(catCmd, vcfName, " | head -5000 | grep \"##\" | wc -l"), T))
    vcf  = read.table( gzfile(vcfName), skip=skipNum, header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)

    sampleName = names(vcf)[10]

    tmp = vcf[[sampleName]]
    field = strsplit(as.character(tmp),":")

    tmpCovStr = unlist(lapply(field, `[[`, ADFieldIndex))
    tmpCov = strsplit(as.character(tmpCovStr),",")

    refCount = as.numeric(unlist(lapply(tmpCov, `[[`, 1)))
    altCount = as.numeric(unlist(lapply(tmpCov, `[[`, 2)))

    return ( data.frame( CHROM = vcf[,1],
                         POS = vcf[,2],
                         refCount = refCount,
                         altCount = altCount )
           )
}


