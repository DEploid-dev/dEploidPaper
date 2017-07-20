#!/usr/bin/env Rscript
rm(list=ls());
#library(DEploid)
source("~/DEploid/utilities/DEploidR.R")
#args = (commandArgs(TRUE))

#experimentID = as.numeric(args[1])
#experimentID = 1 #as.numeric(args[1])


write.coverage <- function(p1, p2){
    WSAF1 = hap1 * p1/100 + hap2 * p2/100
    includeErrorWSAF1 = WSAF1*(1-err)+(1-WSAF1)*err
    altCount1 = rbinom(n.loci, tmpTotalCoverage, includeErrorWSAF1)
    refCount1 = tmpTotalCoverage - altCount1
    write.table(data.frame(CHROM = coverage$CHROM, POS = coverage$POS, REF = refCount1),
        file = paste("experimentCoverage20/",prefix, ".", p1, "v", p2, ".ref", sep = ""), row.names = F, quote=F, col.names=T, sep="\t")
    write.table(data.frame(CHROM = coverage$CHROM, POS = coverage$POS, ALT = altCount1),
        file = paste("experimentCoverage20/",prefix, ".", p1, "v", p2, ".alt", sep = ""), row.names = F, quote=F, col.names=T, sep="\t")
}

groups = c("asia1",
           "asia2",
           "asia3",
           "africa1",
           "africa2",
           "africa3",
           "africa4")

is.0s = c()
for ( group in groups ){

#    prefix = paste(group, "experiment", experimentID, sep ="")
#    set.seed(experimentID)
    #chrom14 = "Pf3D7_01_v3"
    fullPanelFile = read.csv(paste("panel14/", group, "_Pf3D7_14_v3.csv", sep=""), header=T, stringsAsFactors = F, check.names = F)
    err<-0.01;

    fullSampleNames = names(fullPanelFile)[-c(1,2)]

    useSamples = sample(fullSampleNames, 22, replace = F)

    plafFile = read.table(paste("plaf14/", group, "_Pf3D7_14_v3_PLAF.txt", sep=""), header=T, stringsAsFactors = F)

#    plaf.is.0.at = which(plafFile$PLAF == 0)
is.0s = cbind(is.0s, plafFile$PLAF == 0)
#    write.table(data.frame(CHROM = plafFile$CHROM[plaf.is.0.at],
#        POS = plafFile$POS[plaf.is.0.at]),
#        file = paste(group, "excludeAt.txt", sep=""), row.names = F, quote = F, sep = "\t")

}

at.least.oneOf.plaf.is.0.at = which(rowSums(is.0s) > 0)

write.table(data.frame(CHROM = plafFile$CHROM[at.least.oneOf.plaf.is.0.at],
    POS = plafFile$POS[at.least.oneOf.plaf.is.0.at]),
    file = paste("excludeAt.txt", sep=""), row.names = F, quote = F, sep = "\t")

