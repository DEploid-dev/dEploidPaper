rm(list=ls())
source("common.r")

compute.errors <- function(population, case, experiments = c(1:100)){
    excludeAt = read.table(paste(population,"excludeAt.txt", sep=""), header = T, stringsAsFactors = F)

    switches = c()
    muts = c()
    dropError = c()
    for ( experimentID in experiments ){
        print(paste("experiment", experimentID))
    #    hap = read.table(paste("experimentOut/experiment", experimentID, ".25v75.hap", sep=""), header=T, stringsAsFactors = F)
        hap = read.table(paste("experimentOut/", population, "experiment", experimentID, ".", case, "panel20.hap", sep=""), header=T, stringsAsFactors = F)
#experimentOut/asia1experiment99.45v55panel20.hap
        truths = read.table(paste("truth20/", population, "experiment", experimentID, ".true.hap", sep = ""), header=T, stringsAsFactors = F)

        high.lighting = which(! truths$POS %in% excludeAt$POS)
        truths = truths[high.lighting,]
         a =  measure.error.joe.with.drop.2strain (t(hap[,3:4]), t(truths[,3:4]), rel.cost.switch=2)

        switches = c(switches, sum(a$switchError))
        muts = c(muts, sum(a$mutError)/dim(hap)[1])
        dropError = c(dropError, a$dropError / dim(hap)[1])
    }

    return (list(switches = switches,
                 muts = muts,
                 dropError = dropError) )
}



case.asia1.15 = compute.errors("asia1", "15v85", c(1:100))
case.asia1.25 = compute.errors("asia1", "25v75", c(1:100))
case.asia1.35 = compute.errors("asia1", "35v65", c(1:100))
case.asia1.45 = compute.errors("asia1", "45v55", c(1:100))
save("case.asia1.15", "case.asia1.25", "case.asia1.35", "case.asia1.45", file = "asia1.errors.Rdata")



case.asia2.15 = compute.errors("asia2", "15v85", c(1:100))
case.asia2.25 = compute.errors("asia2", "25v75", c(1:100))
case.asia2.35 = compute.errors("asia2", "35v65", c(1:100))
case.asia2.45 = compute.errors("asia2", "45v55", c(1:100))
save("case.asia2.15", "case.asia2.25", "case.asia2.35", "case.asia2.45", file = "asia2.errors.Rdata")


case.asia3.15 = compute.errors("asia3", "15v85", c(1:100))
case.asia3.25 = compute.errors("asia3", "25v75", c(1:100))
case.asia3.35 = compute.errors("asia3", "35v65", c(1:100))
case.asia3.45 = compute.errors("asia3", "45v55", c(1:100))
save("case.asia3.15", "case.asia3.25", "case.asia3.35", "case.asia3.45", file = "asia3.errors.Rdata")



case.africa3.15 = compute.errors("africa3", "15v85", c(1:100))
case.africa3.25 = compute.errors("africa3", "25v75", c(1:100))
case.africa3.35 = compute.errors("africa3", "35v65", c(1:100))
case.africa3.45 = compute.errors("africa3", "45v55", c(1:100))
save("case.africa3.15", "case.africa3.25", "case.africa3.35", "case.africa3.45", file = "africa3.errors.Rdata")

case.africa4.15 = compute.errors("africa4", "15v85", c(1:100))
case.africa4.25 = compute.errors("africa4", "25v75", c(1:100))
case.africa4.35 = compute.errors("africa4", "35v65", c(1:100))
case.africa4.45 = compute.errors("africa4", "45v55", c(1:100))
save("case.africa4.15", "case.africa4.25", "case.africa4.35", "case.africa4.45", file = "africa4.errors.Rdata")

case.africa1.15 = compute.errors("africa1", "15v85", c(1:100))
case.africa1.25 = compute.errors("africa1", "25v75", c(1:100))
case.africa1.35 = compute.errors("africa1", "35v65", c(1:100))
case.africa1.45 = compute.errors("africa1", "45v55", c(1:100))
save("case.africa1.15", "case.africa1.25", "case.africa1.35", "case.africa1.45", file = "africa1.errors.Rdata")








computeMeansAndMedian <- function (population, case, experiments){
    means = c()
    medians = c()

    for ( experimentID in experiments ){
        print(paste("experiment", experimentID))
    #    hap = read.table(paste("experimentOut/experiment", experimentID, ".25v75.hap", sep=""), header=T, stringsAsFactors = F)
        POS = read.table(paste("experimentCoverage20/",population,"experiment", experimentID, ".", case, ".ref", sep=""), header=T, stringsAsFactors = F)[,2]
        ref = read.table(paste("experimentCoverage20/",population,"experiment", experimentID, ".", case, ".ref", sep=""), header=T, stringsAsFactors = F)[,3]
        alt = read.table(paste("experimentCoverage20/",population,"experiment", experimentID, ".", case, ".alt", sep=""), header=T, stringsAsFactors = F)[,3]
        coverage = ref+alt
#        truths = read.table(paste("truth20/experiment", experimentID, ".true.hap", sep = ""), header=T, stringsAsFactors = F)

        high.lighting = which(! POS %in% excludeAt$POS)
        coverage = coverage[high.lighting]
        means = c(means, mean(coverage))
        medians = c(medians, median(coverage))
    }
    return(list(means =means, medians=medians))
}

coverage.asia1.15 = computeMeansAndMedian("asia1", "15v85", c(1:100))
