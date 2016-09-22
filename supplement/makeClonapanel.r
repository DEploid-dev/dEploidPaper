rm ( list = ls() )
library(scatterplot3d)

load("labStrains.v2.cov.Rdata")

n.loci = dim(cov.alt)[1]
threeD7 = rep(0, n.loci)

CHROM = cov.alt$CHROM
POS = cov.alt$POS

dd2gt.from.regression = rep(0, n.loci)

#> names(cov.alt)
# [1] "CHROM"     "POS"       "PG0389-C"  "PG0390-C"  "PG0391-C"  "PG0392-C"
# [7] "PG0393-C"  "PG0394-C"  "PG0395-C"  "PG0396-C"  "PG0397-C"  "PG0398-C"
#[13] "PG0399-C"  "PG0400-C"  "PG0401-C"  "PG0402-C"  "PG0403-C"  "PG0404-C"
#[19] "PG0405-C"  "PG0406-C"  "PG0407-C"  "PG0408-C"  "PG0409-C"  "PG0410-C"
#[25] "PG0411-C"  "PG0412-C"  "PG0413-C"  "PG0414-C"  "PG0415-C"  "PG0051-C"
#[31] "PG0052-C"  "PG0008-CW" "PG0004-CW" "7G8"

samples.with.dd2 = c(3:8, 32)
for ( marker_i in c(1:n.loci) ){
    tmpy = as.numeric(cov.alt[marker_i,samples.with.dd2])
    if (sum(tmpy) > 0){
        tmpx = as.numeric(c(0.1, 0.2, 0.33, 0.67, 0.8, 0.9, 1)*(cov.alt[marker_i,samples.with.dd2] + cov.ref[marker_i,samples.with.dd2]))
        mylm = lm(tmpy~tmpx)
        tmpSummary = summary(mylm)
        if ( runif(1)<0.001 ){
            mymain = paste("Coverage ALT/REF at", CHROM[marker_i], POS[marker_i] )

            png(paste("dd2marker",marker_i,".png",sep=""))
            plot(tmpx, tmpy, xlab = "ALT+REF", ylab="ALT", main = mymain)
            abline(mylm)
            dev.off()

        }

        if ( dim(tmpSummary$coefficients)[1] < 2 | is.nan(tmpSummary$coefficients[2,4]) ){
            next
        }
        if ( tmpSummary$coefficients[2,1] > 0 && tmpSummary$coefficients[2,4] < 0.001 ){ # p-value
            print(marker_i)
            dd2gt.from.regression[ marker_i ] = 1
        }
    }
}

samples.with.hb3.7g8 = c (12:29, 31, 33, 34)

HB3gt.from.regression = rep(0, n.loci)
sevenG8gt.from.regression = rep(0, n.loci)

tmp.cov.alt = cov.alt[,samples.with.hb3.7g8]
tmp.cov.ref = cov.ref[,samples.with.hb3.7g8]

prop_array = c(1, .99, .95, .90, .85, .80, .75, .70, .6, .5, .4, .3, .25, .2, .15, .05, .01, 0, 1, 1, 0) # Proportion of Hb3

assertthat::are_equal( length(samples.with.hb3.7g8), length(prop_array) )

for ( marker_i in c(1:n.loci) ){
    tmpy = as.numeric(tmp.cov.alt[marker_i,])
    if (sum(tmpy) > 0){
        x = as.numeric(tmp.cov.ref[marker_i,]+tmp.cov.alt[marker_i,])
        x1 = x*prop_array
        x2 = x*(1-prop_array)
        mylm = lm(tmpy~x1+x2)
        tmpSummary = summary(mylm)

        if ( runif(1)<0.001 ){
            mymain = paste("Coverage ALT/REF at", CHROM[marker_i], POS[marker_i])
            png(paste("marker",marker_i,".png",sep=""))
            s3d <-scatterplot3d(x1,x2,tmpy, pch=16, highlight.3d=TRUE,
                type="h", main=mymain, xlab = "HB3", ylab = "7G8", zlab = "ALT")
            s3d$plane3d(mylm)
            dev.off()
        }

        if ( dim(tmpSummary$coefficients)[1] < 3 ){
            next
        }
        if ( tmpSummary$coefficients[2,1] > 0 && tmpSummary$coefficients[2,4] < 0.001 ){ # p-value
            print(marker_i)
            HB3gt.from.regression[ marker_i ] = 1
        }
        if ( tmpSummary$coefficients[3,1] > 0 && tmpSummary$coefficients[3,4] < 0.001 ){ # p-value
            print(marker_i)
            sevenG8gt.from.regression[ marker_i ] = 1
        }
    }
}

combinedPanel = cbind(CHROM, POS, threeD7, dd2gt.from.regression, HB3gt.from.regression, sevenG8gt.from.regression)

write.csv(combinedPanel, file = "labStrainsPanelFinal.csv", quote = F, row.names = FALSE)

