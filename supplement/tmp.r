rm ( list = ls() )
library(scatterplot3d)
set.seed(1)
load("labStrains.v2.cov.Rdata")

n.loci = dim(cov.alt)[1]
threeD7 = rep(0, n.loci)

CHROM = cov.alt$CHROM
POS = cov.alt$POS


dd2gt.from.regression = rep(0, n.loci)
HB3gt.from.regression = rep(0, n.loci)
sevenG8gt.from.regression = rep(0, n.loci)

tmp.cov.alt = cov.alt[,c(c(3:8, 32),c(12:29, 31, 33, 34))]
tmp.cov.ref = cov.ref[,c(c(3:8, 32),c(12:29, 31, 33, 34))]

prop_dd2 = c(c(0.1, 0.2, 0.33, 0.67, 0.8, 0.9, 1), rep(0,21))
prop_hb3 = c(rep(0,7), c(1, .99, .95, .90, .85, .80, .75, .70, .6, .5, .4, .3, .25, .2, .15, .05, .01, 0, 1, 1, 0)) # Proportion of Hb3
prop_7g8 = c(rep(0,7), 1-c(1, .99, .95, .90, .85, .80, .75, .70, .6, .5, .4, .3, .25, .2, .15, .05, .01, 0, 1, 1, 0))
tmp=c(1:n.loci)
#tmp=c(19409:20000)
for ( marker_i in tmp ){
    tmpy = as.numeric(tmp.cov.alt[marker_i,])
    if (sum(tmpy) > 0){
        x = as.numeric(tmp.cov.ref[marker_i,]+tmp.cov.alt[marker_i,])
        xdd2 = x*prop_dd2
        xhb3 = x*prop_hb3
        x7g8 = x*prop_7g8
        mylm = lm(tmpy~xdd2+xhb3+x7g8)
        tmpSummary = summary(mylm)

        if ( runif(1)<0.001 ){
            mymain = paste("Coverage ALT/REF at", CHROM[marker_i], POS[marker_i] )

            png(paste("dd2marker",marker_i,".png",sep=""))
            plot(xdd2, tmpy, xlab = "ALT+REF", ylab="ALT", main = mymain)
            abline(lm(tmpy~xdd2))
            dev.off()

            mymain = paste("Coverage ALT/REF at", CHROM[marker_i], POS[marker_i])
            png(paste("marker",marker_i,".png",sep=""))
            s3d <-scatterplot3d(xhb3,x7g8,tmpy, pch=16, highlight.3d=TRUE,
                type="h", main=mymain, xlab = "HB3", ylab = "7G8", zlab = "ALT")
            s3d$plane3d(lm(tmpy~xhb3+x7g8))
            dev.off()
        }

        if ( dim(tmpSummary$coefficients)[1] < 3 ){
            next
        }
        if ( tmpSummary$coefficients[2,1] > 0 && tmpSummary$coefficients[2,4] < 0.001 ){ # p-value
            print(marker_i)
            dd2gt.from.regression[ marker_i ] = 1
        }
        if ( tmpSummary$coefficients[3,1] > 0 && tmpSummary$coefficients[3,4] < 0.001 ){ # p-value
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

write.csv(combinedPanel, file = "labStrainsPanelNew.csv", quote = F, row.names = FALSE)

