rm(list = ls())
true.prop = read.csv("true.prop.4strains.csv", header=T)
deploid.prop = read.csv("deploid.prop.4strains.csv", header=T)
pfmix.prop = read.csv("pfmix.prop.4strains.csv", header=T)

actualPlot <- function (){
    for ( j in 1:2 ){
        points(true.prop[,j][1:6], deploid.prop[,j][1:6], pch = j, col = "red", cex=1.4,lwd=2)
        points(true.prop[,j][1:6], pfmix.prop[,j][1:6], pch = j, col = "blue", cex=1.4, lwd=2)
    }

    for ( j in 3:4 ){
        points(true.prop[,j][7:24], deploid.prop[,j][7:24], pch = j, col = "red", cex=1.4, lwd=2)
        points(true.prop[,j][7:24], pfmix.prop[,j][7:24], pch = j, col = "blue", cex=1.4, lwd=2)
    }

    for ( j in 2:4 ){
        points(true.prop[,j][25:27], deploid.prop[,j][25:27], pch = j, col = "red", cex=1.4,lwd=2)
        points(true.prop[,j][25:27], pfmix.prop[,j][25:27], pch = j, col = "blue", cex=1.4, lwd=2)
    }

    legend("bottomright", legend = c("3D7", "Dd2",  "HB3",  "X7G8"), pch = c(1,2,3,4), cex=1.4, pt.lwd=2)
    legend("topleft", legend = c("DEploid", "pfmix"), text.col = c("red", "blue"), cex=1.4)


}

pdf("trueVsInferred4.pdf", width = 12, height = 12)
par(mfrow = c(2,2))
plot(c(0,0.23), c(0,0.23), type="n", xlab = "True proportions", ylab = "Inferred proportions", cex.lab= 1.4)
lines(c(0.001,1), c(0.001,1), lty=2)
actualPlot()


plot(c(0.23,0.42), c(0.23,0.42), type="n", xlab = "True proportions", ylab = "Inferred proportions", cex.lab= 1.4)
lines(c(0.001,1), c(0.001,1), lty=2)
actualPlot()


plot(c(0.48,0.77), c(0.48,0.77), type="n", xlab = "True proportions", ylab = "Inferred proportions", cex.lab= 1.4)
lines(c(0.001,1), c(0.001,1), lty=2)
actualPlot()


plot(c(0.78,1), c(0.78,1), type="n", xlab = "True proportions", ylab = "Inferred proportions", cex.lab= 1.4)
lines(c(0.001,1), c(0.001,1), lty=2)
actualPlot()
dev.off()

pdf("trueVsInferred.pdf", width = 8, height = 8)
plot(c(0,1), c(0,1), type="n", xlab = "True proportions", ylab = "Inferred proportions", cex.lab= 1.4)
lines(c(0,1), c(0,1), lty=2)
lines(c(0,1), c(0.02,1.02), lty=2, col="grey", lwd = .5)
lines(c(0,1), c(-0.02,0.98), lty=2, col="grey", lwd = .5)
actualPlot()
dev.off()
