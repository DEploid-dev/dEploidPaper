rm(list=ls())

compute.effective.k <- function (prop){
    return(1/sum(prop^2))
}

jitterOff = 50

true.prop = list()
# 3d7, dd2, 7g8, hb3
true.prop[["PG0389-C"]] = c(0.9, 0.1, 0, 0)
true.prop[["PG0390-C"]] = c(0.8, 0.2, 0, 0)
true.prop[["PG0391-C"]] = c(0.67, 0.33, 0, 0)
true.prop[["PG0392-C"]] = c(0.33, 0.67, 0, 0)
true.prop[["PG0393-C"]] = c(0.2, 0.8,0,0)
true.prop[["PG0394-C"]] = c(0.1, 0.9,0,0)
true.prop[["PG0395-C"]] = c(0, 0.33, 0.34, 0.33)
true.prop[["PG0396-C"]] = c(0, 0.25, 0.5, 0.25)
true.prop[["PG0397-C"]] = c(0, 0.143, 0.714, 0.143)
true.prop[["PG0398-C"]] = c(0, 0, 0, 1)
true.prop[["PG0399-C"]] = c(0, 0, 0.01, 0.99)
true.prop[["PG0400-C"]] = c(0, 0, 0.05, 0.95)
true.prop[["PG0401-C"]] = c(0, 0, 0.1, 0.9)
true.prop[["PG0402-C"]] = c(0, 0, 0.15, 0.85)
true.prop[["PG0403-C"]] = c(0, 0, 0.2, 0.8)
true.prop[["PG0404-C"]] = c(0, 0, 0.25, 0.75)
true.prop[["PG0405-C"]] = c(0, 0, 0.3, 0.7)
true.prop[["PG0406-C"]] = c(0, 0, 0.4, 0.6)
true.prop[["PG0407-C"]] = c(0, 0, 0.5, 0.5)
true.prop[["PG0408-C"]] = c(0, 0, 0.6, 0.4)
true.prop[["PG0409-C"]] = c(0, 0, 0.7, 0.3)
true.prop[["PG0410-C"]] = c(0, 0, 0.75, 0.25)
true.prop[["PG0411-C"]] = c(0, 0, 0.80, 0.2)
true.prop[["PG0412-C"]] = c(0, 0, 0.85, 0.15)
true.prop[["PG0413-C"]] = c(0, 0, 0.95, 0.05)
true.prop[["PG0414-C"]] = c(0, 0, 0.99, 0.01)
true.prop[["PG0415-C"]] = c(0, 0, 1,0)

sampleNames = read.table("../labSampleNames", header=F, stringsAsFactors = F)$V1

pdf("eff_k.pdf", width = 8, height = 12)

#sampleNames = c("PG0395-C")
par(mar=c(5.1,21.1,2.1,2.1))
plot(c(0,1), c(0,1), type="n", ylim = c(1, 27), xlim = c(0.8,3.5), xlab = "Effective K", ylab = "", yaxt="n", cex.lab=2, cex=2, cex.axis=1.4)

sampleI = 1
true.eff.k = c()
#eff.2 = c()
#eff.3 = c()
#eff.4 = c()
#eff.5 = c()
my_y_lab = c()
cat ("sample p1 p2 p3 p4 p5\n", file = "deploid.prop.csv", append = F)

seeds = read.table("mostTwoSwitchOne.txt", header=T, stringsAsFactors = T)
for (sampleName in sampleNames){
    myk = 2
    true.eff.k = c(true.eff.k, compute.effective.k(true.prop[[sampleName]]))
    a = read.table(paste(sampleName, "k", myk, ".prop", sep=""), header=F)
    for (i in 1:dim(a)[1]){
        prop = a[i,]
        strains = which(prop>0.01)
        points(compute.effective.k(prop[strains]), sampleI - 0.2, col = adjustcolor( "blue", alpha.f = 0.1), pch=16, cex=1.4)
    }

    myk = 3
    a = read.table(paste(sampleName, "k", myk, ".prop", sep=""), header=F)
    for (i in 1:dim(a)[1]){
        prop = a[i,]
        strains = which(prop>0.01)
        points(compute.effective.k(prop[strains]), sampleI - 0.1, col = adjustcolor( "yellow", alpha.f = 0.1), pch=16, cex=1.4)
    }

    myk = 4
    true.eff.k = c(true.eff.k, compute.effective.k(true.prop[[sampleName]]))
    a = read.table(paste(sampleName, "k", myk, ".prop", sep=""), header=F)
    for (i in 1:dim(a)[1]){
        prop = a[i,]
        strains = which(prop>0.01)
        points(compute.effective.k(prop[strains]), sampleI + 0.1, col = adjustcolor( "green", alpha.f = 0.1), pch=16, cex=1.4)
    }

    myk = 5
    true.eff.k = c(true.eff.k, compute.effective.k(true.prop[[sampleName]]))
    a = read.table(paste(sampleName, "k", myk, ".prop", sep=""), header=F)
    for (i in 1:dim(a)[1]){
        prop = a[i,]
        strains = which(prop>0.01)
        points(compute.effective.k(prop[strains]), sampleI + 0.2, col = adjustcolor( "red", alpha.f = 0.1), pch=16, cex=1.4)
        if ( i == seeds[sampleI,2] ){
#            points(compute.effective.k(prop[strains]), sampleI + 0.2, col = "red", pch="X", cex=1.4)
            cat (paste(sampleName, paste(a[i,], collapse = " "), sep = " "), file = "deploid.prop.csv", append = T)
            cat ("\n", file = "deploid.prop.csv", append = T)
        }
    }


#jitter(sampleI, jitterOff)
    points(compute.effective.k(true.prop[[sampleName]]), sampleI, col = "black", pch="X", cex=1.4)
    printingProp = paste( format(true.prop[[sampleName]], nsmall = 2, digits =2), collapse = ", ")
    my_y_lab = c(my_y_lab, paste(sampleName, "(",printingProp, ")"))
#    text(paste(sampleName, "(",printingProp, ")"), x = 0, y = sampleI)
    sampleI = sampleI + 1
}
my_y_lab = c(my_y_lab, paste("         3D7, Dd2, HB3, 7G8  "))
axis(2, at=1:28,labels=my_y_lab, las=2, lwd = 0, cex=2, cex.axis=1.4)
legend("topright", legend = paste("k=", c(2:5)), col = c("blue", "yellow", "green", "red"), pch = 16, cex=1.4)
dev.off()
