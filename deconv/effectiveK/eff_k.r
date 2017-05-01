rm(list=ls())

compute.effective.k <- function (prop){
    return(1/sum(prop^2))
}

jitterOff = 50

true.prop = list()
true.prop[["PG0389-C"]] = c(0.1, 0.9)
true.prop[["PG0390-C"]] = c(0.2, 0.8)
true.prop[["PG0391-C"]] = c(0.33, 0.67)
true.prop[["PG0392-C"]] = c(0.33, 0.67)
true.prop[["PG0393-C"]] = c(0.2, 0.8)
true.prop[["PG0394-C"]] = c(0.1, 0.9)
true.prop[["PG0395-C"]] = c(0.33, 0.33, 0.34)
true.prop[["PG0396-C"]] = c(0.25, 0.25, 0.5)
true.prop[["PG0397-C"]] = c(0.143, 0.143, 0.714)
true.prop[["PG0398-C"]] = c(1)
true.prop[["PG0399-C"]] = c(0.01, 0.99)
true.prop[["PG0400-C"]] = c(0.05, 0.95)
true.prop[["PG0401-C"]] = c(0.1, 0.9)
true.prop[["PG0402-C"]] = c(0.15, 0.85)
true.prop[["PG0403-C"]] = c(0.2, 0.8)
true.prop[["PG0404-C"]] = c(0.25, 0.75)
true.prop[["PG0405-C"]] = c(0.3, 0.7)
true.prop[["PG0406-C"]] = c(0.4, 0.6)
true.prop[["PG0407-C"]] = c(0.5, 0.5)
true.prop[["PG0408-C"]] = c(0.6, 0.4)
true.prop[["PG0409-C"]] = c(0.7, 0.3)
true.prop[["PG0410-C"]] = c(0.75, 0.25)
true.prop[["PG0411-C"]] = c(0.80, 0.2)
true.prop[["PG0412-C"]] = c(0.85, 0.15)
true.prop[["PG0413-C"]] = c(0.95, 0.05)
true.prop[["PG0414-C"]] = c(0.99, 0.01)
true.prop[["PG0415-C"]] = c(1)

sampleNames = read.table("../labSampleNames", header=F, stringsAsFactors = F)$V1

pdf("eff_k.pdf", width = 8, height = 12)

#sampleNames = c("PG0395-C")
par(mar=c(6.1,20.1,4.1,2.1))
plot(c(0,1), c(0,1), type="n", ylim = c(1, 27), xlim = c(0.8,3.5), xlab = "Effective K", ylab = "", yaxt="n", cex.lab=2, cex=2, cex.axis=1.4)

sampleI = 1
true.eff.k = c()
my_y_lab = c()
for (sampleName in sampleNames){
    myk = 5
    true.eff.k = c(true.eff.k, compute.effective.k(true.prop[[sampleName]]))
    a = read.table(paste(sampleName, "k", myk, ".prop", sep=""), header=F)
    for (i in 1:dim(a)[1]){
        prop = a[i,]
        strains = which(prop>0.01)
        points(compute.effective.k(prop[strains]), sampleI + 0.2, col = adjustcolor( "red", alpha.f = 0.1), pch=16, cex=1.4)
    }

    myk = 4
    true.eff.k = c(true.eff.k, compute.effective.k(true.prop[[sampleName]]))
    a = read.table(paste(sampleName, "k", myk, ".prop", sep=""), header=F)
    for (i in 1:dim(a)[1]){
        prop = a[i,]
        strains = which(prop>0.01)
        points(compute.effective.k(prop[strains]), sampleI + 0.1, col = adjustcolor( "green", alpha.f = 0.1), pch=16, cex=1.4)
    }

    myk = 3
    a = read.table(paste(sampleName, "k", myk, ".prop", sep=""), header=F)
#    ef.k = c()
    for (i in 1:dim(a)[1]){
        prop = a[i,]
        strains = which(prop>0.01)
        points(compute.effective.k(prop[strains]), sampleI - 0.1, col = adjustcolor( "yellow", alpha.f = 0.1), pch=16, cex=1.4)
    }

    myk = 2
    true.eff.k = c(true.eff.k, compute.effective.k(true.prop[[sampleName]]))
    a = read.table(paste(sampleName, "k", myk, ".prop", sep=""), header=F)
    for (i in 1:dim(a)[1]){
        prop = a[i,]
        strains = which(prop>0.01)
        points(compute.effective.k(prop[strains]), sampleI - 0.2, col = adjustcolor( "blue", alpha.f = 0.1), pch=16, cex=1.4)
    }

#jitter(sampleI, jitterOff)
    points(compute.effective.k(true.prop[[sampleName]]), sampleI, col = "black", pch="X", cex=1.4)
    printingProp = paste( sort(true.prop[[sampleName]]), collapse = ", ")
    my_y_lab = c(my_y_lab, paste(sampleName, "(",printingProp, ")"))
#    text(paste(sampleName, "(",printingProp, ")"), x = 0, y = sampleI)
    sampleI = sampleI + 1
}

axis(2, at=1:27,labels=my_y_lab, las=2, lwd = 0, cex=2, cex.axis=1.4)
legend("topright", legend = paste("k=", c(2:5)), col = c("blue", "yellow", "green", "red"), pch = 16, cex=1.4)
dev.off()
