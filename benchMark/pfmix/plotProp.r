rm(list=ls())
pdf("tmp.pdf", width=8, height=8)
plot(c(0.01,.99),c(.99,0.01), xlab = "Proportion component 1", ylab = "Proportion component 2", type="n")#\\, log="xy"
#lines(c(0,1),c(1,0), lty=2)
#lines(c(0,1),c(0.5,0), lty=2)

pchArray = c("*", "#", "+")
names3 = read.table("labSampleNames3Strains", stringsAsFactors=F, header=F)$V1

true.prop = read.csv("true.prop.csv", header=F)
#for (i in 1:27){
for (i in 1:15){
    prop = true.prop[i,]
    if (length(which(prop>0.01))==3){
        prop = sort(prop, decreasing=T)
    }
#    print(prop)

    if ( i > 24 ){
        points(prop[1], prop[2], col="green",cex=2, pch = pchArray[i-24])
    } else {
        points(prop[1], prop[2], col="green",cex=2,lwd=2, pch=i)
    }
}


names2 = read.table("labSampleNames2Strains", stringsAsFactors=F, header=F)$V1
#for ( i in 1:24 ){
for ( i in 1:15 ){
    prop = read.table(paste(names2[i],".pfmix.2.prop", sep=""),head=F)$V1
    if ( i %in% c(1,2,3,7:16) ){
        prop = sort(prop, decreasing=T)
    }

    points(prop[1], prop[2], col="red",cex=2,lwd=2, pch = i)
}

#for (i in 1:3){
#    prop = sort(read.table(paste(names3[i],".pfmix.100.prop", sep=""),head=F)$V1, decreasing=T)
#    print(prop)
#    points(prop[1], prop[2], col="red",cex=2,lwd=2, pch = pchArray[i])
#}

deploid.prop = read.csv("deploid.prop.csv", header=F)
#for (i in 1:27){
for (i in 1:15){
    prop = deploid.prop[i,]
    if (length(which(prop>0.01))==3){
        prop = sort(prop, decreasing=T)
    }
#    print(prop)
    if ( i > 24 ){
        points(prop[1], prop[2], col="blue",cex=2,lwd=2, pch = pchArray[i-24])
    } else {
        points(prop[1], prop[2], col="blue",cex=2,lwd=2 , pch = i)
    }
}

legend("topright", legend = names2[1:12], pch = 1:12, cex=1.4, pt.lwd=2)
legend("bottomleft", legend = names2[13:24], pch = 13:24, cex=1.4, pt.lwd=2)
legend("bottom", legend = c(names3), pch = pchArray, cex=1.4, pt.lwd=2)

legend("top", legend = c("True proportion", "Proportion by pfmix", "Proportion by DEploid"), text.col = c("green", "red", "blue"), cex=1.4)

dev.off()
