rm(list=ls())
pdf("compareToJack.pdf", width=16, height=16)
par(mfrow = c(2,2))

names2 = read.table("labSampleNames2Strains", stringsAsFactors=F, header=F)$V1
true.prop = read.csv("true.prop.csv", header=F)
deploid.prop = read.csv("deploid.prop.csv", header=F)
markerArray = c(0:4, 6, 19:15)
plotGroup <- function (group){
    groupIndex = 1:length(group)
    for (gi in groupIndex){
        i = group[gi]
        prop = true.prop[i,]
        if (length(which(prop>0.01))==3){
            prop = sort(prop, decreasing=T)
        }
        points(prop[1], prop[2], col="blue",cex=4,lwd=2, pch=markerArray[gi])
    }


    for (gi in groupIndex){
        i = group[gi]
        prop = read.table(paste(names2[i],".pfmix.prop", sep=""),head=F)$V1
        if ( i %in% c(1,2,3,7:16) ){
            prop = sort(prop, decreasing=T)
        }
#        print(prop)
        points(prop[1], prop[2], col="green",cex=3,lwd=2, pch = markerArray[gi])
    }


    for (gi in groupIndex){
        i = group[gi]
        prop = deploid.prop[i,]
        if (length(which(prop>0.01))==3){
            prop = sort(prop, decreasing=T)
        }

        points(prop[1], prop[2], col="red",cex=2, lwd=2 , pch = markerArray[gi])
    }

    legend("topright", legend = names2[group], pch = markerArray[groupIndex], cex=1.4, pt.lwd=2)
    legend("top", legend = c("True proportion", "Proportion by pfmix", "Proportion by DEploid"), text.col = c("blue", "green", "red"), cex=1.4)
}
###################################################
plot(c(-0.01,.21),c(1.01,.79), xlab = "Proportion component 1", ylab = "Proportion component 2", type="n")#\\, log="xy"
group1 = c(5,6,21:24)
plotGroup(group1)

###################################################
plot(c(0.16,.82),c(0.84,0.18), xlab = "Proportion component 1", ylab = "Proportion component 2", type="n")#\\, log="xy"
group2 = c(3,4,12:20)
plotGroup(group2)

###################################################
plot(c(0.79,1.01),c(0.21, -0.01), xlab = "Proportion component 1", ylab = "Proportion component 2", type="n")#\\, log="xy"
group3 = c(1,2,7:11)
plotGroup(group3)

###################################################
plot(c(0.3,.75),c(.4,0.1), xlab = "Proportion component 1", ylab = "Proportion component 2", type="n")#\\, log="xy"
group4 = c(25:27)

groupIndex = 1:length(group4)
for (gi in groupIndex){
    i = group4[gi]
    prop = true.prop[i,]
    if (length(which(prop>0.01))==3){
        prop = sort(prop, decreasing=T)
    }
#        points(prop[1], prop[2], col="blue",cex=3,lwd=2, pch=markerArray[gi])
    points(prop[1], prop[2], col="blue",cex=4, lwd =2, pch = markerArray[gi])
}


for (gi in groupIndex){
    i = group4[gi]
    prop = deploid.prop[i,]
    if (length(which(prop>0.01))==3){
        prop = sort(prop, decreasing=T)
    }

    points(prop[1], prop[2], col="red",cex=3, lwd=2 , pch = markerArray[gi])
}

names3 = read.table("labSampleNames3Strains", stringsAsFactors=F, header=F)$V1
for (i in 1:3){
    prop = sort(read.table(paste(names3[i],".pfmix.prop", sep=""),head=F)$V1, decreasing=T)
    print(prop)
    points(prop[1], prop[2], col="green",cex=2,lwd=2, pch = markerArray[i])
}


legend("topright", legend = names3, pch = markerArray, cex=1.4, pt.lwd=2)
legend("top", legend = c("True proportion", "Proportion by pfmix", "Proportion by DEploid"), text.col = c("blue", "green", "red"), cex=1.4)

#legend("bottomleft", legend = names2[13:24], pch = 13:24, cex=1.4, pt.lwd=2)
#

#legend("top", legend = c("True proportion", "Proportion by pfmix", "Proportion by DEploid"), text.col = c("green", "red", "blue"), cex=1.4)




dev.off()
