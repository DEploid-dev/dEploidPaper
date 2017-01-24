rm(list= ls())
source("common.r")



dataDir = "./"
panel = read.table(paste(dataDir,"labStrains.14.panel.txt",sep=""),header=T)

endAt = cumsum(table(panel[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

Ref1 = panel[,5] # HB3
Ref1Name = "HB3"
Ref2 = panel[,6] # 7G8
Ref2Name = "7G8"
Ref3 = panel[,4] # Dd2
Ref3Name = "Dd2"


sampleName = "PG0396-C"
cases = paste(sampleName, c(".14.noPanel",
                            ".14.asiaAfirca",
                            ".14.asiaAfirca_hb3",
                            ".14.asiaAfirca_hb3_7g8",
                            ".14.asiaAfirca_hb3_7g8_dd2"), sep="")
#cases = paste(sampleName, c(".14.noPanel",
#                            ".14.asiaAfirca",
#                            ".14.asiaAfirca_hb3",
#                            ".14.asiaAfirca_hb3_7g8",
#                            ".14.asiaAfirca_hb3_7g8_dd2",
#                            ".14.3d7_hb3_7g8_dd2"
#                            ,".14.asiaAfrica2"
#                            ), sep="")

#png(paste("differentPanelForSample.", sampleName, "Mid.png", sep=""), width = 4800, height = 2700, res=170)
#pdf(paste("differentPanelForSample.", sampleName, ".pdf", sep=""), width = 18, height = 13)
pdf(paste("Fig2.pdf", sep=""), width = 18, height = 13)
#png(paste("differentPanelForSample.", sampleName, "Hi.png", sep=""), width = 9600, height = 5400, res=350)
#png(paste("differentPanelForSample.", sampleName, "Hi2.png", sep=""), width = 9600, height = 7560, res=350)

par ( mfrow = c(length(cases),1), oma=c(6, 6, 2, 2), las=1)

panelNames = c("No LD,",
               "Panel I,",
               "Panel II,",
               "Panel III,",
               "Panel IV,")
#panelNames = c("No panel,",
#               "Panel I,",
#               "Panel II,",
#               "Panel III,",
#               "Panel IV,",
#               "Panel V,"
#               ,"Panel VI,"
#               )
#errors = c()
paneli=1
for ( prefix in cases ){
print(prefix)
    tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])

    hap = as.matrix(read.table(paste(prefix,".hap",sep=""), header=T)[,c(-1,-2)])

    colIndex = which(prop>0.01)

    prop.corrected = prop[colIndex]
    hap.corrected = hap[,colIndex,drop=FALSE]


    switchError = 0
    mutError = 0
    for ( chrom in 1:length(beginAt)){
        tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
        tmpProp = prop.corrected
        tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
        tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]
        tmpRef3 = Ref3[beginAt[chrom]:endAt[chrom]]

        if ( length(prop.corrected) == 3 ){
            rearranged.Index = getIndex3(hap.corrected, tmpRef1, tmpRef2, tmpRef3, tmpProp)
            tmpHap = tmpHap[,rearranged.Index]
            tmpProp = tmpProp = prop.corrected[rearranged.Index]
        }

        hapAndError = fun.computeErrors3( tmpHap, tmpRef1, tmpRef2, tmpRef3)
        cat(prefix, " ", hapAndError$switchError, " ", hapAndError$mutError, "\n")
#        tmpTitle = paste(prefix, rownames(table(panel[,1]))[chrom], hapAndError$switchError, "switch errors", hapAndError$mutError, "miss copy errors")
        tmpTitle = paste(panelNames[paneli]," total switch errors:", sum(hapAndError$switchError), ", total genotype errors:", sum(hapAndError$mutError))
#        tmpTitle = ""
        fun.plotHapWithProp (hapAndError$hap, tmpProp,
             tmpTitle,
             max(chromLength))
        switchError = switchError + hapAndError$switchError
        mutError = mutError + hapAndError$mutError
    }
    write.table(c(switchError, mutError), file = paste(prefix,".errorCount", sep=""), quote = F, row.names=F, col.names=F)
    paneli= paneli+1
}

mtext("Marker index", 1, 3, outer=TRUE, cex = 2)
mtext("Strain Proportion", 2, 3, outer=TRUE, cex = 2, las=0)

dev.off()

#tiff(paste("differentPanelForSample.", sampleName, "Hi.tif", sep=""), width = 9600, height = 5400, res=350)

#par ( mfrow = c(length(cases),1))

#panelNames = c("No panel,",
#               "Panel I,",
#               "Panel II,",
#               "Panel III,",
#               "Panel IV,")
##errors = c()
#paneli=1
#for ( prefix in cases ){
#print(prefix)
#    tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
#    prop = as.numeric(tmpProp[dim(tmpProp)[1],])

#    hap = as.matrix(read.table(paste(prefix,".hap",sep=""), header=T)[,c(-1,-2)])

#    colIndex = which(prop>0.01)

#    prop.corrected = prop[colIndex]
#    hap.corrected = hap[,colIndex,drop=FALSE]


#    switchError = 0
#    mutError = 0
#    for ( chrom in 1:length(beginAt)){
#        tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
#        tmpProp = prop.corrected
#        tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
#        tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]
#        tmpRef3 = Ref3[beginAt[chrom]:endAt[chrom]]

#        if ( length(prop.corrected) == 3 ){
#            rearranged.Index = getIndex3(hap.corrected, tmpRef1, tmpRef2, tmpRef3, tmpProp)
#            tmpHap = tmpHap[,rearranged.Index]
#            tmpProp = tmpProp = prop.corrected[rearranged.Index]
#        }

#        hapAndError = fun.computeErrors3( tmpHap, tmpRef1, tmpRef2, tmpRef3)
#        cat(prefix, " ", hapAndError$switchError, " ", hapAndError$mutError, "\n")
##        tmpTitle = paste(prefix, rownames(table(panel[,1]))[chrom], hapAndError$switchError, "switch errors", hapAndError$mutError, "miss copy errors")
#        tmpTitle = paste(panelNames[paneli]," total switch errors:", sum(hapAndError$switchError), ", total genotype errors:", sum(hapAndError$mutError))
##        tmpTitle = ""
#        fun.plotHapWithProp (hapAndError$hap, tmpProp,
#             tmpTitle,
#             max(chromLength))
#        switchError = switchError + hapAndError$switchError
#        mutError = mutError + hapAndError$mutError
#    }
#    write.table(c(switchError, mutError), file = paste(prefix,".errorCount", sep=""), quote = F, row.names=F, col.names=F)
#    paneli= paneli+1
#}

#dev.off()
