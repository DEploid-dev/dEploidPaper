rm(list=ls())
plaftab=read.table("labStrains.wg.PLAF.txt",header=T)
nonZeroIndex = which(plaftab$PLAF>0)
write.table (plaftab[nonZeroIndex,], file = "labStrains.eg.PLAF.txt", quote=F, sep="\t", row.names=F)

sampleName = "PG0394-C"
#sampleName = "PG0389-C"
#sampleName = "PG0390-C"
#sampleName = "PG0402-C"
vcfName = paste(sampleName, ".wg.vcf", sep="")
skipNum = as.numeric(system(paste("cat ", vcfName, " | head -500 | grep \"##\" | wc -l"), T))
vcf  = read.table( vcfName, skip=skipNum, header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
names(vcf)[1] = "#CHROM"

egVcfName = paste(sampleName, ".eg.vcf", sep="")
system ( paste("grep \"##\"", vcfName, ">", egVcfName) )
write.table(vcf[nonZeroIndex,], file = egVcfName, append = T, sep = "\t", quote = F, row.names = F)

testIndex = which(plaftab[nonZeroIndex,]$CHROM == "Pf3D7_14_v3")
testVcfName = paste(sampleName, ".14.vcf", sep="")
system ( paste("grep \"##\"", vcfName, ">", testVcfName) )
write.table(vcf[nonZeroIndex,][testIndex,], file = testVcfName, append = T, sep = "\t", quote = F, row.names = F)

panel = read.table("labStrains.wg.panel.txt", header=T)
write.table (panel[nonZeroIndex,], file = "labStrains.eg.panel.txt", quote=F, sep="\t", row.names=F)

write.table (plaftab[nonZeroIndex,][testIndex,], file = "labStrains.14.PLAF.txt", quote=F, sep="\t", row.names=F)
write.table (panel[nonZeroIndex,][testIndex,], file = "labStrains.14.panel.txt", quote=F, sep="\t", row.names=F)

asiaPanel = read.csv("asiaGroup1PanelMostDiverse10.wg.csv", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
#write.table (asiaPanel[nonZeroIndex,][testIndex,], file = "asia1.14.panel.txt", quote=F, sep="\t", row.names=F)

africaPanel = read.csv("africaGroup1PanelMostDiverse10.wg.csv", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
#write.table (africaPanel[nonZeroIndex,][testIndex,], file = "africa1.14.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel = cbind( asiaPanel[nonZeroIndex,1:7][testIndex,], africaPanel[nonZeroIndex,3:7][testIndex,])
write.table(asiaAfricaPanel, file = "asiaAfirca.14.panel.txt", quote=F, sep="\t", row.names=F)

panel = read.table(paste("labStrains.14.panel.txt",sep=""),header=T)

asiaAfricaPanel[["HB3gt.from.regression"]] = panel[["HB3gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfirca_hb3.14.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel[["sevenG8gt.from.regression"]] = panel[["sevenG8gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfirca_hb3_7g8.14.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel[["dd2gt.from.regression"]] = panel[["dd2gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfirca_hb3_7g8_dd2.14.panel.txt", quote=F, sep="\t", row.names=F)


asiaPlus = read.table("asia1.14.panel.txt", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
asiaPlus[["PH1000-C"]] <- NULL
asiaPlus[["PH0848-C"]] <- NULL

panel = read.table(paste(dataDir,"labStrains.14.panel.txt",sep=""),header=T)
asiaPlus[["threeD7"]] = panel[["threeD7"]]
asiaPlus[["HB3gt.from.regression"]] = panel[["HB3gt.from.regression"]]
write.table (asiaPlus, file = "asiaPlus.14.panel.txt", quote=F, sep="\t", row.names=F)

asiaPlus[["sevenG8gt.from.regression"]] = panel[["sevenG8gt.from.regression"]]
asiaPlus[["HB3gt.from.regression"]] = panel[["HB3gt.from.regression"]]
write.table (asiaPlus, file = "asiaPlus2.14.panel.txt", quote=F, sep="\t", row.names=F)
