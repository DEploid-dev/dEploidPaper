rm(list=ls())
plaftab=read.table("labStrains.wg.PLAF.txt",header=T)
nonZeroIndex = which(plaftab$PLAF>0)

sampleName = "PG0412-C"
vcfName = paste(sampleName, ".wg.vcf", sep="")
skipNum = as.numeric(system(paste("cat ", vcfName, " | head -500 | grep \"##\" | wc -l"), T))
vcf  = read.table( vcfName, skip=skipNum, header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
names(vcf)[1] = "#CHROM"

#egVcfName = paste(sampleName, ".eg.vcf", sep="")
#system ( paste("grep \"##\"", vcfName, ">", egVcfName) )
#write.table(vcf[nonZeroIndex,], file = egVcfName, append = T, sep = "\t", quote = F, row.names = F)

asiaPanel = read.csv("asiaGroup1PanelMostDiverse10.wg.csv", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
#write.table (asiaPanel[nonZeroIndex,][testIndex,], file = "asia1.14.panel.txt", quote=F, sep="\t", row.names=F)

africaPanel = read.csv("africaGroup1PanelMostDiverse10.wg.csv", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
#write.table (africaPanel[nonZeroIndex,][testIndex,], file = "africa1.14.panel.txt", quote=F, sep="\t", row.names=F)


testIndex = which(plaftab[nonZeroIndex,]$CHROM == "Pf3D7_14_v3")
testVcfName = paste(sampleName, ".14.vcf", sep="")
system ( paste("grep \"##\"", vcfName, ">", testVcfName) )
write.table(vcf[nonZeroIndex,][testIndex,], file = testVcfName, append = T, sep = "\t", quote = F, row.names = F)


#write.table (panel[nonZeroIndex,], file = "labStrains.eg.panel.txt", quote=F, sep="\t", row.names=F)

write.table (plaftab[nonZeroIndex,][testIndex,], file = "labStrains.14.PLAF.txt", quote=F, sep="\t", row.names=F)

panel = read.table("labStrains.wg.panel.txt", header=T)
panel = panel[nonZeroIndex,][testIndex,]

asiaAfricaPanel = cbind( asiaPanel[nonZeroIndex,1:7][testIndex,], africaPanel[nonZeroIndex,3:7][testIndex,])
write.table(asiaAfricaPanel, file = "asiaAfrica.14.panel.txt", quote=F, sep="\t", row.names=F)


asiaAfricaPanel[["HB3gt.from.regression"]] = panel[["HB3gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3.14.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel[["sevenG8gt.from.regression"]] = panel[["sevenG8gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3_7g8.14.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel[["dd2gt.from.regression"]] = panel[["dd2gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3_7g8_dd2.14.panel.txt", quote=F, sep="\t", row.names=F)



asiaPanel = read.csv("asiaGroup1PanelMostDiverse10.wg.csv", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
#write.table (asiaPanel[nonZeroIndex,][testIndex,], file = "asia1.14.panel.txt", quote=F, sep="\t", row.names=F)

africaPanel = read.csv("africaGroup1PanelMostDiverse10.wg.csv", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)

testIndex = which(plaftab[nonZeroIndex,]$CHROM %in% c("Pf3D7_13_v3", "Pf3D7_14_v3"))
testVcfName = paste(sampleName, ".1314.vcf", sep="")
system ( paste("grep \"##\"", vcfName, ">", testVcfName) )
write.table(vcf[nonZeroIndex,][testIndex,], file = testVcfName, append = T, sep = "\t", quote = F, row.names = F)

write.table (plaftab[nonZeroIndex,][testIndex,], file = "labStrains.1314.PLAF.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel = cbind( asiaPanel[nonZeroIndex,1:7][testIndex,], africaPanel[nonZeroIndex,3:7][testIndex,])
write.table(asiaAfricaPanel, file = "asiaAfrica.1314.panel.txt", quote=F, sep="\t", row.names=F)

panel = read.table("labStrains.wg.panel.txt", header=T)
panel = panel[nonZeroIndex,][testIndex,]

asiaAfricaPanel[["HB3gt.from.regression"]] = panel[["HB3gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3.1314.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel[["sevenG8gt.from.regression"]] = panel[["sevenG8gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3_7g8.1314.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel[["dd2gt.from.regression"]] = panel[["dd2gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3_7g8_dd2.1314.panel.txt", quote=F, sep="\t", row.names=F)


asiaPanel = read.csv("asiaGroup1PanelMostDiverse10.wg.csv", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
#write.table (asiaPanel[nonZeroIndex,][testIndex,], file = "asia1.14.panel.txt", quote=F, sep="\t", row.names=F)

africaPanel = read.csv("africaGroup1PanelMostDiverse10.wg.csv", header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)

testIndex = which(plaftab[nonZeroIndex,]$CHROM %in% c("Pf3D7_12_v3", "Pf3D7_13_v3", "Pf3D7_14_v3"))
testVcfName = paste(sampleName, ".121314.vcf", sep="")
system ( paste("grep \"##\"", vcfName, ">", testVcfName) )
write.table(vcf[nonZeroIndex,][testIndex,], file = testVcfName, append = T, sep = "\t", quote = F, row.names = F)

write.table (plaftab[nonZeroIndex,][testIndex,], file = "labStrains.121314.PLAF.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel = cbind( asiaPanel[nonZeroIndex,1:7][testIndex,], africaPanel[nonZeroIndex,3:7][testIndex,])
write.table(asiaAfricaPanel, file = "asiaAfrica.121314.panel.txt", quote=F, sep="\t", row.names=F)

panel = read.table("labStrains.wg.panel.txt", header=T)
panel = panel[nonZeroIndex,][testIndex,]

asiaAfricaPanel[["HB3gt.from.regression"]] = panel[["HB3gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3.121314.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel[["sevenG8gt.from.regression"]] = panel[["sevenG8gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3_7g8.121314.panel.txt", quote=F, sep="\t", row.names=F)

asiaAfricaPanel[["dd2gt.from.regression"]] = panel[["dd2gt.from.regression"]]
write.table(asiaAfricaPanel, file = "asiaAfrica_hb3_7g8_dd2.121314.panel.txt", quote=F, sep="\t", row.names=F)


