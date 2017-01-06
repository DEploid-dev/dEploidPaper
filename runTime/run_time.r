rm ( list = ls () )
#if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2",repos="http://ftp.ussg.iu.edu/CRAN/")}
library(ggplot2)

#PG0412-C.121314.asiaAfirca_hb3_7g8_dd2.time  PG0412-C.1314.asiaAfirca_hb3.time
#PG0412-C.121314.asiaAfirca_hb3_7g8.time      PG0412-C.1314.asiaAfrica.time
#PG0412-C.121314.asiaAfirca_hb3.time          PG0412-C.14.asiaAfirca_hb3_7g8_dd2.time
#PG0412-C.121314.asiaAfrica.time              PG0412-C.14.asiaAfirca_hb3_7g8.time
#PG0412-C.1314.asiaAfirca_hb3_7g8_dd2.time    PG0412-C.14.asiaAfirca_hb3.time
#PG0412-C.1314.asiaAfirca_hb3_7g8.time        PG0412-C.14.asiaAfrica.time


#panels = paste("asiaAfrica", c("", "_hb3", "_hb3_7g8", "_hb3_7g8_dd2"), sep="")
panels = c("labStrains", "asiaAfrica", "asiaAfrica2")
chroms = c("14", "1314", "121314")
#     5992 PG0412-C.121314.vcf
#     4627 PG0412-C.1314.vcf
#     2591 PG0412-C.14.vcf
seqlength = as.character(c(2591, 4627, 5992))
#panelSize = as.character(c(10, 11, 12, 13))
panelSize = as.character(c(4, 10, 16))
prefix = "PG0412-C."

Time_data = data.frame( Panel = character(), Chrom = character(), PanelSize = character(), SeqLength = character(), Runtime = numeric(0), se = numeric(0) , stringsAsFactors=FALSE)
entry = 1

chromI = 1
for (chrom in chroms){
    paneli = 1
    for (panel in panels){
        current_job = c()
        fileName = paste(prefix, chrom, ".", panel, ".time", sep="")
        dummy_time = NaN
        print(fileName)
        if (class(try(read.table(fileName, skip=118, header=F, comment.char=""))) != "try-error"){
            dummy_time = read.table (fileName, skip=118, header=F, comment.char="")$V2[2]
            print(dummy_time)
        }
        current_job = c(current_job, dummy_time)
        Time_data[entry,]$Panel = panel
        Time_data[entry,]$Chrom   = chrom
        Time_data[entry,]$Runtime  = mean(current_job)
        Time_data[entry,]$PanelSize  = panelSize[paneli]
        Time_data[entry,]$SeqLength  = seqlength[chromI]
        print(Time_data[entry,])
        entry =  entry + 1
        paneli = paneli + 1
    }
    chromI = chromI + 1
}

Time_data$Panel <- factor(Time_data$Panel)
Time_data$Chrom <- factor(Time_data$Chrom)

Time_data$PanelSize <- factor(Time_data$PanelSize, levels = panelSize)
Time_data$SeqLength <- factor(Time_data$SeqLength, levels = seqlength)

print(Time_data)

library(scales)     # Need the scales package

myBar = ggplot(Time_data, aes(x=PanelSize, y=Runtime, fill=SeqLength)) +
    ylab("Time / seconds")+
    geom_bar(position=position_dodge(), stat="identity") +
    scale_y_continuous(limits=c(40,1200),oob = rescale_none)
png("runTime.png", width = 800, height = 800)
myBar+ theme(text = element_text(size=20))
dev.off()

png("runTimeHi.png", width = 4000, height = 4000, res=350)
myBar+ theme(text = element_text(size=20))
dev.off()

tiff("runTimeHi.tif", width = 4000, height = 4000, res=350)
myBar+ theme(text = element_text(size=20))
dev.off()


#    myBar+ geom_line(aes(x=as.numeric(PanelSize), y=Runtime), colour="blue")

#plot(as.numeric(as.character(Time_data$SeqLength))[c(1,5,9)], Time_data$Runtime[c(1,5,9)],type="l")
#lines(as.numeric(as.character(Time_data$SeqLength))[c(1,9)], Time_data$Runtime[c(1,9)],type="l")

#index = c(2,5,10)
#plot(as.numeric(as.character(Time_data$SeqLength))[index], Time_data$Runtime[index],type="l")
#lines(as.numeric(as.character(Time_data$SeqLength))[c(2,10)], Time_data$Runtime[c(2,10)],type="l")

#index = c(2,5,10)+2
#plot(as.numeric(as.character(Time_data$SeqLength))[index], Time_data$Runtime[index],type="l")
#lines(as.numeric(as.character(Time_data$SeqLength))[c(2,10)+2], Time_data$Runtime[c(2,10)+2],type="l")
