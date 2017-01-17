rm(list=ls())
# DESCRIPTION:
#
# USAGE:
#    R --slave "--args sampleName 2 seed" < run_pfmix.r

args = (commandArgs(TRUE))

sample.name = args[1] #"PG0389-C"
k = as.numeric(args[2])
seed = as.numeric(args[3])

#sample.name = "PG0396-C"
#k = 3
#seed = 2

set.seed(seed)

source("mixture.funcs.r")

library(VGAM)
non = read.table(paste("altAndRef/", sample.name, ".eg.alt", sep=""), header=T, comment.char = "")[,3]
ref = read.table(paste("altAndRef/", sample.name, ".eg.ref", sep=""), header=T, comment.char = "")[,3]
allele.freq = read.table("labStrains.eg.PLAF.txt", header=T)[,3]

sample.data <-  cbind(non, non+ref)
#num.iter <- 4000;
num.iter <- 500;
#num.iter <- 100;
thin <- 10

model.out <- run.mcmc(num.test = k,allele.freq,sample.data,num.iter,thin, 1)
#current = model.out[[1]][[1]]
#print(current$p)
#prop = sort(current$p[2:(2+current$num.comp)])[1:current$num.comp]
propChain = c()
numSample = 0
for ( i in 1:(num.iter/thin) ){
    numSample = numSample+1
    propChain = rbind(propChain, sort(model.out[[i]][[1]]$p[2:(1+model.out[[i]][[1]]$num.comp)], decreasing=F))
}

prop = print(colMeans(propChain[(numSample-10):numSample,]))
#write.table(prop, file = paste(sample.name, ".pfmix.",seed,".prop", sep = ""), sep = "\t", row.names = F, col.names = F)
write.table(prop, file = paste(sample.name, ".pfmix",".prop", sep = ""), sep = "\t", row.names = F, col.names = F)
png(paste(sample.name, ".png", sep=""))
barplot(t(propChain))
dev.off()
