rm(list=ls())

fun.divide.to.seg <- function(hapLength, by = 3001){
#    return(floor(seq(1, hapLength, length.out =numSeg)))
    myseq = seq(1, hapLength, by = by)
    if ((hapLength-1)%%by != 0){
        myseq = c(myseq, hapLength)
    }
    return(myseq)
}

#getIndex3 <- function (mixedSample, ref1, ref2, ref3, prop){
##1	1	2	2	3	3
##2	3	1	3	2	1
##3	2	3	1	1	2
##    prop = c(.25,.25,.5)
#    prop = c(1,1,1)
#    mySums = c ( sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
#                 sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] ,
#                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
#                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
#                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
#                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] )
##print(mySums)
##print("")
#    case = which.min(mySums)
#    if ( case == 1 ){
#        return (c(1,2,3))
#    } else if ( case == 2 ){
#        return (c(1,3,2))
#    } else if ( case == 3 ){
#        return (c(2,1,3))
#    } else if ( case == 4 ){
#        return (c(2,3,1))
#    } else if ( case == 5 ){
#        return (c(3,2,1))
#    } else if ( case == 6 ){
#        return (c(3,1,2))
#    }
#}


getIndex2 <- function (mixedSample, ref1, ref2){
    sum1 = sum(mixedSample[,1] != ref1) + sum(mixedSample[,2] != ref2)
    sum2 = sum(mixedSample[,2] != ref1) + sum(mixedSample[,1] != ref2)
    if ( sum1 < sum2 ){
        return (c(1,2))
    } else {
        return (c(2,1))
    }
}


#exportSites <- function ( CHROM, POS, tmpIndex, logFileName){
#    for ( index in tmpIndex ){
#        cat ( as.character(CHROM[index]), ",", POS[index], "\n", file = logFileName, append = T)
#    }
#}


fun.calc.f.samp<-function( h, prop, tol=1e-6 ) {
  tmp<-h %*% prop;
  return (tmp);
}


fun.computeErrors2 <- function(hap, ref1, ref2, tmpPathEvent){
    switchError = c(0, 0)
    mutError = c(0, 0)

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    local.n.loci = length(ref1)

#    index.of.seg = fun.divide.to.seg(haplength)
    countSwitch = matrix(NA, ncol = local.n.loci, nrow = 2)

    for ( tmpIndex in 1:(local.n.loci) ){
#        tmpIndex = c(index.of.seg[i]:index.of.seg[i+1])
        # k = 1
        if ( tmpPathEvent[tmpIndex] == 3 ){
            hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
            hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
            countSwitch[1,tmpIndex] = 1

            hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
            hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
            countSwitch[2,tmpIndex] = 2


        } else { # tmpPathEvent[tmpIndex] == 2
            hap[tmpIndex,1][hap[tmpIndex,1] != ref2[tmpIndex]] = 2
            hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 7
            countSwitch[1,tmpIndex] = 2

            hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
            hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
            countSwitch[2,tmpIndex] = 1

        }
    }
    strain = cbind(c(NA,NA), countSwitch)
    strainNext = cbind(countSwitch, c(NA,NA))
#    print(strain)
    switchError = rowSums(strain != strainNext, na.rm=T)
    print(switchError)
    mutError = colSums(hap==2)

    return ( list ( hap = hap,
                    switchError = switchError,
                    mutError = mutError ))

}


fun.plotHapWithProp <- function( hap, prop, fig.title, max.at ){

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, max.at)
    yrange = c(0, 1)
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 2, cex.axis=2)

    xleft = 0:(haplength-1)
    xright = xleft+1
    ycum = c(0, cumsum(prop))

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }
}


fun.sample.two.hap.LS.emiss <- function ( logemiss, miss.copy.rate = 0.01 ){
#    emiss <- exp(logemiss);

#    lk.00 = emiss[,1]
#    lk.10 = emiss[,2]
#    lk.01 = emiss[,3]
#    lk.11 = emiss[,4]

#    emiss = cbind( lk.00 * (1-miss.copy.rate)^2 + (lk.10 + lk.01) * (miss.copy.rate)*(1-miss.copy.rate) + lk.11 * (miss.copy.rate)^2,
#                   lk.10 * (1-miss.copy.rate)^2 + (lk.00 + lk.11) * (miss.copy.rate)*(1-miss.copy.rate) + lk.01 * (miss.copy.rate)^2,
#                   lk.01 * (1-miss.copy.rate)^2 + (lk.00 + lk.11) * (miss.copy.rate)*(1-miss.copy.rate) + lk.10 * (miss.copy.rate)^2,
#                   lk.11 * (1-miss.copy.rate)^2 + (lk.10 + lk.01) * (miss.copy.rate)*(1-miss.copy.rate) + lk.00 * (miss.copy.rate)^2)

    # omu stands for 1 minus u, 1 minus miss.copy.rate
    llk.00 = logemiss[,1]
    llk.10 = logemiss[,2]
    llk.01 = logemiss[,3]
    llk.11 = logemiss[,4]

    log.omu = log(1-miss.copy.rate)
    log.u = log(miss.copy.rate)

    tmp.00.1 = llk.00 + log.omu + log.omu
    tmp.00.2 = llk.10 + log.omu + log.u
    tmp.00.3 = llk.01 + log.omu + log.u
    tmp.00.4 = llk.11 + log.u + log.u

    tmp.10.1 = llk.10 + log.omu + log.omu
    tmp.10.2 = llk.00 + log.omu + log.u
    tmp.10.3 = llk.11 + log.omu + log.u
    tmp.10.4 = llk.01 + log.u + log.u

    tmp.01.1 = llk.01 + log.omu + log.omu
    tmp.01.2 = llk.00 + log.omu + log.u
    tmp.01.3 = llk.11 + log.omu + log.u
    tmp.01.4 = llk.10 + log.u + log.u

    tmp.11.1 = llk.11 + log.omu + log.omu
    tmp.11.2 = llk.10 + log.omu + log.u
    tmp.11.3 = llk.01 + log.omu + log.u
    tmp.11.4 = llk.00 + log.u + log.u

    tmp.max = apply(cbind(tmp.00.1, tmp.00.2, tmp.00.3, tmp.00.4,
                          tmp.10.1, tmp.10.2, tmp.10.3, tmp.10.4,
                          tmp.01.1, tmp.01.2, tmp.01.3, tmp.01.4,
                          tmp.11.1, tmp.11.2, tmp.11.3, tmp.11.4), 1, max)
    emiss = cbind ( exp( tmp.00.1-tmp.max ) + exp( tmp.00.2-tmp.max ) + exp( tmp.00.3-tmp.max ) + exp( tmp.00.4-tmp.max ),
                    exp( tmp.10.1-tmp.max ) + exp( tmp.10.2-tmp.max ) + exp( tmp.10.3-tmp.max ) + exp( tmp.10.4-tmp.max ),
                    exp( tmp.01.1-tmp.max ) + exp( tmp.01.2-tmp.max ) + exp( tmp.01.3-tmp.max ) + exp( tmp.01.4-tmp.max ),
                    exp( tmp.11.1-tmp.max ) + exp( tmp.11.2-tmp.max ) + exp( tmp.11.3-tmp.max ) + exp( tmp.11.4-tmp.max ))

#    emiss<-emiss/apply(emiss, 1, sum); # sum of the row

    return (emiss)
}


fun.llk<-function(ref, alt, f.samp, err=0.01, fac=100) {
  cov.ref = ref
  cov.alt = alt

    f.samp<-f.samp+err*(1-2*f.samp);
    llk<-lbeta(cov.alt+f.samp*fac, cov.ref+(1-f.samp)*fac)-lbeta(f.samp*fac,(1-f.samp)*fac);
#    llk<-lgamma(fac*f.samp+cov.alt)+lgamma(fac*(1-f.samp)+cov.ref)-lgamma(fac*f.samp)-lgamma(fac*(1-f.samp));
    return(llk);
}

fun.computeRecProb <- function(coordinates, num_markers_per_chromosome, average_centimorgan_distance = 15000, Ne = 10) {
  # Compute genetic distances in Morgans.
  genetic_distances <- double(length(coordinates)-1)
  average_morgan_distance <- average_centimorgan_distance*100
  genetic_distances <- c((coordinates[-1] - coordinates[-length(coordinates)])/average_morgan_distance, Inf)

  # Set the amount of recombination between chromosomes to infinite.
  genetic_distances[cumsum(num_markers_per_chromosome)] <- Inf

  # Compute the population scaled recombination rate (i.e. 2*Ne*r)
  rho <- sapply(genetic_distances, function(r){ return(2*Ne*r)})

  # Finally, compute the probability of recombination within the segment separating each pair of SNPs.
  prob_recombination <- (1-(exp(-rho)))

#  if (use.const.recomb.prob){
#      prob_recombination[] = 0.01
#      prob_recombination[cumsum(num_markers_per_chromosome)] = 1
#  }
#  write.table(prob_recombination, file = "recomb.prob", sep="\t", row.names = FALSE, col.names = FALSE)
  return(prob_recombination)
}

fun.normalize.bymax <- function( f ){
    return (f/max(f))
}



fun.sample.two.hap.LS.cal.fwd <- function ( ref.panel, emiss, recomb.probs ){
  fun.buildEmission <-function(ref.panel.at.j, emiss.at.j){
      #> ref.panel.at.j = c(1.1,2.2,3.3,4.4,5.5)
      #> obs1 = matrix(ref.panel.at.j[rowIndex], nrow = n.panel, byrow = F )
      #>       obs2 = matrix(ref.panel.at.j[colIndex], nrow = n.panel, byrow = F )
      #> obs1
      #      [,1] [,2] [,3] [,4] [,5]
      #[1,]  1.1  1.1  1.1  1.1  1.1
      #[2,]  2.2  2.2  2.2  2.2  2.2
      #[3,]  3.3  3.3  3.3  3.3  3.3
      #[4,]  4.4  4.4  4.4  4.4  4.4
      #[5,]  5.5  5.5  5.5  5.5  5.5
      #> obs2
      #      [,1] [,2] [,3] [,4] [,5]
      #[1,]  1.1  2.2  3.3  4.4  5.5
      #[2,]  1.1  2.2  3.3  4.4  5.5
      #[3,]  1.1  2.2  3.3  4.4  5.5
      #[4,]  1.1  2.2  3.3  4.4  5.5
      #[5,]  1.1  2.2  3.3  4.4  5.5
      obs1 = matrix(ref.panel.at.j[rowIndex], nrow = n.panel, byrow = F )
      obs2 = matrix(ref.panel.at.j[colIndex], nrow = n.panel, byrow = F )
      #> obs1 = c(0,1,0,1)
      #> obs2 = c(0,0,1,1)
      #> obs1*1 + obs2*2+1
      #[1] 1 2 3 4
      obs = obs1*1 + obs2*2+1   # should be this ...
      # obs is a n.panel by n.panel matrix,
      # so is emiss.at.j.mat
      emiss.at.j.mat = matrix(emiss.at.j[obs], nrow = n.panel, byrow = F )
      return (emiss.at.j.mat)
  }

    n.panel<-nrow(ref.panel);
    # Let row index denote the first reference sequence,
    # and column index denote the second reference sequence
    #  > n.panel = 5
    #  > rowIndex = matrix(rep(1:n.panel,n.panel), nrow = n.panel, ncol = n.panel, byrow=F)
    #  > colIndex = matrix(rep(1:n.panel,n.panel), nrow = n.panel, ncol = n.panel, byrow=T)
    #  > rowIndex
    #       [,1] [,2] [,3] [,4] [,5]
    #  [1,]    1    1    1    1    1
    #  [2,]    2    2    2    2    2
    #  [3,]    3    3    3    3    3
    #  [4,]    4    4    4    4    4
    #  [5,]    5    5    5    5    5
    #  > colIndex
    #       [,1] [,2] [,3] [,4] [,5]
    #  [1,]    1    2    3    4    5
    #  [2,]    1    2    3    4    5
    #  [3,]    1    2    3    4    5
    #  [4,]    1    2    3    4    5
    #  [5,]    1    2    3    4    5
    rowIndex = matrix(rep(1:n.panel,n.panel), nrow = n.panel, ncol = n.panel, byrow=F)
    colIndex = matrix(rep(1:n.panel,n.panel), nrow = n.panel, ncol = n.panel, byrow=T)
    n.loci<-ncol(ref.panel);

    fwdList = list()
    fwdList[[1]] = fun.buildEmission(ref.panel[,1], emiss[1,])
    diag( fwdList[[1]] ) <- 0
    fwdList[[1]] = fun.normalize.bymax(fwdList[[1]])
    for (j in 2:n.loci) {
        marginalOfRows = matrix(rep(rowSums(fwdList[[j-1]]), n.panel), nrow = n.panel, ncol = n.panel, byrow = T )
        marginalOfCols = matrix(rep(colSums(fwdList[[j-1]]), n.panel), nrow = n.panel, ncol = n.panel, byrow = F )

        # Compute transitions based on recombination map.
        p.rec <- recomb.probs[j-1]
        p.rec.each.hap <- p.rec/n.panel
        p.no.rec <- (1-p.rec)

        rec.rec = p.rec.each.hap * p.rec.each.hap
        rec.norec = p.rec.each.hap * p.no.rec
        norec.norec = p.no.rec * p.no.rec

        fwdList[[j]] =  ( rec.rec * sum(fwdList[[j-1]]) +
                           norec.norec * fwdList[[j-1]] +
                           rec.norec * marginalOfRows * marginalOfCols ) *
                         fun.buildEmission( ref.panel[,j], emiss[j,] )
        diag( fwdList[[j]] ) <- 0
        fwdList[[j]] = fun.normalize.bymax(fwdList[[j]])

    }
    return (fwdList)
}


asia1 = read.table("asia1.14.panel.txt", header=T, check.names=F)
endAt = cumsum(table(asia1[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

coordinates <- as.integer(asia1$POS)
num.markers.chrm = table(asia1$CHROM)
# Compute the recombination probabilities.
recomb.probs <- fun.computeRecProb(coordinates, num.markers.chrm)


Ref1 = asia1[["PH0064-C"]]
Ref1Name = "PH0064-C"
Ref2 = asia1[["PH0193-C"]]
Ref2Name = "PH0193-C"

#Ref2 = asia1[["PH0064-C"]]
#Ref2Name = "PH0064-C"
#Ref1 = asia1[["PH0193-C"]]
#Ref1Name = "PH0193-C"


#if ( currentDir %in% c("PG0389-C", "PG0390-C", "PG0391-C", "PG0392-C", "PG0393-C", "PG0394-C") ){
#}

#cases = c("PG0412-C.14.noPanel",
#          "PG0412-C.14.asia",
#          "PG0412-C.14.asiaPlus",
#          "PG0412-C.14.labPanel")

#suffix = ".noError"
suffix = ".withError"
#sampleName = "25v75"
#sampleName = "75v25"
sampleName = "55v45"
#cases = paste(sampleName, "panel", c(1,2,3), sep="")
#cases = paste(sampleName, "panel", c(1,2,3),".nodis", sep="")
cases = c(paste(sampleName, "panel", c(1,2,3), suffix, sep=""),
          paste(sampleName, "noPanel", suffix, sep = ""))


ref = read.table(paste("mixedFieldSamplePH0063-PH0193",suffix,".",sampleName,".ref", sep=""), header=T)$REF
alt = read.table(paste("mixedFieldSamplePH0063-PH0193",suffix,".",sampleName,".alt", sep=""), header=T)$ALT




png(paste("differentPanelForSample.", sampleName, suffix, ".LS.png", sep=""), width = 1920, height = 1080)
par ( mfrow = c(length(cases),1))

for ( prefix in cases ){

    tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])

    hap = as.matrix(read.table(paste(prefix,".hap",sep=""), header=T)[,c(-1,-2)])

#    colIndex = which(prop>0.01)
colIndex = getIndex2(hap, Ref1, Ref2)
    prop.corrected = prop[colIndex]
    hap.corrected = hap[,colIndex,drop=FALSE]

expected.WSAF = fun.calc.f.samp(hap.corrected, prop.corrected)
llk = fun.llk (ref,alt,expected.WSAF)
  ws=c(1,2)
  expected.WSAF.00 <- expected.WSAF-(prop.corrected[ws[1]]*hap.corrected[,ws[1]] + prop.corrected[ws[2]]*hap.corrected[,ws[2]]);
  expected.WSAF.10 <- expected.WSAF.00 + prop.corrected[ws[1]];
  expected.WSAF.01 <- expected.WSAF.00 + prop.corrected[ws[2]];
  expected.WSAF.11 <- expected.WSAF.00 + prop.corrected[ws[1]] + prop.corrected[ws[2]];

  llk.00<-fun.llk(ref,alt, expected.WSAF.00)
  llk.10<-fun.llk(ref,alt, expected.WSAF.10)
  llk.01<-fun.llk(ref,alt, expected.WSAF.01)
  llk.11<-fun.llk(ref,alt, expected.WSAF.11)

llks<-cbind(llk.00, llk.10, llk.01, llk.11);

emiss = fun.sample.two.hap.LS.emiss ( llks )

ref.panel = rbind(Ref1, Ref2)

fwd.p = fun.sample.two.hap.LS.cal.fwd ( ref.panel, emiss, recomb.probs )

pathEvent = c()
for ( i in 1:2425 ){
    pathEvent = c(pathEvent, which.max(fwd.p[[i]]))
}
    switchError = 0
    mutError = 0
    for ( chrom in 1:length(beginAt)){
        tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
        tmpProp = prop.corrected
        tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
        tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]

#        if ( length(prop.corrected) == 2 ){
#            rearranged.Index = getIndex2(hap.corrected, Ref1, Ref2)
#            tmpHap = tmpHap[,rearranged.Index,drop=FALSE]
#            tmpProp = tmpProp = prop.corrected[rearranged.Index]
#        }

        hapAndError = fun.computeErrors2( tmpHap, tmpRef1, tmpRef2, pathEvent[beginAt[chrom]:endAt[chrom]])

        tmpTitle = paste(prefix, rownames(table(asia1[,1]))[chrom], hapAndError$switchError, "switch errors", hapAndError$mutError, "miss copy errors")

        fun.plotHapWithProp (hapAndError$hap, tmpProp,
             tmpTitle,
             max(chromLength))
        switchError = switchError + hapAndError$switchError
        mutError = mutError + hapAndError$mutError
    }
    write.table(c(switchError, mutError), file = paste(prefix,".errorCount", sep=""), quote = F, row.names=F, col.names=F)
}

dev.off()
