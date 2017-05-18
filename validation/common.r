fun.divide.to.seg <- function(hapLength, by = 50){
#    return(floor(seq(1, hapLength, length.out =numSeg)))
    myseq = seq(1, hapLength, by = by)
    if ((hapLength-1)%%by != 0){
        myseq = c(myseq, hapLength)
    }
    return(myseq)
}

getIndex3 <- function (mixedSample, ref1, ref2, ref3, prop){
#1	1	2	2	3	3
#2	3	1	3	2	1
#3	2	3	1	1	2
#    prop = c(.25,.25,.5)
    prop = c(1,1,1)
    mySums = c ( sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
                 sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] ,
                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] )
#print(mySums)
#print("")
    case = which.min(mySums)
    if ( case == 1 ){
        return (c(1,2,3))
    } else if ( case == 2 ){
        return (c(1,3,2))
    } else if ( case == 3 ){
        return (c(2,1,3))
    } else if ( case == 4 ){
        return (c(2,3,1))
    } else if ( case == 5 ){
        return (c(3,2,1))
    } else if ( case == 6 ){
        return (c(3,1,2))
    }
}


getIndex2 <- function (mixedSample, ref1, ref2){
    sum1 = sum(mixedSample[,1] != ref1) + sum(mixedSample[,2] != ref2)
    sum2 = sum(mixedSample[,2] != ref1) + sum(mixedSample[,1] != ref2)
    if ( sum1 < sum2 ){
        return (c(1,2))
    } else {
        return (c(2,1))
    }
}


exportSites <- function ( CHROM, POS, tmpIndex, logFileName){
    for ( index in tmpIndex ){
        cat ( as.character(CHROM[index]), ",", POS[index], "\n", file = logFileName, append = T)
    }
}


fun.computeErrors2 <- function(hap, ref1, ref2){
    switchError = c(0, 0)
    mutError = c(0, 0)

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    if ( nhap == 2 ){

        index.of.seg = fun.divide.to.seg(haplength)
        countSwitch = matrix(NA, ncol = length(index.of.seg)-1, nrow = 2)

        for ( i in 1:(length(index.of.seg)-1) ){
            tmpIndex = c(index.of.seg[i]:(index.of.seg[i+1]-1))
            # k = 1
            if ( sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref2[tmpIndex]) ){
                hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
                countSwitch[1,i] = 1
            } else {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref2[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 7
                countSwitch[1,i] = 2
            }
#            } else if ( sum(hap[tmpIndex,1] != ref3[tmpIndex]) <= sum(hap[tmpIndex,1] != ref2[tmpIndex]) ) {
#                hap[tmpIndex,1][hap[tmpIndex,1] != ref3[tmpIndex]] = 2
##                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 8
#                countSwitch[1,i] = 3
#            }

            # k = 2
            if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) ){
                hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                countSwitch[2,i] = 2
            } else {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
#                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 8
                countSwitch[2,i] = 1
            }
#            } else if ( sum(hap[tmpIndex,2] != ref1[tmpIndex]) <= sum(hap[tmpIndex,2] != ref3[tmpIndex]) ) {
#                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
#                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
#                countSwitch[2,i] = 1
#            }
        }
        strain = cbind(c(NA,NA), countSwitch)
        strainNext = cbind(countSwitch, c(NA,NA))
        switchError = rowSums(strain != strainNext, na.rm=T)
        mutError = colSums(hap==2)
    }

    return ( list ( hap = hap,
                    switchError = switchError,
                    mutError = mutError ))

}


fun.plotHapWithProp <- function( hap, prop, fig.title, max.at ){

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, max.at)
    yrange = c(0, 1)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 4, cex.axis=2)
    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 3, cex.axis=2, xaxt="n")

    xleft = 0:(haplength-1)
    xright = xleft+1
#    ycum = as.numeric(c(0, prop[1], 1))
    ycum = as.numeric(c(0, cumsum(as.numeric(prop))))
#cat(ycum, "\n")

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }
    lab.at = c(1, 400, 800, 1200, 1600, 2000, 2369)
    axis(1, at=lab.at, labels=lab.at, las=1, lwd = 0, cex=2, cex.axis=2.4)

}

fun.computeErrors3 <- function(hap, ref1, ref2, ref3){
    switchError = c(0, 0, 0)
    mutError = c(0, 0, 0)

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    if ( nhap == 3 ){
#        hap[hap[,1] != ref1, 1 ] = 2
#        hap[hap[,2] != ref2, 2 ] = 2
#        hap[hap[,3] != ref3, 3 ] = 2

        index.of.seg = fun.divide.to.seg(haplength)
        countSwitch = matrix(NA, ncol = length(index.of.seg)-1, nrow = 3)

        for ( i in 1:(length(index.of.seg)-1) ){
            tmpIndex = c(index.of.seg[i]:(index.of.seg[i+1]-1))
            # k = 1
            if ( sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref2[tmpIndex]) & sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref3[tmpIndex]) ){
                hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
                hap[tmpIndex,1][hap[tmpIndex,1] == 1] = 5 # turn off the black color
                countSwitch[1,i] = 1
            } else if ( sum(hap[tmpIndex,1] != ref2[tmpIndex]) < sum(hap[tmpIndex,1] != ref3[tmpIndex]) ) {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref2[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 7
                hap[tmpIndex,1][hap[tmpIndex,1] == 1] = 7 # turn off the black color
                countSwitch[1,i] = 2
            } else if ( sum(hap[tmpIndex,1] != ref3[tmpIndex]) <= sum(hap[tmpIndex,1] != ref2[tmpIndex]) ) {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref3[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 1] = 0 # turn off the black color
                countSwitch[1,i] = 3
            }

            # k = 2
            if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) & sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref3[tmpIndex]) ){
                hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                hap[tmpIndex,2][hap[tmpIndex,2] == 1] = 7 # turn off the black color
                countSwitch[2,i] = 2
            } else if ( sum(hap[tmpIndex,2] != ref3[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) ) {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref3[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 1] = 0 # turn off the black color
                countSwitch[2,i] = 3
            } else if ( sum(hap[tmpIndex,2] != ref1[tmpIndex]) <= sum(hap[tmpIndex,2] != ref3[tmpIndex]) ) {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
                hap[tmpIndex,2][hap[tmpIndex,2] == 1] = 5 # turn off the black color
                countSwitch[2,i] = 1
            }

            # k = 3
            if ( sum(hap[tmpIndex,3] != ref3[tmpIndex]) < sum(hap[tmpIndex,3] != ref2[tmpIndex]) & sum(hap[tmpIndex,3] != ref3[tmpIndex]) < sum(hap[tmpIndex,3] != ref1[tmpIndex]) ){
                hap[tmpIndex,3][hap[tmpIndex,3] != ref3[tmpIndex]] = 2
                hap[tmpIndex,3][hap[tmpIndex,3] == 1] = 0 # turn off the black color
                countSwitch[3,i] = 3
            } else if ( sum(hap[tmpIndex,3] != ref2[tmpIndex]) < sum(hap[tmpIndex,3] != ref1[tmpIndex]) ) {
                hap[tmpIndex,3][hap[tmpIndex,3] != ref2[tmpIndex]] = 2
                hap[tmpIndex,3][hap[tmpIndex,3] == 0] = 7
                hap[tmpIndex,3][hap[tmpIndex,3] == 1] = 7 # turn off the black color
                countSwitch[3,i] = 2
            } else if ( sum(hap[tmpIndex,3] != ref1[tmpIndex]) <= sum(hap[tmpIndex,3] != ref2[tmpIndex]) ) {
                hap[tmpIndex,3][hap[tmpIndex,3] != ref1[tmpIndex]] = 2
                hap[tmpIndex,3][hap[tmpIndex,3] == 0] = 5
                hap[tmpIndex,3][hap[tmpIndex,3] == 1] = 5 # turn off the black color
                countSwitch[3,i] = 1
            }
        }
        strain = cbind(c(NA,NA,NA), countSwitch)
        strainNext = cbind(countSwitch, c(NA,NA,NA))
        switchError = rowSums(strain != strainNext, na.rm=T)
        mutError = colSums(hap==2)
    }

    return ( list ( hap = hap,
                    switchError = switchError,
                    mutError = mutError ))

}


measure.error.joe<-function(h.pair, h.pair.true, rel.cost.switch=2, do.plot=FALSE) {
    l <- ncol(h.pair);
    n.hap <- nrow(h.pair)
    possible.permn = combinat::permn(1:n.hap)
    n.permn = length(possible.permn)
    v<-rep(0, n.permn);
    vn<-v;

    tb<-array(0, c(n.permn, l));

    for ( j in 1:n.permn){
        v[j] = sum(h.pair[,1]!=h.pair.true[possible.permn[[j]],1]);
    }

    ee <- rep(0, n.permn)
    for (i in 2:l) {
        for ( j in 1:n.permn){
            ee[j] = sum(h.pair[,i]!=h.pair.true[possible.permn[[j]],i]);
            ones = rep(1, n.permn)
            ones[j] = 0
            tmp <- v + rel.cost.switch * ones
            vn[j] <- min(tmp) + ee[j]
            tb[j, i] <- which.min(tmp)
        }
        v<-vn;
    }

    #decode
    wm<-which.min(v);
    op<-array(0, l);
    n.s<-0;

    if (wm!=0){
        n.gt<-sum(h.pair[,l]!=h.pair.true[possible.permn[[wm]],l]);
    }


    op[l]<-wm;
    for (i in l:2) {
        wmp<-tb[wm,i];
        if (wmp!=wm) n.s<-n.s+1;

        if (wmp!=0){
            n.gt <- n.gt + sum(h.pair[,i-1] != h.pair.true[possible.permn[[wmp]],i-1]);
        }

        wm<-wmp;
        op[i-1]<-wm;
    }
	cat("\nDecoding gives:\nNo. switches:\t", n.s, "\nNo. GT errs:\t", n.gt, "\n\n");

	if (do.plot) {
		plot(0,0,type="n", xlab="Position", ylab="", yaxt="n", xlim=c(0,ncol(h.pair)), bty="n", ylim=c(-0.5,1.5));
		image(x=1:ncol(h.pair), z=t(rbind(h.pair, rep(0, ncol(h.pair)), h.pair.true)), col=c("white", "black"),
			add=T);
		del<-which(diff(op)!=0);
		if (length(del)>0) points(x=del, y=rep(0.5, length(del)), pch="|", col="red");

		w2<-which(op==2);
		for (i in w2) h.pair[,i]<-h.pair[2:1,i];
		wd1<-which(h.pair[1,] != h.pair.true[1,]);
		wd2<-which(h.pair[2,] != h.pair.true[2,]);
		if (length(wd1)>0) points(x=wd1, y=rep(1.2, length(wd1)), pch=25, col="blue", cex=0.5);
		if (length(wd2)>0) points(x=wd2, y=rep(1.2, length(wd2)), pch=25, col="green", cex=0.5);
	}

#	return(c(n.s, n.gt));
        return (list(switchError = n.s,
                 mutError = n.gt,
                 op = op) )
}

measure.error.joe.2<-function(h.pair, h.pair.true, rel.cost.switch=2, do.plot=FALSE) {
    hapAndError = measure.error.joe(h.pair, h.pair.true, rel.cost.switch, do.plot)

    n.hap = nrow(h.pair)
#    cat("n.hap = ", n.hap, "\n")
    switchError = rep(0, n.hap)
    mutError = rep(0, n.hap)
    possible.permn = combinat::permn(1:n.hap)

    l = ncol(h.pair.true)
    hap = array(0, c(nrow(h.pair), l));
    for ( i in 1:l ){
        hap[,i] = possible.permn[[hapAndError$op[i]]]
    }

    for ( j in 1:n.hap ){
        switchError[j] = sum((hap[j,-l] - hap[j, -1]) != 0)
    }

    for ( i in 1:l ){
        hap[h.pair[,i] != h.pair.true[hap[,i],i], i ] = 0
    }

    for ( j in 1:n.hap ){
        mutError[j] = sum(hap[j,] == 0)
    }

    return ( list ( hap = t(hap),
                    switchError = switchError,
                    mutError = mutError ))

}

