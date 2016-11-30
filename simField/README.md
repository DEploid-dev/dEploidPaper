simulate field mixture of two asia samples PH0064-C/PH0193-C, mixing with proportions (25/75) and (75/25), chromosome 14, both strains are presented in panel asia1.14.panel.txt, but only PH0064-C is presented in the asia1.14.panel.txt

Deconvolute using panels asia1.14.panel.txt  asiaAfirca.14.panel.txt  labStrains.14.panel.txt, and assess the quality

Use asiaGroup1 plaf, first trim plaf to chrom 14.
```R
rm(list= ls())
plaf = read.table("asiaGroup1_PLAF.txt",header=T)
plafIndex = paste(plaf$CHROM, plaf$POS, sep="-")
asia1 = read.table("asia1.14.panel.txt", header=T, check.names=F)
asia2 = asia1
asia2[["PH0064-C"]] = NULL
asia2[["PH0193-C"]] = NULL
write.table(asia2, "asia1-2.14.panel.txt", row.names = F, quote=F, col.names=T, sep="\t")
asia1Index = paste(asia1$CHROM, asia1$POS, sep="-")
keep.at = which(plafIndex %in% asia1Index)
write.table(plaf[keep.at,], file = "asiaGroup1_PLAF.14.txt", row.names = F, quote=F, col.names=T, sep="\t")

PH0064ref = read.table("PH0064-C_ref.txt", header=T, check.names=F)
PH0064alt = read.table("PH0064-C_alt.txt", header=T, check.names=F)
tmpTotalCoverage = PH0064ref[["PH0064-C"]][keep.at] +  PH0064alt[["PH0064-C"]][keep.at]

PH0064 = asia1[["PH0064-C"]]
PH0193 = asia1[["PH0193-C"]]

set.seed(1)
err<-0.01;
if ( err == 0 ){
    suffix = "noError"
} else {
    suffix = "withError"
}
n.loci = length(keep.at)

WSAF1 = PH0064 * 0.25 + PH0193 *0.75
includeErrorWSAF1 = WSAF1*(1-err)+(1-WSAF1)*err
altCount1 = rbinom(n.loci, tmpTotalCoverage, includeErrorWSAF1)
refCount1 = tmpTotalCoverage - altCount1

write.table(data.frame(CHROM = asia1$CHROM, POS = asia1$POS, REF = refCount1),
    file = paste("mixedFieldSamplePH0063-PH0193.", suffix, ".25v75.ref", sep = ""), row.names = F, quote=F, col.names=T, sep="\t")

write.table(data.frame(CHROM = asia1$CHROM, POS = asia1$POS, ALT = altCount1),
    file = paste("mixedFieldSamplePH0063-PH0193.", suffix, ".25v75.alt", sep = ""), row.names = F, quote=F, col.names=T, sep="\t")

WSAF2 = PH0064 * 0.75 + PH0193 *0.25
includeErrorWSAF2 = WSAF2*(1-err)+(1-WSAF2)*err
altCount2 = rbinom(n.loci, tmpTotalCoverage, includeErrorWSAF2)
refCount2 = tmpTotalCoverage - altCount2

write.table(data.frame(CHROM = asia1$CHROM, POS = asia1$POS, REF = refCount2),
    file = paste("mixedFieldSamplePH0063-PH0193.", suffix, ".75v25.ref", sep = ""), row.names = F, quote=F, col.names=T, sep="\t")

write.table(data.frame(CHROM = asia1$CHROM, POS = asia1$POS, ALT = altCount2),
    file = paste("mixedFieldSamplePH0063-PH0193.", suffix, ".75v25.alt", sep = ""), row.names = F, quote=F, col.names=T, sep="\t")

```

```
seed=1
suffix="noError"
dEploidCommon="-seed ${seed} -nSample 250 -rate 8 -burn 0.67 -k 2 -exportPostProb"
common1="-ref mixedFieldSamplePH0063-PH0193.${suffix}.25v75.ref -alt mixedFieldSamplePH0063-PH0193.${suffix}.25v75.alt -plaf asiaGroup1_PLAF.14.txt"
common2="-ref mixedFieldSamplePH0063-PH0193.${suffix}.75v25.ref -alt mixedFieldSamplePH0063-PH0193.${suffix}.75v25.alt -plaf asiaGroup1_PLAF.14.txt"
panel1="asia1.14.panel.txt"
panel2="asia1-2.14.panel.txt"
panel3="labStrains.14.panel.txt"

prefix1panel1="25v75panel1.${suffix}"
prefix1panel2="25v75panel2.${suffix}"
prefix1panel3="25v75panel3.${suffix}"

prefix2panel1="75v25panel1.${suffix}"
prefix2panel2="75v25panel2.${suffix}"
prefix2panel3="75v25panel3.${suffix}"

dEploid ${common1} ${dEploidCommon} -panel ${panel1} -o ${prefix1panel1}
R --slave "--args ${common1} -dEprefix ${prefix1panel1} -o ${prefix1panel1}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common1} ${dEploidCommon} -panel ${panel2} -o ${prefix1panel2}
R --slave "--args ${common1} -dEprefix ${prefix1panel2} -o ${prefix1panel2}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common1} ${dEploidCommon} -panel ${panel3} -o ${prefix1panel3}
R --slave "--args ${common1} -dEprefix ${prefix1panel3} -o ${prefix1panel3}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common2} ${dEploidCommon} -panel ${panel1} -o ${prefix2panel1}
R --slave "--args ${common2} -dEprefix ${prefix2panel1} -o ${prefix2panel1}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common2} ${dEploidCommon} -panel ${panel2} -o ${prefix2panel2}
R --slave "--args ${common2} -dEprefix ${prefix2panel2} -o ${prefix2panel2}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common2} ${dEploidCommon} -panel ${panel3} -o ${prefix2panel3}
R --slave "--args ${common2} -dEprefix ${prefix2panel3} -o ${prefix2panel3}" < ~/DEploid/utilities/interpretDEploid.r
```
