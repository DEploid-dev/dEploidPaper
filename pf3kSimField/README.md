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

WSAF3 = PH0064 * 0.55 + PH0193 *0.45
includeErrorWSAF3 = WSAF3*(1-err)+(1-WSAF3)*err
altCount3 = rbinom(n.loci, tmpTotalCoverage, includeErrorWSAF3)
refCount3 = tmpTotalCoverage - altCount3

write.table(data.frame(CHROM = asia1$CHROM, POS = asia1$POS, REF = refCount3),
    file = paste("mixedFieldSamplePH0063-PH0193.", suffix, ".55v45.ref", sep = ""), row.names = F, quote=F, col.names=T, sep="\t")

write.table(data.frame(CHROM = asia1$CHROM, POS = asia1$POS, ALT = altCount3),
    file = paste("mixedFieldSamplePH0063-PH0193.", suffix, ".55v45.alt", sep = ""), row.names = F, quote=F, col.names=T, sep="\t")

```

```bash
seed=2
suffix="withError"
dEploidCommon="-seed ${seed} -nSample 250 -rate 8 -burn 0.67 -k 2"
common1="-ref mixedFieldSamplePH0063-PH0193.${suffix}.25v75.ref -alt mixedFieldSamplePH0063-PH0193.${suffix}.25v75.alt -plaf asiaGroup1_PLAF.14.txt"
common2="-ref mixedFieldSamplePH0063-PH0193.${suffix}.75v25.ref -alt mixedFieldSamplePH0063-PH0193.${suffix}.75v25.alt -plaf asiaGroup1_PLAF.14.txt"
common3="-ref mixedFieldSamplePH0063-PH0193.${suffix}.55v45.ref -alt mixedFieldSamplePH0063-PH0193.${suffix}.55v45.alt -plaf asiaGroup1_PLAF.14.txt"
panel1="asia1.14.panel.txt"
panel2="asia1-2.14.panel.txt"
panel3="labStrains.14.panel.txt"


prefix2panel1="75v25panel1.${suffix}"
prefix2panel2="75v25panel2.${suffix}"
prefix2panel3="75v25panel3.${suffix}"
prefix2panel4="75v25noPanel.${suffix}"

dEploid ${common2} ${dEploidCommon} -noPanel -o ${prefix2panel4}
R --slave "--args ${common2} -dEprefix ${prefix2panel4} -o ${prefix2panel4}" < ~/DEploid/utilities/interpretDEploid.r


dEploid ${common2} ${dEploidCommon} -panel ${panel1} -o ${prefix2panel1} -exportPostProb
R --slave "--args ${common2} -dEprefix ${prefix2panel1} -o ${prefix2panel1}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common2} ${dEploidCommon} -panel ${panel2} -o ${prefix2panel2} -exportPostProb
R --slave "--args ${common2} -dEprefix ${prefix2panel2} -o ${prefix2panel2}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common2} ${dEploidCommon} -panel ${panel3} -o ${prefix2panel3} -exportPostProb
R --slave "--args ${common2} -dEprefix ${prefix2panel3} -o ${prefix2panel3}" < ~/DEploid/utilities/interpretDEploid.r


dEploidCommon="-seed ${seed} -nSample 250 -rate 8 -burn 0.67 -k 2 -miss 0.001"
prefix2panel1="75v25panel1.${suffix}miss0001"
prefix2panel2="75v25panel2.${suffix}miss0001"
prefix2panel3="75v25panel3.${suffix}miss0001"
prefix2panel4="75v25noPanel.${suffix}miss0001"

dEploid ${common2} ${dEploidCommon} -noPanel -o ${prefix2panel4}
R --slave "--args ${common2} -dEprefix ${prefix2panel4} -o ${prefix2panel4}" < ~/DEploid/utilities/interpretDEploid.r


dEploid ${common2} ${dEploidCommon} -panel ${panel1} -o ${prefix2panel1} -exportPostProb
R --slave "--args ${common2} -dEprefix ${prefix2panel1} -o ${prefix2panel1}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common2} ${dEploidCommon} -panel ${panel2} -o ${prefix2panel2} -exportPostProb
R --slave "--args ${common2} -dEprefix ${prefix2panel2} -o ${prefix2panel2}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common2} ${dEploidCommon} -panel ${panel3} -o ${prefix2panel3} -exportPostProb
R --slave "--args ${common2} -dEprefix ${prefix2panel3} -o ${prefix2panel3}" < ~/DEploid/utilities/interpretDEploid.r





prefix3panel1="55v45panel1.${suffix}"
prefix3panel2="55v45panel2.${suffix}"
prefix3panel3="55v45panel3.${suffix}"
prefix3panel4="55v45noPanel.${suffix}"

dEploid ${common3} ${dEploidCommon} -noPanel -o ${prefix3panel4}
R --slave "--args ${common3} -dEprefix ${prefix3panel4} -o ${prefix3panel4}" < ~/DEploid/utilities/interpretDEploid.r


dEploid ${common3} ${dEploidCommon} -panel ${panel1} -o ${prefix3panel1} -exportPostProb
R --slave "--args ${common3} -dEprefix ${prefix3panel1} -o ${prefix3panel1}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common3} ${dEploidCommon} -panel ${panel2} -o ${prefix3panel2} -exportPostProb
R --slave "--args ${common3} -dEprefix ${prefix3panel2} -o ${prefix3panel2}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common3} ${dEploidCommon} -panel ${panel3} -o ${prefix3panel3} -exportPostProb
R --slave "--args ${common3} -dEprefix ${prefix3panel3} -o ${prefix3panel3}" < ~/DEploid/utilities/interpretDEploid.r

prefix1panel1="25v75panel1.${suffix}"
prefix1panel2="25v75panel2.${suffix}"
prefix1panel3="25v75panel3.${suffix}"

dEploid ${common1} ${dEploidCommon} -panel ${panel1} -o ${prefix1panel1} -exportPostProb
R --slave "--args ${common1} -dEprefix ${prefix1panel1} -o ${prefix1panel1}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common1} ${dEploidCommon} -panel ${panel2} -o ${prefix1panel2} -exportPostProb
R --slave "--args ${common1} -dEprefix ${prefix1panel2} -o ${prefix1panel2}" < ~/DEploid/utilities/interpretDEploid.r

dEploid ${common1} ${dEploidCommon} -panel ${panel3} -o ${prefix1panel3} -exportPostProb
R --slave "--args ${common1} -dEprefix ${prefix1panel3} -o ${prefix1panel3}" < ~/DEploid/utilities/interpretDEploid.r
```


```R
rm(list=ls())

ref = read.table("mixedFieldSamplePH0063-PH0193.noError.75v25.ref", header=T)$REF
alt = read.table("mixedFieldSamplePH0063-PH0193.noError.75v25.alt", header=T)$ALT

panel1 = read.table("asia1.14.panel.txt", header=T, check.names=F)
PH0064 = panel1[["PH0064-C"]]
PH0193 = panel1[["PH0193-C"]]

panel1haps = read.table("75v25panel1.withError.hap", header=T)
PH0064.wrong.at = which(panel1haps[,3] != PH0064)
PH0193.wrong.at = which(panel1haps[,4] != PH0193)

plaf = read.table("asiaGroup1_PLAF.14.txt", header=T)$PLAF

ref[PH0064.wrong.at]
alt[PH0064.wrong.at]
plaf[PH0064.wrong.at]
PH0064[PH0064.wrong.at]
PH0193[PH0064.wrong.at]

panel1[PH0064.wrong.at,]

ref[PH0193.wrong.at]
alt[PH0193.wrong.at]
plaf[PH0193.wrong.at]
PH0064[PH0193.wrong.at]
PH0193[PH0193.wrong.at]
panel1[PH0193.wrong.at,]


> ref[PH0064.wrong.at]
[1] 2
> alt[PH0064.wrong.at]
[1] 0
> plaf[PH0064.wrong.at]
[1] 0.592652
> PH0064[PH0064.wrong.at]
[1] 0
>
> ref[PH0193.wrong.at]
[1]  4  0 11
> alt[PH0193.wrong.at]
[1] 0 2 0
> plaf[PH0193.wrong.at]
[1] 0.0000000 0.1091019 0.3166564
> PH0193[PH0193.wrong.at]
[1] 0 0 1

```
