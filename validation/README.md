This will contain comparison of four panels,
1. panel free
2. panel of 5 asiaGroup1 + 5 africaGroup1
3. panel of panel 2 + hb3
4. panel of panel 3 + 7g8
5. panel of panel 4 + dd2
On chrom 14 only, to show differences.

Use lines in `untitled.r`, to preprocess data.


Use `compareToPanel.r` to compare the results with the panel
```
dEploid -vcf PG0407-C.14.vcf -plaf labStrains.14.PLAF.txt -panel asiaPlus2.14.panel.txt -o PG0407-C.14.asiaPlus2 -exportPostProb

R --slave "--args -vcf PG0407-C.14.vcf -plaf labStrains.14.PLAF.txt -dEprefix PG0407-C.14.asiaPlus2 -o PG0407-C.14.asiaPlus2 " < ~/DEploid/utilities/interpretDEploid.r
```

