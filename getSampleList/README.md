```bash
zcat SNP_INDEL_Pf3D7_01_v3.high_quality_biallelic_snps.union.clonal.vcf.gz | head -400 | grep "#CHROM" > clonals
sed -e "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t//g" -e "s/.1\t/\n/g" -e "s/SenV092.05.1/SenV092.05/g" clonals > names

zcat SNP_INDEL_Pf3D7_01_v3.high_quality_biallelic_snps.union.mixed.vcf.gz | head -400 | grep "#CHROM" > mixeds
sed -e "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t//g" -e "s/.1\t/\n/g" -e "s/.2\t/\n/g" -e "s/.3\t/\n/g" -e "s/.4\t/\n/g" -e "s/.5\t/\n/g"  -e "s/SenT185.10.2/SenT185.10/g" mixeds >> names
```

```R
rm(list=ls())
a = read.table("names", header=F)
write.table(unique(a), file = "uniqueNames", row.names = F, col.names = F, quote = F)
```

```bash
cp ~/pf3k_mixed_infection/fieldSamples/clusters/a*_samples .
cat *_samples | sort > samples.sorted
sort uniqueNames > uniqueNames.sorted
comm -3 samples.sorted uniqueNames.sorted > samples_didnot_DEploid
```
