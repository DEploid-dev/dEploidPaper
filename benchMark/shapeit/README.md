```
cp ../../validation/PG0*eg.vcf .

for file in ./*.eg.vcf ; do bgzip ${file};  done
for file in ./*.eg.vcf.gz ; do tabix -p vcf ${file}; done

vcf-merge PG0389-C.eg.vcf.gz PG0390-C.eg.vcf.gz PG0391-C.eg.vcf.gz PG0392-C.eg.vcf.gz PG0393-C.eg.vcf.gz PG0394-C.eg.vcf.gz PG0395-C.eg.vcf.gz PG0396-C.eg.vcf.gz PG0397-C.eg.vcf.gz PG0398-C.eg.vcf.gz PG0399-C.eg.vcf.gz PG0400-C.eg.vcf.gz PG0401-C.eg.vcf.gz PG0402-C.eg.vcf.gz PG0403-C.eg.vcf.gz PG0404-C.eg.vcf.gz PG0405-C.eg.vcf.gz PG0406-C.eg.vcf.gz PG0407-C.eg.vcf.gz PG0408-C.eg.vcf.gz PG0409-C.eg.vcf.gz PG0410-C.eg.vcf.gz PG0411-C.eg.vcf.gz PG0412-C.eg.vcf.gz PG0413-C.eg.vcf.gz PG0414-C.eg.vcf.gz PG0415-C.eg.vcf.gz | bgzip -c > lab.vcf.gz


gunzip -d lab.vcf.gz
head -n 10000 lab.vcf | grep "^#" >header
#split into chunks by chromosome
grep -v "^#" lab.vcf | cut -f 1 | sort | uniq | while read i;do
cat header > $i.vcf
tabix lab.vcf.gz $i >> $i.vcf
done
rm -f header


on cluster

while read chr; do echo ${chr}; /apps/well/shapeit/2.r790/shapeit -V tmp/${chr}.vcf -O tmp/${chr}; done < tmp/vcfFiles


while read chr; do echo ${chr}; cat ${chr}.haps >> shape2.haps; done < vcfFiles


cut -d ' ' -f 1,3,6,7 shape2.haps > PG0389-C.haps
cut -d ' ' -f 1,3,8,9 shape2.haps > PG0390-C.haps
cut -d ' ' -f 1,3,10,11 shape2.haps > PG0391-C.haps
cut -d ' ' -f 1,3,12,13 shape2.haps > PG0392-C.haps
cut -d ' ' -f 1,3,14,15 shape2.haps > PG0393-C.haps
cut -d ' ' -f 1,3,16,17 shape2.haps > PG0394-C.haps
cut -d ' ' -f 1,3,18,19 shape2.haps > PG0395-C.haps
cut -d ' ' -f 1,3,20,21 shape2.haps > PG0396-C.haps
cut -d ' ' -f 1,3,22,23 shape2.haps > PG0397-C.haps
cut -d ' ' -f 1,3,24,25 shape2.haps > PG0398-C.haps
cut -d ' ' -f 1,3,26,27 shape2.haps > PG0399-C.haps
cut -d ' ' -f 1,3,28,29 shape2.haps > PG0400-C.haps
cut -d ' ' -f 1,3,30,31 shape2.haps > PG0401-C.haps
cut -d ' ' -f 1,3,32,33 shape2.haps > PG0402-C.haps
cut -d ' ' -f 1,3,34,35 shape2.haps > PG0403-C.haps
cut -d ' ' -f 1,3,36,37 shape2.haps > PG0404-C.haps
cut -d ' ' -f 1,3,38,39 shape2.haps > PG0405-C.haps
cut -d ' ' -f 1,3,40,41 shape2.haps > PG0406-C.haps
cut -d ' ' -f 1,3,42,43 shape2.haps > PG0407-C.haps
cut -d ' ' -f 1,3,44,45 shape2.haps > PG0408-C.haps
cut -d ' ' -f 1,3,46,47 shape2.haps > PG0409-C.haps
cut -d ' ' -f 1,3,48,49 shape2.haps > PG0410-C.haps
cut -d ' ' -f 1,3,50,51 shape2.haps > PG0411-C.haps
cut -d ' ' -f 1,3,52,53 shape2.haps > PG0412-C.haps
cut -d ' ' -f 1,3,54,55 shape2.haps > PG0413-C.haps
cut -d ' ' -f 1,3,56,57 shape2.haps > PG0414-C.haps
cut -d ' ' -f 1,3,58,59 shape2.haps > PG0415-C.haps


```
