#!/bin/bash


#grab the header
head -n 10000 $1 | grep "^#" >header
#split into chunks by chromosome
grep -v "^#" $1 | cut -f 1 | sort | uniq | while read i;do
cat header > $i.vcf
tabix $1.gz $i > tmp1.vcf
done
rm -f header


