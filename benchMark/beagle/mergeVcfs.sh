#!/bin/bash

for file in ./*.vcf ; do bgzip ${file}; done
sleep 60

for file in ./*.vcf.gz ; do bcftools index ${file} -f ; done
sleep 60

# split the process, then merge together
bcftools merge $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g') -Oz -o labMixedVcf.vcf.gz
bcftools index labMixedVcf.vcf.gz
