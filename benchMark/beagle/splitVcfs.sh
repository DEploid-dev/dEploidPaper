#!/bin/bash

while read sample;
do
bcftools view -s ${sample} labMixedV.gl.out.vcf.gz | bgzip > ${sample}.gt.vcf.gz
done < labSampleNames
