#!/bin/bash

while read sample;
do
./rewriteGT.py ${sample}.gt.vcf.gz
done < labSampleNames

