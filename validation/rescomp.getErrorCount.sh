#!/bin/bash


while read sample ; do
for seed in $(seq 1 15);
do
prefix="/well/mcvean/joezhu/pf3k/pf3k_5_1_final/dEploidOut/"${sample}/${sample}_seed${seed}k2
R --slave "--args ${prefix}" < tmp2.r
done
done < labSampleNames2Strains

while read sample ; do
for seed in $(seq 1 15);
do
prefix="/well/mcvean/joezhu/pf3k/pf3k_5_1_final/dEploidOut/"${sample}/${sample}_seed${seed}k3
R --slave "--args ${prefix}" < tmp3.r
done
done < labSampleNames3Strains
