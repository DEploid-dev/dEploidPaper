#!/bin/bash
#$ -cwd
#$ -V
#$ -P mcvean.prjc -q short.qc
#$ -e ErrFiles
#$ -o OutFiles
#$ -N subSample
#$ -t 1-15
seed=${SGE_TASK_ID}
root="/well/mcvean/joezhu/pf3k/pf3k_5_1_final/subSampleCov/"
common="-panel labStrains.eg.panel.txt -plaf labStrains.eg.PLAF.txt -k 2 -seed ${seed}"

sample="PG0402-C"
#sample="PG0406-C"

dEploid ${common} -vcf ${sample}.subSample20.vcf.gz -o ${root}${sample}.seed${seed}.subSample20.lab.out
dEploid ${common} -vcf ${sample}.subSample50.vcf.gz -o ${root}${sample}.seed${seed}.subSample50.lab.out
dEploid ${common} -vcf ${sample}.subSample80.vcf.gz -o ${root}${sample}.seed${seed}.subSample80.lab.out
dEploid ${common} -vcf ${sample}.subSample100.vcf.gz -o ${root}${sample}.seed${seed}.subSample100.lab.out


common="-panel asiaAfirca.eg.panel.txt -plaf labStrains.eg.PLAF.txt -k 2 -seed ${seed}"

dEploid ${common} -vcf ${sample}.subSample20.vcf.gz -o ${root}${sample}.seed${seed}.subSample20.asiaAfirca.out
dEploid ${common} -vcf ${sample}.subSample50.vcf.gz -o ${root}${sample}.seed${seed}.subSample50.asiaAfirca.out
dEploid ${common} -vcf ${sample}.subSample80.vcf.gz -o ${root}${sample}.seed${seed}.subSample80.asiaAfirca.out
dEploid ${common} -vcf ${sample}.subSample100.vcf.gz -o ${root}${sample}.seed${seed}.subSample100.asiaAfirca.out
