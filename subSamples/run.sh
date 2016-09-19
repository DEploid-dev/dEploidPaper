#!/bin/bash
#$ -cwd
#$ -V
#$ -P mcvean.prjc -q short.qc
#$ -e ErrFiles
#$ -o OutFiles
#$ -N subSample
#$ -t 1-15
root="/well/mcvean/joezhu/pf3k/pf3k_5_1_final/subSampleCov/"
common="-panel labStrains.14.panel.txt -plaf labStrains.14.PLAF.txt -k 2 -seed ${SGE_TASK_ID}"

#sample="PG0402-C"
sample="PG0406-C"

dEploid ${common} -vcf ${sample}.subSample20.vcf -o ${root}${sample}.subSample20.lab.out
dEploid ${common} -vcf ${sample}.subSample50.vcf -o ${root}${sample}.subSample50.lab.out
dEploid ${common} -vcf ${sample}.subSample80.vcf -o ${root}${sample}.subSample80.lab.out
dEploid ${common} -vcf ${sample}.subSample100.vcf -o ${root}${sample}.subSample100.lab.out

