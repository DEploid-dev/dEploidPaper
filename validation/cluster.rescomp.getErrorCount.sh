#!/bin/bash
#$ -cwd
#$ -V
#$ -P mcvean.prjc -q short.qc
#$ -e ErrFiles
#$ -o OutFiles
#$ -N ${sample}Ex
#$ -t 1-30

while read sample ; do
do
prefix="/well/mcvean/joezhu/pf3k/pf3k_5_1_final/dEploidOut/"${sample}/${sample}_seed${SGE_TASK_ID}k2
R --slave "--args ${prefix}" < getErrorCount.r
done < labSampleNames2Strains

