#!/bin/bash
#$ -cwd
#$ -V
#$ -e ErrFiles
#$ -o OutFiles
#$ -N sample3s
#$ -t 1-3

sampleName=$( head -${SGE_TASK_ID} labSampleNames3Strains | tail -1  )

#for seed in $(seq 1 5);
#do
#R --slave "--args ${sampleName} 3 ${seed}" < run_pfmix.r
#done
R --slave "--args ${sampleName} 3 2" < run_pfmix.r
