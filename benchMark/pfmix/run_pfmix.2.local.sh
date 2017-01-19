#!/bin/bash
#$ -cwd
#$ -V
#$ -e ErrFiles
#$ -o OutFiles
#$ -N sample2s
#$ -t 1-24

sampleName=$( head -${SGE_TASK_ID} labSampleNames2Strains | tail -1  )

#for seed in $(seq 1 5);
#do
#R --slave "--args ${sampleName} 2 ${seed}" < run_pfmix.r
#done
#R --slave "--args ${sampleName} 2 2" < run_pfmix.r
R --slave "--args ${sampleName} 2 ${SGE_TASK_ID}" < run_pfmix.r
