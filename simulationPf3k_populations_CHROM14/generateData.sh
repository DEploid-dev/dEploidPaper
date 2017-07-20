#!/bin/bash
#$ -cwd
#$ -V
#$ -P mcvean.prjc -q short.qc
#$ -e ErrFiles
#$ -o OutFiles
#$ -N generateData
#$ -t 1-100

./generateData.r ${SGE_TASK_ID}
