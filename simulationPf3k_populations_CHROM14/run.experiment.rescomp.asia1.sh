#!/bin/bash
#$ -cwd
#$ -V
#$ -P mcvean.prjc -q short.qc
#$ -e ErrFiles
#$ -o OutFiles
#$ -N asia1
#$ -t 1-100

panelSize=20
group="asia1"
common="-plaf plaf14/${group}_Pf3D7_14_v3_PLAF.txt -exclude ${group}excludeAt.txt"
prefix=${group}experiment${SGE_TASK_ID}
dEploidCommon="${common} -panel panels${panelSize}/${group}experiment${SGE_TASK_ID}.panel -seed ${SGE_TASK_ID}"

#################################################################################
suffix="15v85"
outPrefix="experimentOut/${prefix}.${suffix}panel${panelSize}"
speificCommon="-ref experimentCoverage20/${prefix}.${suffix}.ref -alt experimentCoverage20/${prefix}.${suffix}.alt -o ${outPrefix}"

dEploid ${speificCommon} ${dEploidCommon} -initialP 0.15 0.85
initialProp=$( cat ${outPrefix}.prop | tail -1 | sed -e "s/\t/ /g" )
dEploid ${speificCommon} ${common} -painting ${outPrefix}.hap -o ${outPrefix} -initialP ${initialProp} -panel truth20/${group}experiment${SGE_TASK_ID}.true.hap
interpretDEploid.r ${speificCommon} ${common} -dEprefix ${outPrefix}

#################################################################################
suffix="25v75"
outPrefix="experimentOut/${prefix}.${suffix}panel${panelSize}"
speificCommon="-ref experimentCoverage20/${prefix}.${suffix}.ref -alt experimentCoverage20/${prefix}.${suffix}.alt -o ${outPrefix}"

dEploid ${speificCommon} ${dEploidCommon} -initialP 0.25 0.75
initialProp=$( cat ${outPrefix}.prop | tail -1 | sed -e "s/\t/ /g" )
dEploid ${speificCommon} ${common} -painting ${outPrefix}.hap -o ${outPrefix} -initialP ${initialProp} -panel truth20/${group}experiment${SGE_TASK_ID}.true.hap
interpretDEploid.r ${speificCommon} ${common} -dEprefix ${outPrefix}

#################################################################################
suffix="35v65"
outPrefix="experimentOut/${prefix}.${suffix}panel${panelSize}"
speificCommon="-ref experimentCoverage20/${prefix}.${suffix}.ref -alt experimentCoverage20/${prefix}.${suffix}.alt -o ${outPrefix}"

dEploid ${speificCommon} ${dEploidCommon} -initialP 0.35 0.65
initialProp=$( cat ${outPrefix}.prop | tail -1 | sed -e "s/\t/ /g" )
dEploid ${speificCommon} ${common} -painting ${outPrefix}.hap -o ${outPrefix} -initialP ${initialProp} -panel truth20/${group}experiment${SGE_TASK_ID}.true.hap
interpretDEploid.r ${speificCommon} ${common} -dEprefix ${outPrefix}

#################################################################################
suffix="45v55"
outPrefix="experimentOut/${prefix}.${suffix}panel${panelSize}"
speificCommon="-ref experimentCoverage20/${prefix}.${suffix}.ref -alt experimentCoverage20/${prefix}.${suffix}.alt -o ${outPrefix}"

dEploid ${speificCommon} ${dEploidCommon} -initialP 0.45 0.55
initialProp=$( cat ${outPrefix}.prop | tail -1 | sed -e "s/\t/ /g" )
dEploid ${speificCommon} ${common} -painting ${outPrefix}.hap -o ${outPrefix} -initialP ${initialProp} -panel truth20/${group}experiment${SGE_TASK_ID}.true.hap
interpretDEploid.r ${speificCommon} ${common} -dEprefix ${outPrefix}

