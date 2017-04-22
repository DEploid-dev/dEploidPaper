#!/bin/bash
currentDir="~/dEploidPaper/deconv"
root="/well/mcvean/joezhu/pf3k/pf3k_5_1_final/"
function run_dEploid {
    mkdir ${root}"dEploidOut/"${sample}
    cd ${root}"dEploidOut/"${sample}
    #rm -r ErrFiles OutFiles
    mkdir ErrFiles OutFiles

echo "
#!/bin/bash
#$ -cwd
#$ -V
#$ -P mcvean.prjc -q short.qc
#$ -e ErrFiles
#$ -o OutFiles
#$ -N ${sample}
#$ -t 1-15


plaf=${root}/labStrains_PLAF.txt
panel=/well/mcvean/joezhu/pf3k/pf3k_5_1_final/panels/labStrainsPanelFinal.csv
exludeAt=/users/mcvean/joezhu/pf3k_mixed_infection/fieldSamples/clusters/labStrainsExclude.txt
vcf=${root}/vcf/${sample}.wg.vcf.gz

prefix=${sample}_seed\${SGE_TASK_ID}k$@
common=\"-vcf \${vcf} -plaf \${plaf} -exclude \${exludeAt} -o \${prefix}\"
dEploidCommon=\"\${common} -panel \${panel} -seed \${SGE_TASK_ID} -nSample 250 -rate 8 -burn 0.67\"
rCommon=\"\${common} -dEprefix \${prefix}\"

(time dEploid \${dEploidCommon} -k $@) &> ${root}/dEploidOut/\${sample}/\${prefix}.time

initialProp=\$( cat \${prefix}.prop | tail -1 | sed -e \"s/\t/ /g\" )
dEploid \${common} -panel \${panel} -painting \${prefix}.hap -o \${prefix} -initialP \${initialProp}

R --slave \"--args \${rCommon} \" < ~/DEploid/utilities/interpretDEploid.r

" > ${sample}k$@.sh
    qsub ${sample}k$@.sh
}

while read sample ; do
    run_dEploid 3
    run_dEploid 5
done < labSampleNames


