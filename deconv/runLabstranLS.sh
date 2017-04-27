#!/bin/bash
currentDir="~/dEploidPaper/deconv/"
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
#$ -N ${sample}Ex
#$ -t 1-30


#plaf=${root}/labStrains_PLAF.txt
#panel=/well/mcvean/joezhu/pf3k/pf3k_5_1_final/panels/labStrainsPanelFinal.csv
#excludeAt=/users/mcvean/joezhu/pf3k_mixed_infection/fieldSamples/clusters/labStrainsExclude.txt
#vcf=${root}/vcf/${sample}.wg.vcf.gz

plaf=${currentDir}labStrains.eg.PLAF.txt
panel=${currentDir}labStrains.eg.panel.txt
excludeAt=${currentDir}exclude.txt
vcf=${currentDir}${sample}.vcf.gz

prefix=${sample}_seed\${SGE_TASK_ID}k$@
common=\"-vcf \${vcf} -plaf \${plaf} -exclude \${excludeAt} -o \${prefix}\"
dEploidCommon=\"\${common} -seed \${SGE_TASK_ID} -nSample 500 -rate 8 -burn 0.67\"
rCommon=\"\${common} -dEprefix \${prefix}\"

(time dEploid \${dEploidCommon} -panel \${panel} -k $@) &> ${root}/dEploidOut/\${sample}/\${prefix}.time
#(time dEploid \${dEploidCommon} -noPanel -k $@) &> ${root}/dEploidOut/\${sample}/\${prefix}.time
#initialProp=\$( cat \${prefix}.prop | tail -1 | sed -e \"s/\t/ /g\" )
#(time dEploid \${dEploidCommon} -panel \${panel} -initialP \${initialProp} -k $@) &> ${root}/dEploidOut/\${sample}/\${prefix}.time

dEploid \${common} -panel \${panel} -painting \${prefix}.hap -o \${prefix} -initialP \${initialProp}

R --slave \"--args \${rCommon} \" < ~/DEploid/utilities/interpretDEploid.r

" > ${sample}k$@.sh
    qsub ${sample}k$@.sh
}

while read sample ; do
    run_dEploid 2
    run_dEploid 3
    run_dEploid 4
    run_dEploid 5
done < labSampleNames


