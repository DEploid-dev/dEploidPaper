# This experiment is pretty much the same as the DEploid paper figure 4 experiment.

Population groups (7), proportion pair (4), experiment (100)


Download panels

```bash

while read group; do
scp rescomp:/well/mcvean/joezhu/pf3k/pf3k_5_1_final/${group}.gz .
done < groups

while read group; do
zcat fullPanels/${group}.gz | head -1 > ${group}_Pf3D7_14_v3.csv
zcat fullPanels/${group}.gz | grep Pf3D7_14_v3 >> ${group}_Pf3D7_14_v3.csv
done < groups

while read group; do
cat plafs/${group}_PLAF.txt | head -1 > ${group}_Pf3D7_14_v3_PLAF.txt
cat plafs/${group}_PLAF.txt | grep Pf3D7_14_v3 >> ${group}_Pf3D7_14_v3_PLAF.txt
done < groups

```
