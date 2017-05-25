#!/bin/bash
seed=1

#plaf=labStrains.14.PLAF.txt
#vcf=PG0412-C.14.vcf
#common="-vcf ${vcf} -plaf ${plaf} -seed ${seed} -k 2"

#(time -p dEploid ${common} -panel labStrains.14.panel.txt -o PG0412-C.14.labStrains) &> PG0412-C.14.labStrains.time
#(time -p dEploid ${common} -panel asiaAfrica2.14.panel.txt -o PG0412-C.14.asiaAfrica2) &> PG0412-C.14.asiaAfrica2.time
#(time -p dEploid ${common} -panel asiaAfrica.14.panel.txt -o PG0412-C.14.asiaAfrica) &> PG0412-C.14.asiaAfrica.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3.14.panel.txt -o PG0412-C.14.asiaAfrica_hb3) &> PG0412-C.14.asiaAfrica_hb3.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3_7g8.14.panel.txt -o PG0412-C.14.asiaAfrica_hb3_7g8) &> PG0412-C.14.asiaAfrica_hb3_7g8.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3_7g8_dd2.14.panel.txt -o PG0412-C.14.asiaAfrica_hb3_7g8_dd2) &> PG0412-C.14.asiaAfrica_hb3_7g8_dd2.time

plaf=labStrains.1314.PLAF.txt
vcf=PG0412-C.1314.vcf
common="-vcf ${vcf} -plaf ${plaf} -seed ${seed} -k 2"
(time -p dEploid ${common} -panel labStrains.1314.panel.txt -o PG0412-C.1314.labStrains) &> PG0412-C.1314.labStrains.time
echo "(time -p dEploid ${common} -panel labStrains.1314.panel.txt -o PG0412-C.1314.labStrains) &> PG0412-C.1314.labStrains.time"
#(time -p dEploid ${common} -panel asiaAfrica2.1314.panel.txt -o PG0412-C.1314.asiaAfrica2) &> PG0412-C.1314.asiaAfrica2.time
#(time -p dEploid ${common} -panel asiaAfrica.1314.panel.txt -o PG0412-C.1314.asiaAfrica) &> PG0412-C.1314.asiaAfrica.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3.1314.panel.txt -o PG0412-C.1314.asiaAfrica_hb3) &> PG0412-C.1314.asiaAfrica_hb3.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3_7g8.1314.panel.txt -o PG0412-C.1314.asiaAfrica_hb3_7g8) &> PG0412-C.1314.asiaAfrica_hb3_7g8.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3_7g8_dd2.1314.panel.txt -o PG0412-C.1314.asiaAfrica_hb3_7g8_dd2) &> PG0412-C.1314.asiaAfrica_hb3_7g8_dd2.time

#plaf=labStrains.121314.PLAF.txt
#vcf=PG0412-C.121314.vcf
#common="-vcf ${vcf} -plaf ${plaf} -seed ${seed} -k 2"
#(time -p dEploid ${common} -panel labStrains.121314.panel.txt -o PG0412-C.121314.labStrains) &> PG0412-C.121314.labStrains.time
#(time -p dEploid ${common} -panel asiaAfrica2.121314.panel.txt -o PG0412-C.121314.asiaAfrica2) &> PG0412-C.121314.asiaAfrica2.time
#(time -p dEploid ${common} -panel asiaAfrica.121314.panel.txt -o PG0412-C.121314.asiaAfrica) &> PG0412-C.121314.asiaAfrica.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3.121314.panel.txt -o PG0412-C.121314.asiaAfrica_hb3) &> PG0412-C.121314.asiaAfrica_hb3.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3_7g8.121314.panel.txt -o PG0412-C.121314.asiaAfrica_hb3_7g8) &> PG0412-C.121314.asiaAfrica_hb3_7g8.time
#(time -p dEploid ${common} -panel asiaAfrica_hb3_7g8_dd2.121314.panel.txt -o PG0412-C.121314.asiaAfrica_hb3_7g8_dd2) &> PG0412-C.121314.asiaAfrica_hb3_7g8_dd2.time

