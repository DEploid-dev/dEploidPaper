#!/bin/bash
while read sample ; do
R --slave "--args ${sample}" < tmp2.r
done < labSampleNames2Strains
