#!/usr/bin/env python
# USAGE: convert_pf3k_vcf_to_beagle_v3.py VCFFILE.vcf/VCFFILE.vcf.gz
# for chr in $(seq 1 14); do ./convert_pf3k_vcf_to_beagle_impute.py /well/mcvean/joezhu/pf3k/Ghana/ByChrom/Ghana_Chrom_${chr}.vcf.gz;  if [ $? -ne 0 ]; then echo " converting to beagle impute done" | mail -s "chr ${chr} failed" joezhu@well.ox.ac.uk ; else echo " converting to beagle impute done" | mail -s "chr ${chr} worked" joezhu@well.ox.ac.uk; fi done;
"""
============================================================
v3 ( 20150306)
line 144, assert ( category_count[3] > 0 ) # There is missing data at every sites
This was true for Ghana sample, but not necessary for NEC LAOS combined data

v3 ( 201502..)
skip homo sites, and singletons.
Count the frequence of genenotype 0, 0/1, 1, and '.'
skip SNP of f(0) + f(0/1) < 2 # This will exclude sites are like [0, 1, 472, 28], [1, 0, 472, 28], [0, 0, 472, 28]
and f(1) + f(0/1) < 2 # This will exclude sites are like [484, 0, 0, 17], [479, 1, 0, 21], [472, 0, 1, 28]

============================================================
v2 ( 20150129 -- 20150131 )
beagle still does not like missing sites
skip line with missing sites as the following
Pf3D7_05_v3 822827  .   G   .   .   PASS    .   GT:DP:AD:DF:DR  .:41:41:13:28

However, for the site
Pf3D7_01_v3	98871	.	C	A	.	PASS	.	GT:DP:AD:DF:DR	.:2:2,0:2,0:0,0
This is treated as missing data: 0/0, and the count for missing site is 0

Found missing entry at
Pf3D7_05_v3 822827  .   G   .   .   PASS    .   GT:DP:AD:DF:DR  .:41:41:13:28
 - modify the conversion program, if it is a missing variant. treat as "0/0", and count is the count, then 0.

============================================================
v1 (20150128)
 - IGNORE ALL VARIANT ENTRIES THAT HAVE MULTIPLE ALLELES
 - genotype of "0" is expressed as "0/0"
 - genotype of "0/1" is expressed as "0/1"
 - genotype of "1" is expressed as "1/1"
 - missing genotype "." is treated as "0/0"
"""


import sys
import gzip
from math import log10

def convert_info ( INFO, REF, ALT, check_line = "" ):
    splitted_INFO = INFO.split(":")

    category = -1
    NEW_INFO = "X"
    if ( splitted_INFO[0] == "0/0" ):
        NEW_INFO = REF
        category = 0    # REF site
        #assert ( ref_freq > 0 and alt_freq < 2 ) #assert ( ref_freq > 0 and alt_freq == 0 )
    elif ( splitted_INFO[0] == "0/1" ):
        NEW_INFO = "N"
        category = 1    # Het site
        #assert ( ref_freq > 0 and alt_freq > 0 )
    elif ( splitted_INFO[0] == "1/1" ):
        #GT = "1/1"
        NEW_INFO = ALT
        category = 2    # ALT site
        #assert ( ref_freq < 2 and alt_freq > 0 ) #
    elif ( splitted_INFO[0] == "./." ):
        #GT = "./."  # This is different from ../convert_pf3k_vcf_to_beagle.py
        NEW_INFO = "X"
        category = 3    # Missing site
        #assert ( ref_freq < 5 and alt_freq < 5 )      #assert ( ref_freq == 0 and alt_freq == 0 )
    else:
        print check_line
        raise Exception ( "Undefined GENOTYPE: " + splitted_INFO[0] + ", of " + INFO )

    assert ( category != -1 )  # category must have been assigned with a value

    return NEW_INFO

if __name__ == "__main__":
    #current_VERSION = "." + open("VERSION","r").readline().strip()
    input_file_name = sys.argv[1]
    extra_str = sys.argv[2] if len( sys.argv ) == 3 else ""

    print "Input file : ", input_file_name

    if ( "gz" in input_file_name ):
        input_file = gzip.open( input_file_name, "r" )
        output_file_name = input_file_name.strip(".vcf.gz") + extra_str  + ".txt"
        output_file = open( output_file_name, "w" )
    #else:
        #input_file = open( input_file_name, "r" )
        #filtered_file_name = input_file_name.strip(".vcf") + extra_str  + ".pre.gl.vcf"
        #filtered_output_file = open( filtered_file_name, "w" )
        #output_file_name = input_file_name.strip(".vcf") + extra_str  + ".gl.vcf"
        #output_file = open( output_file_name, "w" )

    #print "Filtered file: ", filtered_file_name
    print "Output file: ", output_file_name

    counter = 0
    nsam = 0
    sampleNames = []
    sequence = []
    for line in input_file:
        fields = line.split()
        # comment lines
        if ( line[0] == "#" ):
            # write the comment lines into the output file
            #output_file.write(line)
            #filtered_output_file.write ( line )
            if ( fields[0] == "#CHROM" ) :
                sampleNames = fields[9:len(fields)]
                nsam =  len(sampleNames)
                for i in range(nsam):
                    sequence.append("")
            continue

        # conversion are carried out in non-comment lines
        # sample starts from the 9th column

        if len(fields[4]) > 1:
            # skip the line with multiple ALT, e.g.:
            # Pf3D7_01_v3     111774  .       G       T,A     .       PASS    .       GT:DP:AD:DF:DR  .:75:56,19,0:30,0,0:26,0,0
            continue
        elif (fields[4] == "."):
            # skip missing data
            # Pf3D7_05_v3 822827  .   G   .   .   PASS    .   GT:DP:AD:DF:DR  .:41:41:13:28
            continue

        category_count = [ 0, 0, 0, 0 ]
        # Sample starts from the 9th field til the end
        for sample_i in range(9, len(fields)):
            sample_i_info = fields[sample_i]
            new_sample_i_info = convert_info (sample_i_info, fields[3], fields[4], line)
            sequence[sample_i-9] += new_sample_i_info

        nsam = len(fields) - 9
        counter += 1

    for sample_i in range(nsam):
        tmp = sampleNames[sample_i] + " " + sequence[sample_i] +"\n"
        output_file.write ( tmp )

    input_file.close()
    output_file.close()

