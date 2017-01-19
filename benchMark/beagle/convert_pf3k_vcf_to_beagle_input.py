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

def convert_info ( INFO, check_line = "" ):

    splitted_INFO = INFO.split(":")
    # The following line does not work is because that there might be multiple ALTs,
    # splitted_INFO[2] hence has multiple reads, in the case of missing data, this splitted_INFO[2] has length 1,
    # which is the read depth of the REF
    #[ ref_freq, alt_freq ] = [ int(x) for x in splitted_INFO[2].split(",") ]
    allele_freq = [ int(x) for x in splitted_INFO[1].split(",") ]
    ref_freq = allele_freq[0]
    alt_freq = allele_freq[1]
    #alt_freq = 0 if missing else allele_freq[1]
    # Assertion check the sum of the allele frequency is equal to the read depth
    #print(ref_freq)
    #print(alt_freq)
    #print(int(splitted_INFO[2]))
    #print(INFO)
    #assert ( ( ref_freq + alt_freq ) == int(splitted_INFO[2]) )
    error_rate = [ .01, .5, .99 ]
    gl_list = [ ( ref_freq * log10(1-x) + alt_freq * log10(x) ) for x in error_rate ]
    max_gl = max (gl_list)
    gl_list = [ round( x, 2 ) for x in gl_list ]
    #0/0:-0.03,-1.21,-5.00

    category = -1
    if ( splitted_INFO[0] == "0/0" ):
        GT = "0/0"
        category = 0    # REF site
        #assert ( ref_freq > 0 and alt_freq < 2 ) #assert ( ref_freq > 0 and alt_freq == 0 )
    elif ( splitted_INFO[0] == "0/1" ):
        GT = "0/1"
        category = 1    # Het site
        #assert ( ref_freq > 0 and alt_freq > 0 )
    elif ( splitted_INFO[0] == "1/1" ):
        GT = "1/1"
        category = 2    # ALT site
        #assert ( ref_freq < 2 and alt_freq > 0 ) #
        #print(ref_freq)
        #print(alt_freq)
        #assert ( ref_freq == 0 and alt_freq > 0 )
    elif ( splitted_INFO[0] == "./." ):
        GT = "./."  # This is different from ../convert_pf3k_vcf_to_beagle.py
        category = 3    # Missing site
        #assert ( ref_freq < 5 and alt_freq < 5 )      #assert ( ref_freq == 0 and alt_freq == 0 )
    else:
        print check_line
        raise Exception ( "Undefined GENOTYPE: " + splitted_INFO[0] + ", of " + INFO )

    assert ( category != -1 )  # category must have been assigned with a value

    NEW_INFO = GT + ":" + `gl_list[0]` + "," + `gl_list[1]` + "," + `gl_list[2]`
    return NEW_INFO, category

if __name__ == "__main__":
    #current_VERSION = "." + open("VERSION","r").readline().strip()
    input_file_name = sys.argv[1]
    extra_str = sys.argv[2] if len( sys.argv ) == 3 else ""

    print "Input file : ", input_file_name

    if ( "gz" in input_file_name ):
        input_file = gzip.open( input_file_name, "r" )
        filtered_file_name = input_file_name.strip(".vcf.gz") + extra_str  + ".pre.gl.vcf.gz"
        filtered_output_file = gzip.open( filtered_file_name, "w" )
        output_file_name = input_file_name.strip(".vcf.gz") + extra_str  + ".gl.vcf.gz"
        output_file = gzip.open( output_file_name, "w" )
    else:
        input_file = open( input_file_name, "r" )
        filtered_file_name = input_file_name.strip(".vcf") + extra_str  + ".pre.gl.vcf"
        filtered_output_file = open( filtered_file_name, "w" )
        output_file_name = input_file_name.strip(".vcf") + extra_str  + ".gl.vcf"
        output_file = open( output_file_name, "w" )

    print "Filtered file: ", filtered_file_name
    print "Output file: ", output_file_name

    for line in input_file:
        # comment lines
        if ( line[0] == "#" ):
            # write the comment lines into the output file
            output_file.write(line)
            filtered_output_file.write ( line )
            continue

        # conversion are carried out in non-comment lines
        # sample starts from the 9th column
        fields = line.split()

        #assert ( fields[8] == "GT:DP:AD:DF:DR" )

        # replace Format by GT:GL
        fields[8] = "GT:GL"

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
            ( new_sample_i_info, current_category ) = convert_info (sample_i_info, line)
            fields[sample_i] = new_sample_i_info
            category_count[current_category] += 1

        nsam = len(fields) - 9

        assert ( sum(category_count) == nsam ) # Check if all the possiblilities were considered
        #assert ( category_count[3] > 0 ) # There is missing data at every sites

        #if ( category_count[0] + category_count[1] < 2 ):
            ## This will exclude sites are like [0, 1, 472, 28], [1, 0, 472, 28], [0, 0, 472, 28]
            ##print "case 1", category_count
            #continue
        #elif ( category_count[2] + category_count[1] < 2 ):
            ## This will exclude sites are like [484, 0, 0, 17], [479, 1, 0, 21], [472, 0, 1, 28]
            ##print "case 2", category_count
            #continue
        #else:
        newline = '\t'.join( fields )+"\n"
        output_file.write ( newline )
        filtered_output_file.write ( line )
            #continue

        #print line, newline
        #raise Exception ( "Should never reach here!!!" )

    input_file.close()
    filtered_output_file.close()
    output_file.close()

