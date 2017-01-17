#!/usr/bin/env python
# USAGE: filteringSNP.py VCFFILE.vcf/VCFFILE.vcf.gz
# eg: ./filtering.py SNP_INDEL_Pf3D7_02_v3.combined.filtered.vcf.gz
# INPUT: release 4 vcf
# OUTPUT:
#
"""
The object of this script is to:
1. Filter out the bad sites, and put into a seperate file.
    1.1 multi alt sites
        # Pf3D7_01_v3     111774  .       G       T,A     .       PASS    .       GT:DP:AD:DF:DR  .:75:56,19,0:30,0,0:26,0,0

    1.2 bad missing calls, such as
        # Pf3D7_05_v3 822827  .   G   .   .   PASS    .   GT:DP:AD:DF:DR  .:41:41:13:28
        However, for the site
        Pf3D7_01_v3	98871	.	C	A	.	PASS	.	GT:DP:AD:DF:DR	.:2:2,0:2,0:0,0
        This is treated as missing data: ., and the count for missing site is 0

2. filtering out homo sites, and singletons.
    2.1 Count the frequence of genenotype 0, 0/1, 1, and '.'
    2.2 skip SNP, such that of f(0) + f(0/1) < 2 # This will exclude sites are like [0, 1, 472, 28], [1, 0, 472, 28], [0, 0, 472, 28]
                  and f(1) + f(0/1) < 2 # This will exclude sites are like [484, 0, 0, 17], [479, 1, 0, 21], [472, 0, 1, 28]
"""

import sys
import gzip
from math import log10

def checkFrequency ( INFO, check_line = "" ):
    splitted_INFO = INFO.split(":")
    # The following line does not work is because that there might be multiple ALTs,
    # splitted_INFO[2] hence has multiple reads, in the case of missing data, this splitted_INFO[2] has length 1,
    # which is the read depth of the REF
    #[ ref_freq, alt_freq ] = [ int(x) for x in splitted_INFO[2].split(",") ]

    # version 2_1 "GT:DP:AD:DF:DR"
    #allele_freq = [ int(x) for x in splitted_INFO[2].split(",") ]
    # version 3_1 "GT:AD:DF:DP:DR"
    #allele_freq = [ int(x) for x in splitted_INFO[1].split(",") ]
    # version 4_1, 5_1 "GT:AD:DP:GQ:PGT:PID:PL"
    allele_freq = [ int(x) for x in splitted_INFO[1].split(",") ]
    ref_freq = allele_freq[0]
    alt_freq = allele_freq[1]
    #alt_freq = 0 if missing else allele_freq[1]

    # Assertion check the sum of the allele frequency is equal to the read depth
    #assert ( ( ref_freq + alt_freq ) == int(splitted_INFO[1]) ) # version 2_1
    #assert ( ( ref_freq + alt_freq ) == int(splitted_INFO[3]) )  # version 3_1

#    assert ( ( ref_freq + alt_freq ) == int(splitted_INFO[2]) )  # version 4_1

    category = -1
    if ( splitted_INFO[0] == "0/0" ): # use "0/0" in 4_1, it was "0" in the older version
        category = 0    # REF site
        #print (splitted_INFO[0])
        #print (ref_freq)
        #print(alt_freq)
        #assert ( ref_freq > 0 and alt_freq < 7 ) #assert ( ref_freq > 0 and alt_freq == 0 )
    elif ( splitted_INFO[0] == "0/1" ):
        category = 1    # Het site
        #assert ( ref_freq > 0 and alt_freq > 0 )
    elif ( splitted_INFO[0] == "1/1" ): # use "1/1" in 4_1, it was "1" in the older version
        category = 2
        #print (splitted_INFO[0])
        #print (ref_freq)
        #print(alt_freq)
        #assert ( ref_freq < 3 and alt_freq > 0 ) #assert ( ref_freq == 0 and alt_freq > 0 )
    elif ( splitted_INFO[0] == "./." ): # use "./." in 4_1, it was "." in the older version
        category = 3    # Missing site
        #print (splitted_INFO[0])
        #print (ref_freq)
        #print(alt_freq)
        # before 4_1 release, the following assertion is valid.
        #assert ( ref_freq < 5 and alt_freq < 5 )      #assert ( ref_freq == 0 and alt_freq == 0 )
        # Now use
        #print (splitted_INFO[0])
        #print (ref_freq)
        #print(alt_freq)
        #print (splitted_INFO)
        #assert ( ref_freq < 23 and alt_freq < 23 )
        #assert ( ref_freq < 24 and alt_freq < 24 )
    else:
        print check_line
        raise Exception ( "Undefined GENOTYPE: " + splitted_INFO[0] + ", of " + INFO )
    if (category == -1):
        print("current ", INFO, splitted_INFO[0])
    assert ( category != -1 )  # category must have been assigned with a value

    return [category, ref_freq, alt_freq ]

if __name__ == "__main__":
    inputFileName = sys.argv[1]
    print (inputFileName)
    print "Input file : ", inputFileName
    #assert ( inputFileName.find("3_1") > 0 ) # Making sure the input file name contains string 3_1

    extra_str = sys.argv[2] if len( sys.argv ) == 3 else ""

    if ( "gz" in inputFileName ):
        inputFile = gzip.open( inputFileName, "r" )
        #filterKeptFileName = inputFileName.strip(".vcf.gz") + extra_str + ".filterKept.vcf.gz"
        #filterKept = gzip.open( filterKeptFileName, "w" )
        refFileName = inputFileName.strip(".vcf.gz") + extra_str + ".ref"
        altFileName = inputFileName.strip(".vcf.gz") + extra_str + ".alt"
        refFile = open(refFileName, "w")
        altFile = open(altFileName, "w")
        gtFileName = inputFileName.strip(".vcf.gz") + extra_str + ".gt"
        gtFile = open(gtFileName, "w")

        #filterOutFileName = inputFileName.strip(".vcf.gz") + extra_str + ".filterOut.vcf.gz"
        #filterOut  = gzip.open( filterOutFileName, "w" )
    else:
        inputFile = open( inputFileName, "r" )
        refFileName = inputFileName.strip(".vcf") + extra_str + ".ref"
        altFileName = inputFileName.strip(".vcf") + extra_str + ".alt"
        refFile = open(refFileName, "w")
        altFile = open(altFileName, "w")
        gtFileName = inputFileName.strip(".vcf") + extra_str + ".gt"
        gtFile = open(gtFileName, "w")
        #filterKeptFileName = inputFileName.strip(".vcf") + extra_str + ".filterKept.vcf"
        #filterKept = open( filterKeptFileName, "w" )
        #filterOutFileName = inputFileName.strip(".vcf") + extra_str + ".filterOut.vcf"
        #filterOut = open( filterOutFileName, "w" )

    #print "Filter kept file: ", filterKeptFileName
    #print "Filtered out file: ", filterOutFileName

    numberCommentLine = 0
    for line in inputFile:
        # Metedata
        if ( line[0] == "#" ):
            numberCommentLine += 1
            # write the comment lines into the output file
            #filterKept.write(line)
            #filterOut.write(line)
            if (line[1] != "#" ):
                fields = line.split()
                fields = fields[:2] + fields[9:]
                #fields[0] = "CHROM"
                refFile.write('\t'.join(fields)+"\n")
                altFile.write('\t'.join(fields)+"\n")
                gtFile.write('\t'.join(fields)+"\n")
            continue

        # conversion are carried out in non-comment lines
        # sample starts from the 9th column
        fields = line.split()

        #assert ( fields[8] == "GT:DP:AD:DF:DR" ) # version 2_1
        #assert ( fields[8] == "GT:AD:DF:DP:DR" ) # version 3_1

        #if ( fields[6].find("VQSLOD") != -1 & fields[6].find("Low_VQSLOD") == -1 ): # if PASS was not found, continue
            #print(fields[0]+" "+fields[1]+" "+fields[6])

        # 20151209 only using passing sites
        assert ( fields[6] == "PASS" )

        # at line 1174, it was #GT:AD:DP:GQ:PL
        field8Split = fields[8].split(":")
        assert ( field8Split[1] == "AD" )
        #if (fields[8] != "GT:AD:DP:GQ:PGT:PID:PL"):
            #filterOut.write ( line )
            #continue

            #print (line)
        #assert ( fields[8] == "GT:AD:DP:GQ:PGT:PID:PL" ) # version 4_1
        if len(fields[4]) > 1:
            filterOut.write ( line )
            # skip the line with multiple ALT, e.g.:
            # Pf3D7_01_v3     111774  .       G       T,A     .       PASS    .       GT:DP:AD:DF:DR  .:75:56,19,0:30,0,0:26,0,0
            continue
        elif (fields[4] == "."):
            filterOut.write ( line )
            # skip missing data
            # Pf3D7_05_v3 822827  .   G   .   .   PASS    .   GT:DP:AD:DF:DR  .:41:41:13:28
            continue
        elif len(fields[3]) > 1:
            filterOut.write ( line )
            # skip indels
            # Pf3D7_01_v3     37      .       GAACCCTAAACCCTGAACCCTA  G
            continue

        category_count = [ 0, 0, 0, 0 ]
        # Sample starts from the 9th field til the end
        refLine = [fields[0], fields[1]]
        altLine = [fields[0], fields[1]]
        gtLine = [fields[0], fields[1]]
        for sample_i in range(9, len(fields)):
            sample_i_info = fields[sample_i]
            tmpCategory = checkFrequency (sample_i_info, line)
            current_category = tmpCategory[0]
            category_count[current_category] += 1
            refLine.append(`tmpCategory[1]`)
            altLine.append(`tmpCategory[2]`)
            #tmpGt = "NA"
            if current_category == 0:
                gtLine.append ("0")
            elif current_category == 2:
                gtLine.append ("1")
            elif current_category == 1:
                gtLine.append ("0/1")
            else:
                gtLine.append (".")

        nsam = len(fields) - 9

        assert ( sum(category_count) == nsam ) # Check if all the possiblilities were considered
        #assert ( category_count[3] > 0 ) # There is missing data at every sites

        #filterKept.write ( line )
        refFile.write( '\t'.join(refLine) + "\n")
        altFile.write( '\t'.join(altLine) + "\n")
        gtFile.write( '\t'.join(gtLine) + "\n")

    #print ( "numberOfSingleton = ", numberOfSingleton )

    inputFile.close()
    #filterKept.close()
    #filterOut.close()
    gtFile.close()
    refFile.close()
    altFile.close()
