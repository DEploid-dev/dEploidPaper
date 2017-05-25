#!/usr/bin/env python


import sys
import gzip
from math import log10


if __name__ == "__main__":
    #current_VERSION = "." + open("VERSION","r").readline().strip()
    input_file_name = sys.argv[1]
    extra_str = sys.argv[2] if len( sys.argv ) == 3 else ""

    print "Input file : ", input_file_name

    if ( "gz" in input_file_name ):
        input_file = gzip.open( input_file_name, "r" )
        output_file_name = input_file_name.strip(".gl.out.vcf.gz") + extra_str  + ".gt.txt"
        output_file = open( output_file_name, "w" )
        #output_file_name = input_file_name.strip(".vcf.gz") + extra_str  + ".gl.vcf.gz"
        #output_file = gzip.open( output_file_name, "w" )

    print "output_file_name file: ", output_file_name
    #print "Output file: ", output_file_name
    output_file.write("CHROM\tPOS\th1\th2\n")
    for line in input_file:
        # comment lines
        if ( line[0] == "#" ):
            # write the comment lines into the output file
            #output_file.write(line)
            #filtered_output_file.write ( line )
            continue

        # conversion are carried out in non-comment lines
        # sample starts from the 9th column
        fields = line.split()
        newline = fields[0] + '\t' + fields[1] + '\t' + "\t".join(fields[9].split(":")[0].split("|")) + "\n"
        output_file.write ( newline )
        ##assert ( fields[8] == "GT:DP:AD:DF:DR" )

        ## replace Format by GT:GL
        #fields[8] = "GT:GL"

        #if len(fields[4]) > 1:
            ## skip the line with multiple ALT, e.g.:
            ## Pf3D7_01_v3     111774  .       G       T,A     .       PASS    .       GT:DP:AD:DF:DR  .:75:56,19,0:30,0,0:26,0,0
            #continue
        #elif (fields[4] == "."):
            ## skip missing data
            ## Pf3D7_05_v3 822827  .   G   .   .   PASS    .   GT:DP:AD:DF:DR  .:41:41:13:28
            #continue

        #category_count = [ 0, 0, 0, 0 ]
        ## Sample starts from the 9th field til the end
        #for sample_i in range(9, len(fields)):
            #sample_i_info = fields[sample_i]
            #( new_sample_i_info, current_category ) = convert_info (sample_i_info, line)
            #fields[sample_i] = new_sample_i_info
            #category_count[current_category] += 1

        #nsam = len(fields) - 9

        #assert ( sum(category_count) == nsam ) # Check if all the possiblilities were considered
        ##assert ( category_count[3] > 0 ) # There is missing data at every sites

        #if ( category_count[0] + category_count[1] < 2 ):
            ## This will exclude sites are like [0, 1, 472, 28], [1, 0, 472, 28], [0, 0, 472, 28]
            ##print "case 1", category_count
            #continue
        #elif ( category_count[2] + category_count[1] < 2 ):
            ## This will exclude sites are like [484, 0, 0, 17], [479, 1, 0, 21], [472, 0, 1, 28]
            ##print "case 2", category_count
            #continue
        #else:
            #newline = '\t'.join( fields )+"\n"
            #output_file.write ( newline )
            #filtered_output_file.write ( line )
            #continue

        #print line, newline
        #raise Exception ( "Should never reach here!!!" )

    input_file.close()
    #filtered_output_file.close()
    output_file.close()
