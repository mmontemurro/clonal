#!/usr/bin/env python

##############################################################
#                                                            #
#   Template script to write VCF files from scratch          #
#   with the vcfpy library                                   #
#   (https://vcfpy.readthedocs.io/en/master/index.html#)     #
#                                                            #
##############################################################

from collections import OrderedDict
import argparse
import itertools
import sys
import time

import vcfpy


def main():
    
    parser = argparse.ArgumentParser(description="vcf writer")
    parser.add_argument("output", metavar='output.vcf', action='store',
                            help='vcf file (without header).', type=str)

    args = parser.parse_args()

    outvcf = args.output

    header = vcfpy.Header(samples=vcfpy.SamplesInfos(["Sample1"]))
    
    # adding format lines 
    header.add_format_line(OrderedDict([("ID", "DP"),("Number", "1"), ("Type","Integer"), ("Description", "Filtered read depth (MAPQ > 30)")]))
    
    # writing the vcf
    with vcfpy.Writer.from_path(outvcf, header) as writer:
        
        # creating one record
        record = vcfpy.Record(
            CHROM="1", POS=1, ID=[], REF="C", ALT=[], QUAL=None, FILTER=[], INFO={}, FORMAT=["DP"],
            calls=[vcfpy.Call("Sample1", OrderedDict([("DP", "47")]))]
       )
        #record.add_format(key="GT")
        #record.calls.append(vcfpy.Call("Sample1", OrderedDict([("GT", "0|1")])))
        writer.write_record(record)
        

if __name__ == "__main__":
    sys.exit(main())
