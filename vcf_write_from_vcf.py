#!/usr/bin/env python

##############################################################
#                                                            #
#   Template script to write VCF files from another          #
#   VCF file with the vcfpy library.                         #
#   (https://vcfpy.readthedocs.io/en/master/index.html#)     #
#                                                            #          
#   Here the FORMAT header lines of the input file are       #
#   turned in INFO header lines of the output file           #
#                                                            #
##############################################################
from datetime import date
from collections import OrderedDict
import argparse
import itertools
import sys
import shlex
import vcfpy


def mapping_to_str(mapping):
    """Convert mapping to string"""
    result = ["<"]
    for i, (key, value) in enumerate(mapping.items()):
        if i > 0:
            result.append(",")
        result += [key, "=", serialize_for_header(key, value)]
    result += [">"]
    return "".join(result)

def str_to_mapping(s):
    # remove leading '<'
    s = s.lstrip('<')
    # remove trailing '>'
    s = s.rstrip('>')
    
    # split on commas, avoiding to split quoted strings
    lexer = shlex.shlex(s, posix=True)
    lexer.whitespace_split = True
    lexer.whitespace = ','

    # pair: key=value 
    return OrderedDict(pair.split('=', 1) for pair in lexer)
    

def main():

    parser = argparse.ArgumentParser(description="vcf writer")
    parser.add_argument("input", metavar='input.vcf', action='store',
                            help='vcf file.', type=str)
    parser.add_argument("output", metavar='output.vcf', action='store',
                            help='vcf file.', type=str)

    args = parser.parse_args()

    outvcf = args.output
    invcf = args.input
    
    
    #########################
    #                       #
    #  creating the header  #
    #                       #
    #########################

    # The header can contain some fixed type lines (INFO, FORMAT, FILTER, etc.) and some general ones
    # In this case, the header will contain a line storing the name of the program which generated 
    # the file. We also add the information about the name of the sample which have been analyzed

    header = vcfpy.Header(lines=[vcfpy.HeaderLine(key="source", value=sys.argv[0]), vcfpy.HeaderLine(key="fileformat", value="VCFv4.3"), vcfpy.HeaderLine(key="fileDate", value=date.today().strftime("%d/%m/%Y")) ], samples=vcfpy.SamplesInfos(["Sample1"]))

    
    # adding format lines 
    header.add_format_line(OrderedDict([("ID", "GT"),("Number", "1"), ("Type","String"), ("Description", "Genotype")]))
    header.add_format_line(OrderedDict([("ID", "DP"),("Number", "1"), ("Type","Integer"), ("Description", "Filtered read depth (MAPQ > 30)")]))

    # read the input vcf
    with vcfpy.Reader.from_path(invcf) as reader:

        # get the FORMAT header lines of the input file
        # and convert them in INFO header lines of the output file 
        format_ids = reader.header.format_ids()
        for format_id in format_ids:
            format_line = reader.header.get_format_field_info(format_id)
            '''
            output example:
        
            FormatHeaderLine('FORMAT', '<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">', {'ID': 'AD', 'Number': 'R', 'Type': 'Integer', 'Description': 'Allelic depths for the ref and alt alleles in the order listed'})

            key = 'FORMAT'
            value = '<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
            '''
            header.add_info_line(str_to_mapping(format_line.value))
            #print(header)
    
    # write the vcf
    with vcfpy.Writer.from_path(outvcf, header) as writer:
        
        # creating one record
        record = vcfpy.Record(
                CHROM="1", POS=1, ID=[], REF="C", ALT=[vcfpy.Substitution(type_="SNV", value="G")], QUAL=None, FILTER=[], INFO={}, FORMAT=["GT", "DP"], calls=[vcfpy.Call("Sample1", OrderedDict([("GT", "0/1"),("DP", "47")]))]
       )
        #print(record)
        writer.write_record(record)

if __name__ == "__main__":
    sys.exit(main())



