#!/usr/bin/env python

##############################################################
#                                                            #
#   Template script to write VCF files from scratch          #
#   with the vcfpy library                                   #
#   (https://vcfpy.readthedocs.io/en/master/index.html#)     #
#                                                            #
##############################################################
from datetime import date
from collections import OrderedDict
import argparse
import itertools
import sys
import time
import vcfpy


def main():
    
    parser = argparse.ArgumentParser(description="vcf writer")
    parser.add_argument("output", metavar='output.vcf', action='store',
                            help='vcf file.', type=str)

    args = parser.parse_args()

    outvcf = args.output

    #########################
    #                       #
    #  creating the header  #
    #                       #
    #########################

    # The header can contain some fixed type lines (INFO, FORMAT, FILTER, etc.) and some general ones
    # In this case, the header will contain a line storing the name of the program which generated 
    # the file. We also add the information about the name of the sample which have been analyzed

    header = vcfpy.Header(lines=[vcfpy.HeaderLine(key="source", value=sys.argv[0]), vcfpy.HeaderLine(key="fileformat", value="VCFv4.3"), vcfpy.HeaderLine(key="fileDate", value=date.today().strftime("%d/%m/%Y")) ], samples=vcfpy.SamplesInfos(["Sample1", "Sample2"]))
    
    # Tuples of valid entries -----------------------------------------------------
    #
    #: valid INFO value types
    # INFO_TYPES = ("Integer", "Float", "Flag", "Character", "String")
    #: valid FORMAT value types
    # FORMAT_TYPES = ("Integer", "Float", "Character", "String")
    #: valid values for "Number" entries, except for integers
    # VALID_NUMBERS = ("A", "R", "G", ".")
    #: header lines that contain an "ID" entry
    # LINES_WITH_ID = ("ALT", "contig", "FILTER", "FORMAT", "INFO", "META", "PEDIGREE", "SAMPLE")
    # Constants for "Number" entries ----------------------------------------------
    #
    #: number of alleles excluding reference
    # HEADER_NUMBER_ALLELES = "A"
    #: number of alleles including reference
    # HEADER_NUMBER_REF = "R"
    #: number of genotypes
    # HEADER_NUMBER_GENOTYPES = "G"
    #: unbounded number of values
    # HEADER_NUMBER_UNBOUNDED = "."

    # adding filter lines
    header.add_filter_line(OrderedDict([("ID", "PASS"),("Description", "All filters passed")]))

    # adding info lines
    header.add_info_line(OrderedDict([("ID", "DP"), ("Number", "1"), ("Type","Integer"), ("Description", "Raw read depth (without mapping quality filters)")]))
    header.add_info_line(OrderedDict([("ID", "MUT"), ("Number", "1"), ("Type","Integer"), ("Description", "States if the record mutation is supported (1) or not (0).")]))

    # adding format lines 
    header.add_format_line(OrderedDict([("ID", "GT"),("Number", "1"), ("Type","String"), ("Description", "Genotype")]))
    header.add_format_line(OrderedDict([("ID", "DP"),("Number", "1"), ("Type","Integer"), ("Description", "Filtered read depth (MAPQ > 30)")]))
    #header.add_format_line(OrderedDict([vcfpy.header.RESERVED_FORMAT["GT"]]))


    # adding contig lines
    header.add_contig_line(OrderedDict([("ID", "chr1"), ("length","248956422")]))
    
    # adding sample lines
    header.add_line(vcfpy.SampleHeaderLine.from_mapping(OrderedDict([("ID", "Sample1"),("Description", "Tumor")])))
    
    # writing the vcf
    with vcfpy.Writer.from_path(outvcf, header) as writer:
        
        # creating one record
        calls = []
        calls.append(vcfpy.Call("Sample1", OrderedDict([("GT", "0/1"),("DP", "47")])))
        calls.append(vcfpy.Call("Sample2", OrderedDict([("GT", "0/1"),("DP", "31")])))

        record = vcfpy.Record(
                CHROM="1", POS=1, ID=[], REF="C", ALT=[vcfpy.Substitution(type_="SNV", value="G")], QUAL=None, FILTER=["PASS"], INFO={"DP":"50", "MUT":0}, FORMAT=["GT","DP"],
            calls=calls
       )
        #record.add_format(key="GT")
        #record.calls.append(vcfpy.Call("Sample1", OrderedDict([("GT", "0|1")])))
        writer.write_record(record)
        

if __name__ == "__main__":
    sys.exit(main())
