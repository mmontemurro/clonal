#!/usr/bin/env python
import os
import sys
import shlex
import vcfpy
import pysam
import argparse
from datetime import date
from collections import OrderedDict


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
    parser = argparse.ArgumentParser(description="Looks for a given set of SNPs whithin a bam file.")


    parser.add_argument("bam", metavar='sample.bam', action='store',
        help='BAM file.', type=str)

    parser.add_argument("vcf", metavar='file.vcf', action='store',
        help="VCF file storing SNPs.", type=str)

    parser.add_argument("sample_name", metavar='sample1', action='store',
                help="Sample identifier.", type=str)



    parser.add_argument("out_prefix", metavar="outdir/sample", action="store",
        help="Output VCF file prefix.", type=str)

    #parser.add_argument("--sample_name2", metavar='sample2', action='store',
    #                            help="Another sample name", type=str)

    args = parser.parse_args()
    bam= args.bam
    invcf = args.vcf
    sample = args.sample_name
    outvcf = args.out_prefix + ".snpseeker.vcf"

    '''
    if args.sample_name2:
    sample_name2 = args.sample_name2
    else:
    sample_name2 = null
    '''

    #read bam file
    samfile = pysam.AlignmentFile(bam, "rb")

    #build the header of the output vcf
    header_out = vcfpy.Header(lines=[vcfpy.HeaderLine(key="fileformat", value="VCFv4.3"), vcfpy.HeaderLine(key="source", value=sys.argv[0]), vcfpy.HeaderLine(key="fileDate", value=date.today().strftime("%d/%m/%Y")) ], samples=vcfpy.SamplesInfos([sample]))

    # sample header lines
    header_out.add_line(vcfpy.HeaderLine(key="SampleName", value=sample))
    '''
    if sample_name2 is not null:
    header_out.add_line(vcfpy.SampleHeaderLine.from_mapping(OrderedDict([("ID", sample_name2),("Description", "Second sample name")])))
    '''
    # info header lines
    header_out.add_info_line(OrderedDict([("ID", "SUPP"), ("Number", "1"), ("Type","Integer"), ("Description", "States if the mutation is supported (1) or not (0).")]))

    # adding format lines 
    header_out.add_format_line(OrderedDict([("ID", "GT"),("Number", "1"), ("Type","String"), ("Description", "Genotype (0/1, 0/0)")]))
    header_out.add_format_line(OrderedDict([("ID", "SDP"),("Number", "1"), ("Type","Integer"), ("Description", "Samtools read depth (secondary alignments, PCR duplicates, unppammed reads and reads not passing vendor QC are filtered)")]))
    header_out.add_format_line(OrderedDict([("ID", "DP"),("Number", "1"), ("Type","Integer"), ("Description", "Filtered read depth (reads with MAPQ < 30, indels and gaps are filtered)")]))
    header_out.add_format_line(OrderedDict([("ID", "RD"),("Number", "1"), ("Type","Integer"), ("Description", "Reference allele read depth")]))
    header_out.add_format_line(OrderedDict([("ID", "AD"),("Number", "1"), ("Type","Integer"), ("Description", "Alternate allele read depth")]))
    header_out.add_format_line(OrderedDict([("ID", "AF"),("Number", "1"), ("Type","Float"), ("Description", "Allele frequency: AD/(RD+AD). Other alleles, in case of mutli-allelic regions, are ignored.")]))

    # read input vcf
    reader = vcfpy.Reader.from_path(invcf)

    format_ids = reader.header.format_ids()
    for format_id in format_ids:
        format_line = reader.header.get_format_field_info(format_id)
        '''
            output example:
        
            FormatHeaderLine('FORMAT', '<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">', {'ID': 'AD', 'Number': 'R', 'Type': 'Integer', 'Description': 'Allelic depths for the ref and alt alleles in the order listed'})
            key = 'FORMAT'
            value = '<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        '''
        mapping = str_to_mapping(format_line.value)
        mapping["Description"] = "(Info about mutation in the original vcf)" + mapping["Description"]
        header_out.add_info_line(str_to_mapping(format_line.value))

    # open the output vcf
    writer = vcfpy.Writer.from_path(outvcf, header_out) 

    #read bam file
    samfile = pysam.AlignmentFile(bam, "rb")

    #for each mutation in the vcf file
    for record_in in reader:
    # filter out indels: only interested in snvs in this analysis phase
        if not record_in.is_snv():
            continue
        chrom = record_in.CHROM
        pos = record_in.POS-1 #to correct on 1-based positions
        ref = record_in.REF
        alt = record_in.ALT[0].value  #record.ALT is a list by construction which contains only one value
                                    # if the mutation is a SNV
        #line += [call.data.get('GT') or './.' for call in record.calls]

        #look for the pileup in the samfile at position (chrom,pos)
        for pileupcolumn in samfile.pileup(chrom, pos, pos+1, stepper='all', truncate=True, max_depth=10000):
            #number of reads at this position
            sdp = pileupcolumn.n
            #number of supporting reads for the alternate base
            ad = 0
            rd = 0
            dp = 0
            af = 0.0
            for base in pileupcolumn.pileups:
                # .is_del -> the base is a deletion?
                # .is_refskip -> the base is a N in the CIGAR string ?
                if not base.is_del and not base.is_refskip and not base.alignment.mapping_quality < 30:
                    dp += 1
                    if base.alignment.query_sequence[base.query_position] == alt:
                        ad += 1
                    elif base.alignment.query_sequence[base.query_position] == ref:
                        rd += 1

        if ad > 0:
            af = ad / (rd + ad)
            supp = 1
            gt = "0/1" #temporary, all the supported mutations are set to 0/1
        else:
            supp = 0
            gt = "0/0" 


        #af = ad / (rd + ad)

        info_d = {}
        info_d['SUPP'] = supp
        for f in record_in.FORMAT:
            info_d[f] = record_in.calls[0].data.get(f)

        record_out = vcfpy.Record(CHROM=chrom, POS=pos+1, ID=[], REF=ref, ALT=[vcfpy.Substitution(type_="SNV", value=alt)], QUAL=None, FILTER=[], INFO=info_d, FORMAT=["GT","SDP","DP","RD","AD","AF"],
                calls=[vcfpy.Call(sample, OrderedDict([("GT", gt), ("SDP",sdp), ("DP", dp), ("RD", rd), ("AD", ad), ("AF", af)]))]
           )
        writer.write_record(record_out)


    reader.close()
    writer.close()
    samfile.close()

if __name__ == "__main__":
    sys.exit(main())

