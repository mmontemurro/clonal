#!/usr/bin/env python
import sys
import shlex
import vcfpy
import pysam
import argparse
from datetime import date
from collections import OrderedDict

def samples_dict(samples):
    d = {}
    for s in samples:
        d[s] = {"gt" : "0/0",
                "dp" : 0,
                "rd" : 0,
                "ad" : 0,
                "af" : 0.0
                }
    return d

def list_to_dict(list_):
    dict_ = {}
    for item in list_:
        dict_[item[0]] = item[1]
    return dict_

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

    parser.add_argument("barcodes", metavar='barcodes.list', action='store',
            help="File containing cell barcodes (the same used in the alignment file to identify cell reads).", type=str)

    parser.add_argument("vcf", metavar='variants.vcf', action='store',
                    help="VCF file storing BULK SNPs.", type=str)

    parser.add_argument("sample_name", metavar='sample1', action='store',
                    help="Sample identifier.", type=str)


    parser.add_argument("out_prefix", metavar="outdir/sample", action="store",
            help="Output VCF file prefix.", type=str)


    parser.add_argument("--gt", metavar='1/1 (0/1)', choices=["0/0", "0/1", "1/1"], action='store',
            help="Genotype filter: considers only mutations with the specified GT in the original vcf file.", type=str)

    args = parser.parse_args()
    bam= args.bam
    barcodes = args.barcodes
    invcf = args.vcf
    sample = args.sample_name
    outvcf = args.out_prefix + ".snpseeker.vcf"

    if args.gt:
        gt_filter = True
        gt = args.gt
    
    else:
        gt_filter = False

    with open(barcodes, "r") as f:
        samples = f.read().splitlines()
    #read bam file
    samfile = pysam.AlignmentFile(bam, "rb")

    #build the header of the output vcf
    header_out = vcfpy.Header(lines=[vcfpy.HeaderLine(key="fileformat", value="VCFv4.3"), vcfpy.HeaderLine(key="source", value=sys.argv[0]), vcfpy.HeaderLine(key="fileDate", value=date.today().strftime("%d/%m/%Y")) ], samples=vcfpy.SamplesInfos(samples))

    # sample header lines
    header_out.add_line(vcfpy.SampleHeaderLine.from_mapping(OrderedDict([("ID", sample),("Description", "Sample name")])))
    
    # filter header lines
    # sample header lines
    header_out.add_filter_line(OrderedDict([("ID", "1/1"),("Number", "1"), ("Description", "Filtered on such GT")]))
    header_out.add_filter_line(OrderedDict([("ID", "0/1"),("Number", "1"), ("Description", "Filtered on such GT")]))
    header_out.add_filter_line(OrderedDict([("ID", "0/0"),("Number", "1"), ("Description", "Filtered on such GT")]))


    #header_out.add_info_line(OrderedDict([("ID", "MUT"), ("Number", "1"), ("Type","Integer"), ("Description", "States if the record mutation is supported (1) or not (0).")]))
     
    # format header lines 
    header_out.add_format_line(OrderedDict([("ID", "GT"),("Number", "1"), ("Type","String"), ("Description", "Genotype (0/1, 0/0)")]))
    header_out.add_format_line(OrderedDict([("ID", "DP"),("Number", "1"), ("Type","Integer"), ("Description", "Filtered read depth (reads with MAPQ < 30, indels and gaps are filtered)")]))
    header_out.add_format_line(OrderedDict([("ID", "RD"),("Number", "1"), ("Type","Integer"), ("Description", "Reference allele read depth")]))
    header_out.add_format_line(OrderedDict([("ID", "AD"),("Number", "1"), ("Type","Integer"), ("Description", "Alternate allele read depth")]))
    header_out.add_format_line(OrderedDict([("ID", "AF"),("Number", "1"), ("Type","Float"), ("Description", "Allele frequency: AD/(RD+AD). Other alleles, in case of mutli-allelic regions, are ignored.")]))

    # read input vcf
    reader = vcfpy.Reader.from_path(invcf)

    # info header lines
    # Use input FORMAT lines as output INFO line 
    header_out.add_info_line(OrderedDict([("ID", "SUPP"), ("Number", "1"), ("Type","Integer"), ("Description", "Number of cells supporting the mutation.")]))
    
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
        mapping["Description"] = "(Info about bulk mutation)" + mapping["Description"]
        header_out.add_info_line(str_to_mapping(format_line.value))



    # open the output vcf
    writer = vcfpy.Writer.from_path(outvcf, header_out) 

    #read bam file
    samfile = pysam.AlignmentFile(bam, "rb")
    
    
    #for each mutation in the vcf file
    for record_in in reader:
        d = samples_dict(samples)    
        supp = 0
        # filter out indels: only interested in snvs in this analysis phase
        if gt_filter:
            if record.calls[0].data.get('GT') != gt:
                continue

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
            for base in pileupcolumn.pileups:
                # .is_del -> the base is a deletion?
                # .is_refskip -> the base is a N in the CIGAR string ?
                if not base.is_del and not base.is_refskip and not base.alignment.mapping_quality < 30:
                    #iterate on cells
                    tags = list_to_dict(base.alignment.tags)
                    if "CB" not in tags.keys():
                        ''' reads with no error-corrected barcode are discarded '''
                        continue
                    elif tags["CB"].split("-")[0] not in samples:
                        ''' The barcode hasn't been labeled has belonging to a cell by cellranger (floating DNA)'''
                        continue
                    cb = tags["CB"].split("-")[0] #10x barcodes
                    #print("barcode {} is a cell barcode ".format(cb))
                    d[cb]['dp'] += 1 #update info for the sample identified by CB
                    if base.alignment.query_sequence[base.query_position] == alt:
                        d[cb]['ad'] += 1
                    elif base.alignment.query_sequence[base.query_position] == ref:
                        d[cb]['rd'] += 1
        for cb in d.keys():
            if d[cb]['ad'] > 0:
                supp += 1
                d[cb]['gt'] = "0/1" #temporary, all the supported mutations are set to 0/1
                d[cb]['af'] =  d[cb]['ad'] / (d[cb]['rd'] + d[cb]['ad'])
     

        # generate calls for each sample/cell
        calls = []
        for cb in d.keys():
            calls.append(vcfpy.Call(cb, OrderedDict([("GT", d[cb]['gt']), ("DP", d[cb]['dp']), ("RD", d[cb]['rd']), ("AD", d[cb]['ad']), ("AF", d[cb]['af'])])))        
        
        
        # create a mapping between each FORMAT entry and the 
        # corresponding value, in the call, in the input vcf file
        # note that the input vcf contains only one sample, so
        # the calls field of each record contains only one entry
        info_d = {}
        info_d['SUPP'] = supp
        for f in record_in.FORMAT:
            info_d[f] = record_in.calls[0].data.get(f)
        
        if gt_filter == True:
            filter_l = [gt]
        else:
            filter_l = []

        # build and write the output record
        
        record_out = vcfpy.Record(CHROM=chrom, POS=pos+1, ID=[], REF=ref, ALT=[vcfpy.Substitution(type_="SNV", value=alt)], QUAL=None, FILTER=filter_l, INFO=info_d, FORMAT=["GT","DP","RD","AD","AF"],
                calls=calls
           )
        writer.write_record(record_out)


    reader.close()
    writer.close()
    samfile.close()


if __name__ == "__main__":
    sys.exit(main())

