#!/usr/bin/env python
from datetime import date
import sys
import argparse
import vcfpy
from collections import OrderedDict
import pandas as pd

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


def main():
    parser = argparse.ArgumentParser(description="From single cell VCF to clones vcf.")
    parser.add_argument("input1", metavar="sample.muts.vcf", action="store", help="Single cell VCF file.", type=str)
    parser.add_argument("input2", metavar="clusters.list", action="store", help="Clusters list.", type=str)
    #parser.add_argument("input_type", choices=["gz", "vcf"], help="VCF input type (vcf/gz).", type=str)
    #parser.add_argument("sample", metavar="sample_name", action="store", help="Sample name", type=str)
    parser.add_argument("outprefix", metavar="out/path/prefix", action="store", help="Output prefix", type=str)

    args = parser.parse_args()

    input1 = args.input1
    input2 = args.input2
    prefix = args.outprefix
    #sample = args.sample
    #input_type = args.input_type

        
    clusters_df = pd.read_csv(input2)
    #clusters_df['cluster'] = clusters_df['a'].apply(lambda x: "{}_{}".format(sample, x))    

    clusters = [str(cluster) for cluster in clusters_df['cluster'].unique()]
    # Create out header
    header_out = vcfpy.Header(lines=[ vcfpy.HeaderLine(key="fileformat", value="VCFv4.3"), vcfpy.HeaderLine(key="source", value=sys.argv[0]), vcfpy.HeaderLine(key="fileDate", value=date.today().strftime("%d/%m/%Y")) ], samples=vcfpy.SamplesInfos(clusters))
     
    # format header lines 
    header_out.add_format_line(OrderedDict([("ID", "GT"),("Number", "1"), ("Type","String"), ("Description", "Genotype (0/1, 0/0)")]))
    header_out.add_format_line(OrderedDict([("ID", "DP"),("Number", "1"), ("Type","Integer"), ("Description", "Filtered read depth (reads with MAPQ < 30, indels and gaps are filtered)")]))
    header_out.add_format_line(OrderedDict([("ID", "RD"),("Number", "1"), ("Type","Integer"), ("Description", "Reference allele read depth")]))
    header_out.add_format_line(OrderedDict([("ID", "AD"),("Number", "1"), ("Type","Integer"), ("Description", "Alternate allele read depth")]))
    header_out.add_format_line(OrderedDict([("ID", "AF"),("Number", "1"), ("Type","Float"), ("Description", "Allele frequency: AD/(RD+AD). Other alleles, in case of mutli-allelic regions, are ignored.")]))
    
    # info header lines
 
    header_out.add_info_line(OrderedDict([("ID", "SUPP"), ("Number", "1"), ("Type","Integer"), ("Description", "Whether the mutation is supported or not.")]))
    
    # read input vcf
    reader = vcfpy.Reader.from_path(input1)
    # open the output vcf
    writer = vcfpy.Writer.from_path(prefix+"_clusters.vcf", header_out)
 
    """
    snps = read_vcf(input1, input_type)
    #Filtering bulk mutations not supported by cells
    snps = snps[~snps['INFO'].str.startswith("SUPP=0")]
    
    #Create mutation id column and set it as index
    snps["mutid"] = snps["CHROM"] + "_"+snps["POS"].map(str) + "_" + snps["REF"] + "_" +snps["ALT"]
    snps = snps.set_index('mutid')
    """

    #for each record in the vcf file
    for record_in in reader:
        d = samples_dict(clusters_df['cluster'].unique())
        supp = 0
        chrom = record_in.CHROM
        pos = record_in.POS-1 #to correct on 1-based positions
        ref = record_in.REF
        alt = record_in.ALT[0].value
        
        #for each cluster compute 'GT:DP:RD:AD:AF' to be provided as call argument
        for c in clusters_df['cluster'].unique():
            #retrieve cell columns for cells in current cluster
            cells = clusters_df['cellid'][clusters_df['cluster'] == c]
            
          
            #retrieve cell data
            calls = [record_in.call_for_sample[cell] for cell in cells]
            #sum total read count, alt read count and ref read count of cells in the cluster
            for call in calls:    
                d[c]['dp'] = d[c]['dp'] + call.data.get('DP') 
                d[c]['rd'] = d[c]['rd'] + call.data.get('RD')
                d[c]['ad'] = d[c]['dp'] + call.data.get('AD')

            if d[c]['ad'] > 0:
                d[c]['gt'] = "0/1"
                d[c]['af'] = d[c]['ad'] / (d[c]['rd'] + d[c]['ad'])
                supp = 1
    
        calls = []
        # create one call for each cluster
        for c in d.keys():
            calls.append(vcfpy.Call(str(c), OrderedDict([("GT", d[c]['gt']), ("DP", d[c]['dp']), ("RD", d[c]['rd']), ("AD", d[c]['ad']), ("AF", d[c]['af'])])))        
        print(calls)
         
        # write new record
        record_out = vcfpy.Record(CHROM=chrom, POS=pos+1, ID=[], REF=ref, ALT=[vcfpy.Substitution(type_="SNV", value=alt)], QUAL=None, FILTER=[], INFO={"SUPP":supp}, FORMAT=["GT","DP","RD","AD","AF"],
                calls=calls
           )
        writer.write_record(record_out)
        
    reader.close()
    writer.close()

if __name__ == "__main__":
    sys.exit(main())
