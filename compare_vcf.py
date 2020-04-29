import io
import sys
import argparse
import pandas as pd


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def main():

    parser = argparse.ArgumentParser(description="VCF intersection.")

    parser.add_argument("vcf_left", metavar="file1.vcf", help="VCF file", action="store", type=str)
    parser.add_argument("vcf_right", metavar="file2.vcf", help="VCF file", action="store", type=str)
