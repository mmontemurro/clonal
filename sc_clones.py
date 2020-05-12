#!/usr/bin/env python

import gzip
import os
import io
import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import rankdata
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from sklearn.decomposition import PCA
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def read_vcf(path, input_type):

    if input_type == 'vcf':
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
    elif input_type == 'gz':
         with gzip.open(path, 'r') as f:
            lines = [l.decode("utf-8") for l in f if not l.startswith(b'##')]
    else:
        print("Unknown type: {}".format(input_type))
        exit(1)
    return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                  'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def main():

    parser = argparse.ArgumentParser(description="hdbscan")


    parser.add_argument("input1", metavar='clusters.csv', action='store',
            help='Sc-CNV clusters file.', type=str)
    parser.add_argument("input2", metavar='snps.vcf.gz/.vcf', action='store',
            help='Bulk-to-sc SNPs file.', type=str)
    parser.add_argument("prefix", metavar="out/path/prefix", action='store',
            help='Output path prefix.', type=str)
    
    parser.add_argument("--nmuts", metavar='100', action='store', default=100,
                                    help='Number of significantly enriched mutations.', type=int)

    parser.add_argument("--input_type", metavar='vcf/gz', choices=["vcf", "gz"], default="gz", action='store',
            help='Specify if VCF file is compressed (gz) or not (vcf). Default = gz.', type=str)
    parser.add_argument("--low_freq_filter", metavar='0.2', action='store',
            help='Filter out mutations with a bulk AF below the specified threshold.', type=float)
    parser.add_argument("--pvals", metavar='pvals.csv', action='store',
                        help='Pvals file.', type=str)


    args = parser.parse_args()

    input1 = args.input1
    input2 = args.input2
    prefix = args.prefix
    nmuts = args.nmuts
    input_type = args.input_type
    freq_filter_flag = False
    compute_pvals = True
    if args.low_freq_filter:
        freq_threshold = args.low_freq_filter
        freq_filter_flag = True

    if args.pvals:
        pvals_df = pd.read_csv(args.pvals, index_col=0)
        compute_pvals = False

    clusters = pd.read_csv(input1, index_col=0)
    snps = read_vcf(input2, input_type)

    #Filtering bulk mutations not supported by cells
    snps = snps[~snps['INFO'].str.startswith("SUPP=0")]
    
    #Create mutation id column and set it as index
    snps["mutid"] = snps["CHROM"] + "_"+snps["POS"].map(str) + "_" + snps["REF"] + "_" +snps["ALT"]
    snps = snps.set_index('mutid')

    #If AF filter specified, filter out rows where snps["INFO"] specifies an AF < threshold
    #e.g. snps["INFO"] = SUPP=6;GT=0/1;SDP=6;DP=6;RD=0;AD=6;AF=1.0
    if freq_filter_flag:
        snps = snps[pd.Series(snps['INFO'].str.split(";", expand=True)[6]).str.split("=", expand=True)[1].astype(float) >= freq_threshold]
    
    #Filtering cells flagged as outliers in cnv analysis
    snps = snps.drop([c for c in snps.columns if c not in clusters.index.values], axis=1)


    #replace each cell genotype (#GT:DP:RD:AD:AF) with a single binary value
    for c in snps.columns:
        #snps[c] = snps[c].str.split(":", expand=True)[4].astype(float)
        snps[c][snps[c].str.split(":", expand=True)[0] == "0/0"] = 0
    for c in snps.columns:
        #snps[c] = snps[c].str.split(":", expand=True)[4].astype(float)
        snps[c][snps[c].str.split(":", expand=True)[0] == "0/1"] = 1
    snps = snps[snps.columns].astype(int)

    #invert rows and columns
    snps = snps.transpose()
    
    # Prepare contingency matrix 
    # contingency = pd.DataFrame(snps.apply(pd.value_counts).sum(axis=1), columns=['expected'])
    
    #append cluster information and sort by cluster label
    muts = snps.columns
    
    snps['cluster'] = clusters

    snps = snps.sort_values('cluster')
    
    if compute_pvals == True:
        pvals_df = pd.DataFrame(columns=muts, index=clusters['cluster'].unique())

        for c in clusters['cluster'].unique():
            #retrieve rows corresponding to current cluster
            df1 = snps[snps['cluster'] == c].drop(['cluster'], axis=1)
            #compute the observed distribution of each mutation for the current cluster
            observed_fract = df1.apply(pd.value_counts).fillna(0)
            
            #all other clusters
            df2 = snps[snps['cluster'] != c].drop(['cluster'], axis=1)
            #compute the observed distribution of each mutation for all other clusters
            expected_fract = df2.apply(pd.value_counts).fillna(0)
            for m in muts:
                contingency = [observed_fract[m].tolist(),expected_fract[m].tolist()]
                _, pvalue = fisher_exact(contingency)
                #print(pvalue)
                pvals_df.loc[c,m] = pvalue
                
        pvals_df.to_csv(prefix+"_fisher.pvals.csv")
        """
        # Transform pvals df to list to perform Benjamini/Hochberg (non-negative) pvalue adjustment
        pvals_list = []
        for idx, row in pvals_df.iterrows():
            pvals_list = np.append(pvals_list, row)
        #print("Pvals pre-adjustment: ")
        #print(pvals_list)

        
        reject, adj_pvals_list, _, _ = multipletests(pvals_list, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        print("Adjusted pvals: ")
        print(adj_pvals_list)
        # generate multiindex to map adjusted pvals to cluster and mutation
        # c1 mut1 pval
        #    mut2 pval
        #    mut3 pval
        # c2 mut1 pval
        #    mut2 pval
        #    mut3 pval
        # ...
        rows = pvals_df.index.values
        columns = pvals_df.columns.values
        
        l1 = []
        for r in rows:
            l1 = np.append(l1, [r] * len(columns))

        l2 = []
        for i in range(len(rows)):
            l2 = np.append(l2, columns)

        levels = [l1, l2] #MultiIndex
        # Generate series and unstack to a dataframe with clusters as idx and muts as columns
        
        adj_pvals_df = pd.Series(adj_pvals_list, index=levels).unstack(level=1) 
        #print(adj_pvals_df)
        adj_pvals_df.to_csv(prefix+"_fisher.pvals.adj.csv")
        
        """

    # Get minimum pval for each mutation, sort them in descending order and save this list into a new series
    pvals_df = pvals_df[pvals_df.columns].astype(float)
    
    ranked_min_pvals = pvals_df.min().sort_values()
    ranked_min_pvals.to_csv(prefix+"_fisher.pvals.no_correct.rank.csv")
    #min_pvals_df['cluster'] = pvals_df.idxmin()
    
    top = ranked_min_pvals.head(nmuts).index.values
    
    snps = snps.set_index('cluster', append=True).swaplevel(0,1)
    
    
    snps = snps.drop([c for c in snps.columns.values if c not in top], axis=1)

    #keep only mutations in top1000
    
    #Plot heatmap
    labels = snps.index.get_level_values('cluster')
    color_palette = sns.color_palette("hls", len(np.unique(labels)))
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in labels]
    cmap = colors.ListedColormap(['whitesmoke', 'k'])
    
    h = sns.clustermap(snps,
            row_cluster=False,
            col_cluster=False,
            yticklabels = False,
            row_colors=cluster_colors,
            cmap=cmap,
            cbar_kws={"drawedges":True, })

    colorbar = h.ax_heatmap.collections[0].colorbar
    colorbar.solids.set_edgecolor("black")
    colorbar.set_ticks([0.25,0.75])
    colorbar.set_ticklabels(['0', '1'])

    plt.setp(h.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)
    #h.cax.set_position([0.05, .2, .03, .45])
    plt.gcf().suptitle('sc-CNV clusters VS bulk SNPs', fontsize=20, fontweight='bold' )
    plt.gcf().set_size_inches(12, 21)
    plt.savefig(prefix+'_snvs_to_cnvs.heatmap.png')
    plt.clf()
    
    
    #cluster on mutations
    
    h = sns.clustermap(snps,
            row_cluster=False,
            col_cluster=True,
            yticklabels = False,
            row_colors=cluster_colors,
            cmap=cmap)

    colorbar = h.ax_heatmap.collections[0].colorbar
    colorbar.solids.set_edgecolor("black")
    colorbar.set_ticks([0.25,0.75])
    colorbar.set_ticklabels(['0', '1'])

    plt.setp(h.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)
    h.cax.set_position([0.05, .2, .03, .45])
    plt.gcf().suptitle('sc-CNV clusters VS bulk SNPs', fontsize=20, fontweight='bold' )
    plt.gcf().set_size_inches(12, 21)
    plt.savefig(prefix+'_snvs_to_cnvs.heatmap.clustered.png')
    plt.clf()
    
    """
    # Filter out mutations which mean AF over cells < 0.6
    #df.loc[:,df.max() > 0]
    #snps = snps.loc[:, snps.mean() >= 0.6]
    
    q_90th = snps.mean().quantile(0.90)
    q_75th = snps.mean().quantile(0.75)
    median = snps.mean().quantile(0.5)
    sns.distplot(snps.mean(), rug=True)
    plt.gcf().suptitle("Mean single-cell AF distribution")
    plt.axvline(x=median, color="red", linestyle="--", label="median")
    plt.axvline(x=q_75th, color="blue", linestyle="--", label="75th quantile")
    plt.axvline(x=q_90th, color="green", linestyle="--", label="90th quantile")
    plt.legend()
    plt.gcf().set_size_inches(20, 12)
    plt.savefig(prefix+".mean_af.distr.png")
    plt.clf()

    
    snps = snps.loc[:, snps.mean() >= q_90th]
    print(snps)
    h = sns.clustermap(snps,
            row_cluster=False,
            col_cluster=False,
            yticklabels = False,
            row_colors=cluster_colors,
            cmap='cividis')

    plt.setp(h.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)
    h.cax.set_position([0.05, .2, .03, .45])
    plt.gcf().suptitle('sc-CNV clusters VS bulk SNPs (mean AF > 90th percentile)', fontsize=20, fontweight='bold' )
    plt.gcf().set_size_inches(37, 21)
    plt.savefig(prefix+'snvs_to_cnvs.heatmap.90_perc.png')
    plt.clf()

    h = sns.clustermap(snps,
            row_cluster=True,
            col_cluster=True,
            yticklabels = False,
            row_colors=cluster_colors,
            cmap='cividis')

    plt.setp(h.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)
    h.cax.set_position([0.05, .2, .03, .45])
    plt.gcf().suptitle('sc-CNV clusters VS bulk SNPs (mean AF > 90th percentile)', fontsize=20, fontweight='bold' )
    plt.gcf().set_size_inches(37, 21)
    plt.savefig(prefix+'snvs_to_cnvs.heatmap.90_perc.clusters.png')
    plt.clf()
    """


if __name__ == "__main__":
    sys.exit(main())

