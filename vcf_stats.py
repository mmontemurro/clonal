from funcs import read_vcf

import sys
import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="vcf stats plotter")


parser.add_argument("input", metavar='file.vcf', action='store',
                help='vcf file (without header).', type=str)

parser.add_argument("sample", metavar='sample1', action='store',
                        help='Sample name.', type=str)
parser.add_argument("outdir", metavar='path/to/output/dir', action='store',
                        help='Filepath of the desidered output directory.', type=str)


args = parser.parse_args()

vcf = args.input
outdir = args.outdir
sample = args.sample

df = read_vcf(vcf)

columns = df['FORMAT'][0].split(":") #expand format string
df_res = df[sample].str.split(":", expand=True)
df_res.columns = columns

df_res.FREQ = pd.to_numeric(df_res.FREQ.str[:-1]) #remove '%'
df_res.AD = pd.to_numeric(df_res.AD)
df_res.RD = pd.to_numeric(df_res.RD)

df_res["read_count"] = df_res.AD + df_res.RD

df_01 = df_res[df_res["GT"] == "0/1"]
df_11 = df_res[df_res["GT"] == "1/1"]

#density plots
ax = sns.distplot(df_res["FREQ"])
ax.set_xlabel("Variant allele frequency")
plt.gcf().suptitle("Variant allele frequency density plot")
plt.savefig(outdir+"/af_density.png")
plt.clf()

ax = sns.distplot(df_01["FREQ"])
ax.set_xlabel("Variant allele frequency (0/1)")
plt.gcf().suptitle("Variant allele frequency density plot (0/1 variants)")
plt.savefig(outdir+"/af_density_01.png")
plt.clf()

ax = sns.distplot(df_11["FREQ"])
ax.set_xlabel("Variant allele frequency (1/1)")
plt.gcf().suptitle("Variant allele frequency density plot (1/1 variants)")
plt.savefig(outdir+"/af_density_11.png")
plt.clf()

#scatter plots af vs variant allele read count
ax = df_01.plot.scatter(x="FREQ", y="AD")
ax.set_xlabel("Variant allele frequency (0/1)")
ax.set_ylabel("Variant allele read count")
plt.gcf().suptitle("Variant allele frequency vs variant allele read count")
plt.gcf().set_size_inches(20,20)
plt.savefig(outdir+"/variant_read_count_scatter_01.png")
plt.clf()

ax = df_11.plot.scatter(x="FREQ", y="AD")
ax.set_xlabel("Variant allele frequency (1/1)")
ax.set_ylabel("Variant allele read count")
plt.gcf().suptitle("Variant allele frequency vs variant allele read count")
plt.gcf().set_size_inches(20,20)
plt.savefig(outdir+"/variant_read_count_scatter_11.png")
plt.clf()

#scatter plots af vs tot read count
ax = df_01.plot.scatter(x="FREQ", y="read_count")
ax.set_xlabel("Variant allele frequency (0/1)")
ax.set_ylabel("Tot read count")
plt.gcf().suptitle("Variant allele frequency vs tot read count")
plt.gcf().set_size_inches(20,20)
plt.savefig(outdir+"/tot_read_count_scatter_01.png")
plt.clf()

ax = df_11.plot.scatter(x="FREQ", y="AD")
ax.set_xlabel("Variant allele frequency (1/1)")
ax.set_ylabel("Tot read count")
plt.gcf().suptitle("Variant allele frequency vs tot read count")
plt.gcf().set_size_inches(20,20)
plt.savefig(outdir+"/tot_read_count_scatter_11.png")
plt.clf()








