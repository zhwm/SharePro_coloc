import pandas as pd
import argparse
import os
from scipy.stats import chi2

parser = argparse.ArgumentParser(description='Match GWAS and QTL summary statistics with bim file using rsid')
parser.add_argument('--rss', type=str, default=None, nargs='+', help='path to GWAS summary stats', required=True)
parser.add_argument('--bim', type=str, default=None, help='path to bim file', required=True)
parser.add_argument('--rsID', type=str, default=None, nargs='+', help='column names for rsid in rss', required=True)
parser.add_argument('--A1', type=str, default=None, nargs='+', help='column names for A1', required=True)
parser.add_argument('--A2', type=str, default=None, nargs='+', help='column names for A2', required=True)
parser.add_argument('--BETA', type=str, default=None, nargs='+', help='column names for effect size', required=True)
parser.add_argument('--SE', type=str, default=None, nargs='+',
                    help='column names for effect size standard error', required=True)
parser.add_argument('--N', type=str, default=None, nargs='+', help='column names for sample size')
parser.add_argument('--EAF', type=str, default=None, nargs='+', help='column names for allele frequency')
parser.add_argument('--save', type=str, default=None, help='path to save output', required=True)
parser.add_argument('--prefix', type=str, default=None, nargs='+', help='prefix of output', required=True)
parser.add_argument('--cols', type=str, default=None, nargs='+',
                    help='columns in output: SNP/A1/A2/CHR/POS/BETA/SE/P/Z', required=True)
parser.add_argument('--header', action="store_true", help='options for adding header')

args = parser.parse_args()
if not os.path.exists(args.save):
    os.makedirs(args.save)
assert len(args.rss) == len(args.rsID) == len(args.A1) == len(args.A2) == len(args.BETA) == len(args.SE), \
  'Please make sure inputs have the same length'
nums = len(args.rss)
rss = [pd.read_csv(args.rss[i], sep='\s+', index_col=args.rsID[i]) for i in range(nums)]
bim = pd.read_csv(args.bim, sep='\s+', index_col=1, header=None)
mid = bim.index.intersection(set.intersection(*[set(i.index.drop_duplicates(keep=False)) for i in rss]))
print('{} SNPs matched'.format(len(mid)))
matchbim = bim.loc[mid]
ss = [i.loc[mid] for i in rss]
for i in range(nums):
    ss[i]['SNP'] = ss[i].index
    ss[i]['A1'] = ss[i][args.A1[i]].str.upper()
    ss[i]['A2'] = ss[i][args.A2[i]].str.upper()
    ss[i]['sign'] = 1
    ss[i].loc[ss[i]['A1'] == matchbim[4], 'sign'] = -1
    ss[i]['A1'] = matchbim[5]
    ss[i]['A2'] = matchbim[4]
    ss[i]['CHR'] = matchbim[0]
    ss[i]['BETA'] = ss[i][args.BETA[i]]*ss[i]['sign']
    ss[i]['SE'] = ss[i][args.SE[i]]
    ss[i]['Z'] = (ss[i]['BETA'])/(ss[i]['SE'])
    ss[i]['P'] = chi2.sf((ss[i]['Z'])**2, 1)
    ss[i]['POS'] = matchbim[3]
    if args.EAF is not None:
        ss[i]['EAF'] = ss[i][args.EAF[i]]
    if args.N is not None:
        ss[i]['N'] = ss[i][args.N[i]]
    ss[i][args.cols].to_csv(os.path.join(args.save, "{}.txt".format(args.prefix[i])),
                            sep='\t', header=args.header, index=False)
