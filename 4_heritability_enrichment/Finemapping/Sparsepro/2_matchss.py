import pandas as pd
import argparse
import os
import sys
import numpy as np

def preprocess(ss,args):
    ss = ss.drop_duplicates()
    if args.A1 is not None and args.A2 is not None:
        ss[args.A1] = ss[args.A1].str.upper()
        ss[args.A2] = ss[args.A2].str.upper()
        ss.rename(columns={args.A1:'A1',args.A2:'A2'},inplace=True)
    else:
        sys.exit("Please provide A1/A2 information or use --ukb")
    if args.log:
        ss[args.BETA]=np.log(ss[args.BETA])
    if args.rsid is not None and args.split:
        ss[args.rsid]=[i.split(':')[0] for i in ss[args.rsid]] # for GIANT
    if args.Z is not None:
        ss.rename(columns={args.Z:'Z'},inplace=True)
    elif args.BETA is not None and args.SE is not None:
        ss['Z']=ss[args.BETA]/ss[args.SE]
    else:
        sys.exit("Please provide either Z score or BETA|SE")
    if args.rsid is not None:
        print("Matching with rsid")
        ss.rename(columns={args.rsid:'rsid'},inplace=True)
        ss=ss.drop_duplicates(subset='rsid',keep=False)
    elif args.CHR is not None and args.POS is not None:
        print("Matching with CHR|POS")
        ss.rename(columns={args.CHR:'CHR',args.POS:'POS'},inplace=True)
        ss=ss.drop_duplicates(subset=['CHR','POS'],keep=False)
        ss['CHR']=[int(i) for i in ss['CHR']]
        ss['POS']=[int(i) for i in ss['POS']]
    return ss

def match(ss,args,ch):
    idx=pd.read_csv(os.path.join(args.idir,"{}.rsid".format(ch)),header=None,sep='\t')
    idx.columns=['SNP','rsid']
    idx['POS']=[int(i.split('.')[1]) for i in idx['SNP']]
    idx['A1']=[i.split('.')[2] for i in idx['SNP']] 
    idx['A2']=[i.split('.')[3] for i in idx['SNP']]
    if args.rsid is not None:
        sub=ss.set_index('rsid').loc[ss['rsid'][ss['rsid'].isin(idx['rsid'])]]
        idx=idx.set_index('rsid').loc[sub.index]
    elif args.CHR is not None and args.POS is not None:
        sub=ss.loc[(ss['CHR']==ch)]
        sub=sub.set_index('POS').loc[sub['POS'][sub['POS'].isin(idx['POS'])]]
        idx=idx.set_index('POS').loc[sub.index]
    ol=len(sub)
    sub=sub.loc[((sub['A2']==idx['A2']) & (sub['A1']==idx['A1'])) | ((sub['A2']==idx['A1']) & (sub['A1']==idx['A2']))]
    idx=idx.loc[sub.index]
    idx['sign']=1
    idx.loc[((sub['A2']==idx['A2']) & (sub['A1']==idx['A1'])),'sign']=-1
    sub['Z']=(sub['Z']*idx['sign']).round(4)
    sub=sub.reset_index()
    idx=idx.reset_index()
    sub['SNP']=idx['SNP'].values
    sub['rsid']=idx['rsid'].values
    sub['CHR']=ch
    sub['POS']=idx['POS'].values
    sub['A2']=idx['A1'].values
    sub['A1']=idx['A2'].values
    sub.index=sub['SNP']
    print("matching {} out of {} variants in chromosome {}".format(len(sub), ol, ch))
    return sub

def preprocess_ukb(ss,args):
    sub=ss.loc[ss.low_confidence_variant==False].copy()
    sub['SNP']=[i.replace(':','.') for i in sub.variant]
    sub.index = sub['SNP']
    return sub

def match_with_ukb(ss,args,ch):
    idx=pd.read_csv(os.path.join(args.idir,"{}.rsid".format(ch)),header=None,sep='\t')
    sub = ss.loc[ss.index.intersection(idx[0]),['SNP','tstat']].copy()
    sub.columns = ['SNP','Z']
    print("matching {} variants in chromosome {}".format(len(sub), ch))
    return sub

def merge(args):
    SS = [pd.read_csv(i,sep=args.delim,compression=args.compression) for i in args.rss]
    if not args.ukb:
        SS = [preprocess(i,args) for i in SS]
        for ch in args.ch:
            matched=[match(i,args,ch) for i in SS]
            mid = pd.concat([i['SNP'] for i in matched],axis=1,join='inner').index
            if len(mid)==0:
                continue
            [matched[i].loc[mid,args.col].to_csv(os.path.join(args.save,"{}_{}.z".format(args.prefix[i],ch)),sep='\t',header=False,index=False) for i in range(len(matched))]
    else:
        SS = [preprocess_ukb(i,args) for i in SS]
        for ch in args.ch:
            matched=[match_with_ukb(i,args,ch) for i in SS]
            mid = pd.concat([i['SNP'] for i in matched],axis=1,join='inner').index
            if len(mid)==0:
                continue
            [matched[i].loc[mid].to_csv(os.path.join(args.save,"{}_{}.z".format(args.prefix[i],ch)),sep='\t',header=False,index=False) for i in range(len(matched))]

parser = argparse.ArgumentParser(description='Match Summary Statistics:')
parser.add_argument('--rss', type=str, default=None, nargs='+', help='path to raw summary stats', required=True)
parser.add_argument('--prefix',type=str, default=None, nargs='+', help='prefix of output', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save output', required=True)
parser.add_argument('--idir', type=str, default=None, help='path to index directory',required=True)
parser.add_argument('--col', type=str, default=['SNP','Z'], nargs='+', help='columns in output')
parser.add_argument('--rsid', type=str, default=None, help='name of the rsid column')
parser.add_argument('--CHR', type=str, default=None, help='name of the chromosome column')
parser.add_argument('--POS', type=str, default=None, help='name of the position column')
parser.add_argument('--A1', type=str, default=None, help='name of the A1 column')
parser.add_argument('--A2', type=str, default=None, help='name of the A2 column')
parser.add_argument('--BETA', type=str, default=None, help='name of the effect column')
parser.add_argument('--SE', type=str, default=None, help='name of standard error of the effect column')
parser.add_argument('--Z', type=str, default=None, help='name of the standardized effect Z column')
parser.add_argument('--compression',type=str, default=None, help='type of compression in GWAS summary statistics')
parser.add_argument('--delim',type=str, default='\t', help='type of delimiter in GWAS summary statistics')
parser.add_argument('--log', action='store_true', help='options for providing OR as effect')
parser.add_argument('--split', action='store_true', help='option for spliting rsid for GIANT')
parser.add_argument('--ukb', action='store_true', help='option for processing UK Biobank summary statistics from Neale lab')
parser.add_argument('--ch', type=int, default=None, nargs='+',help='which chromosome should be included')

args = parser.parse_args()

if not os.path.exists(args.save):
    os.makedirs(args.save)

merge(args)
