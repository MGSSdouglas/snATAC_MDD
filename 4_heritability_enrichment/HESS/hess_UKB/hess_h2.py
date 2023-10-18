import pandas as pd
import numpy as np
import scipy.sparse as sparse
import argparse
import os
from scipy.stats import chi2

#obtain from https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/readme_ld.txt
def load_ld_npz(ld_prefix):
    #load the SNPs metadata
    gz_file = '%s.gz'%(ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps['A1'] + '.' + df_ld_snps['A2']
    #load the LD matrix
    npz_file = '%s.npz'%(ld_prefix)
    try:
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s'%(npz_file))
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps

#HESS extended from Shi et al.,2016
def get_HESS_h2_Z(LD,Z,N,LDthres=0.1):
    '''calculate local heritabilities'''
    idx_retain = []
    idx_exclude = [i for i in range(len(Z))]
    zsquare= Z**2
    while len(idx_exclude)>0:
        maxid = idx_exclude[np.argmax(zsquare[idx_exclude])] #find the idx with smallest p-value
        idx_retain.append(maxid)
        idx_exclude = [i for i in idx_exclude if i not in np.where(np.abs(LD[maxid,:])>LDthres)[0]]
    Indidx = np.sort(idx_retain) #obtain independent signals
    P = len(Indidx)
    LD_id = LD[np.ix_(Indidx,Indidx)]
    R_inv = np.linalg.inv(LD_id)
    vec_id = Z[Indidx]
    h2_hess = (np.dot(np.dot(vec_id.transpose(),R_inv),vec_id)-P)/(N-P)
    return h2_hess

parser = argparse.ArgumentParser(description='HESS Commands:')
parser.add_argument('--bed', type=str, default=None, help='path to bed files')
parser.add_argument('--zdir', type=str, default=None, help='path to zscore files')
parser.add_argument('--bp', type=int, default=None, help='half length of peaks')
parser.add_argument('--LDdir', type=str, default=None, help='path to ld files')
parser.add_argument('--N', type=int, default=None, help='sample size')
parser.add_argument('--save', type=str, default=None, help='path to save')
parser.add_argument('--prefix', type=str, default=None, help='prefix')
args = parser.parse_args()

bed = pd.read_csv(args.bed,sep='\t',header=None)
bed['P']=0
bed['h2']=0
z = pd.read_csv(args.zdir,sep='\t',header=None)
z['CHR']=[int(i.split('.')[0]) for i in z[0]]
z['POS']=[int(i.split('.')[1]) for i in z[0]]
z.index=z[0]

for i in range(0,len(bed)):
    ch, st, ed, name = bed.iloc[i,0:4]
    zsub=z.loc[(z['CHR']==int(ch.replace('chr',''))) & (z['POS']>st-args.bp) & (z['POS']<ed+args.bp)]
    if len(zsub)<50 or np.abs(zsub[1]).max()<4.42: #not enough signals or not enough variants
        continue
    ldfile = '{}_{}_{}'.format(ch,max(0,(st-args.bp)//1000000)*1000000+1,max(0,(st-args.bp)//1000000)*1000000+3000001)
    if not os.path.exists("{}.gz".format(os.path.join(args.LDdir,ldfile))):
        continue
    bed.loc[i,'P']=len(zsub)
    df_R, df_ld_snps = load_ld_npz(os.path.join(args.LDdir,ldfile))
    idx = df_R.index.intersection(zsub.index)
    LD = df_R.loc[idx,idx]
    Z = zsub[1].values
    h2_hess=get_HESS_h2_Z(LD.values,Z,args.N)
    print('{} {} {}'.format(name,len(zsub),h2_hess))
    bed.loc[i,'h2']=h2_hess.round(4)
bed.to_csv(os.path.join(args.save,"{}.h2".format(args.prefix)),sep='\t',header=False,index=False)
