import pandas as pd
from pybedtools import BedTool
import argparse
import os

parser = argparse.ArgumentParser(description='Generate annotations:')
parser.add_argument('--pdir', type=str, default=None, help='path to pip files', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix of pip files', required=True)
parser.add_argument('--bedlst', type=str, default=None, help='path to bed list', required=True)
parser.add_argument('--bdir', type=str, default=None, help='path to bed files', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save results', required=True)

args = parser.parse_args()

if not os.path.exists(args.save):
    os.makedirs(args.save)

blst=pd.read_csv(args.bedlst,header=None)
allPIP = pd.DataFrame()
allannot = pd.DataFrame()

for ch in range(1,23):
  pip = pd.read_csv(os.path.join(args.pdir,'{}_{}.pip'.format(args.prefix,ch)),sep='\t',header=None)
  pbed=BedTool([['chr'+i.split('.')[0], int(i.split('.')[1]), int(i.split('.')[1])+1] for i in pip[0]])
  abed=[BedTool(os.path.join(args.bdir,"{}.bed".format(i))).sort() for i in blst[0]]
  ibed=[pbed.intersect(i) for i in abed]
  ist=[[[i.start] for i in j] for j in ibed]
  annot=pd.DataFrame({'SNP':pip[0]})
  annot.index=[int(i.split('.')[1]) for i in pip[0]]
  for ite in range(len(blst)):
    annot['{}'.format(blst.iloc[ite,0])]=0
    annot.loc[annot.index.intersection(pd.DataFrame(ist[ite])[0]),'{}'.format(blst.iloc[ite,0])]=1
  allannot=pd.concat([allannot,annot])
  allPIP=pd.concat([allPIP,pip])

allPIP.to_csv(os.path.join(args.save,'{}.allpip'.format(args.prefix)),sep='\t',header=False,index=False)
allannot.to_csv(os.path.join(args.save,'{}.allannot'.format(args.prefix)),sep='\t',header=True,index=False)
