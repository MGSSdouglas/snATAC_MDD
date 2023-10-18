import pandas as pd
from pybedtools import BedTool
from pybedtools.cbedtools import Attributes
import argparse

parser = argparse.ArgumentParser(description='Map to nearest gene:')
parser.add_argument('--gtf', type=str, default=None, help='gtf file contain gene info', required=True)
parser.add_argument('--rs', type=str, default=None, help='dictionary for rsid conversion', required=True)
parser.add_argument('--cs', type=str, default=None, help='path to sharepro results', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save results', required=True)
parser.add_argument('--eQTL', type=str, default=None, help='eQTL dictionary', required=True)
parser.add_argument('--sQTL', type=str, default=None, help='sQTL dictionary', required=True)
args = parser.parse_args()

genebed = BedTool(args.gtf).sort()
rsdict = pd.read_csv(args.rs, sep='\t', index_col=0, header=None)[1]

cs = pd.concat([pd.read_csv('{}{}.cs'.format(args.cs, i), sep='\t') for i in range(1, 23)])
cs['rsid'] = [[rsdict.get(i) for i in j.split('/')] for j in cs['cs']]
sqtl = pd.read_csv(args.sQTL, sep='\t', header=None, index_col=0)[1].sort_index()
eqtl = pd.read_csv(args.eQTL, sep='\t', header=None, index_col=0)[1].sort_index()
cs['sQTL'] = ['/'.join(set([elem for i in j for elem in sqtl.get(i, 'NA').split(';') if elem != 'NA']))
              for j in cs['rsid']]
cs['eQTL'] = ['/'.join(set([elem for i in j for elem in eqtl.get(i, 'NA').split(';') if elem != 'NA']))
              for j in cs['rsid']]


def get_nearest_gene(cs, genebed):
    allchr = [[int(j.split('.')[0]) for j in i.split('/')][0] for i in cs['cs']]
    allpos = [[int(j.split('.')[1]) for j in i.split('/')] for i in cs['cs']]
    cs['Top_variant'] = [i.split('/')[0] for i in cs['cs']]
    cs['Chr'] = ['chr' + str(i) for i in allchr]
    cs['Start'] = [min(i) for i in allpos]
    cs['End'] = [max(i) for i in allpos]
    csbed = BedTool(cs[['Chr', 'Start', 'End']].values.tolist()).sort()
    cgene = csbed.closest(genebed, d=True)
    c2g = {'.'.join(i[:3]): Attributes(i[11])['gene_name'] for i in cgene}
    c2d = {'.'.join(i[:3]): int(i[12]) for i in cgene}
    cs.index = cs['Chr'] + '.' + cs['Start'].astype('str') + '.' + cs['End'].astype(str)
    cs.index.name = 'Region'
    cs['Nearest_gene'] = pd.Series(c2g)
    cs['Distance_to_nearest_gene'] = pd.Series(c2d)
    return cs


cs = get_nearest_gene(cs, genebed)
cs['rsid'] = ['/'.join(i) for i in cs['rsid']]
cs[['Top_variant', 'Nearest_gene', 'Distance_to_nearest_gene', 'eQTL', 'sQTL', 'cs',
    'rsid']].to_csv(args.save, sep='\t', header=True, index=False)
