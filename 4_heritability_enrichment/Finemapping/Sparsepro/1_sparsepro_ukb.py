import pandas as pd
import argparse
import time
import os
import numpy as np
from scipy.special import softmax
from scipy.stats import chi2, entropy
import scipy.sparse as sparse
import sys


def title():
    print('**********************************************************************')
    print('* SparsePro for efficient genome-wide fine-mapping                   *')
    print('* Version 3.0.0                                                      *')
    print('* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *')
    print('**********************************************************************')
    print()


# obtained from https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/readme_ld.txt
def load_ld_npz(ld_prefix):
    # load the SNPs metadata
    gz_file = '%s.gz' % (ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid': 'SNP', 'chromosome': 'CHR', 'position': 'BP', 'allele1': 'A1', 'allele2': 'A2'},
                      inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps[
        'A1'] + '.' + df_ld_snps['A2']
    # load the LD matrix
    npz_file = '%s.npz' % (ld_prefix)
    try:
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s' % (npz_file))
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps


class SparsePro(object):

    def __init__(self, P, K, XX, h2, var_b):
        """initialize and set hyperparameters"""
        self.p = P
        self.k = K
        self.gamma = np.zeros((self.p, self.k))
        self.beta_mu = np.zeros((self.p, self.k))
        self.beta_prior_tau = np.tile((1.0 / var_b * np.array([k + 1 for k in range(self.k)])), (self.p, 1))
        self.y_tau = 1.0 / (1.0 - h2)  # var_Y get cancelled eventually
        self.prior_pi = np.ones((self.p,)) * (1 / self.p)
        self.beta_post_tau = np.tile(XX.reshape(-1, 1), (1, self.k)) * self.y_tau + self.beta_prior_tau

    def infer_q_beta(self, ytX, XtX):
        """perform variational updates"""
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            beta_all_k = (self.gamma[:, idxall] * self.beta_mu[:, idxall]).sum(axis=1)
            self.beta_mu[:, k] = (ytX - np.dot(beta_all_k, XtX)) / self.beta_post_tau[:, k] * self.y_tau  # update mu
            u = -0.5 * np.log(self.beta_post_tau[:, k]) + np.log(self.prior_pi.transpose()) + \
                0.5 * self.beta_mu[:, k] ** 2 * self.beta_post_tau[:, k]
            self.gamma[:, k] = softmax(u)  # update gamma

    def get_elbo(self, XX, ytX, XtX):
        beta_all = (self.gamma * self.beta_mu).sum(axis=1)
        ll1 = - 2 * np.dot(beta_all, ytX)
        ll2 = ((self.gamma * self.beta_mu ** 2).sum(axis=1) * XX).sum()
        W = self.gamma * self.beta_mu
        WtRW = np.dot(np.dot(W.transpose(), XtX), W)
        ll3 = WtRW.sum() - np.diag(WtRW).sum()
        ll = -0.5 * self.y_tau * (ll1 + ll2 + ll3)
        betaterm1 = -0.5 * (self.beta_prior_tau * self.gamma * (self.beta_mu ** 2)).sum()
        gammaterm1 = (self.gamma * np.tile(self.prior_pi.reshape(-1, 1), (1, self.k))).sum()
        gammaterm2 = (self.gamma[self.gamma != 0] * np.log(self.gamma[self.gamma != 0])).sum()
        eut = -0.5 * (self.gamma * np.log(self.beta_post_tau)).sum()  # extra unstandardized term
        mkl = betaterm1 + gammaterm1 - gammaterm2 + eut
        elbo = ll + mkl
        return ll, mkl, elbo

    def get_PIP(self):
        return np.max(self.gamma, axis=1)

    def update_pi(self, new_pi):
        self.prior_pi = new_pi

    def get_effect(self, cthres=0.95, ethres=20):
        vidx = np.argsort(-self.gamma, axis=1)
        matidx = np.argsort(-self.gamma, axis=0)
        mat_eff = np.zeros((self.p, self.k)) # effective gamma
        for p in range(self.p):
            mat_eff[p, vidx[p, 0]] = self.gamma[p, vidx[p, 0]]
        csum = mat_eff.sum(axis=0).round(2)
        print("Attainable coverage for effect groups: {}".format(csum))  # statistical evidence
        eff = {}
        eff_gamma = {}
        eff_mu = {}
        for k in range(self.k):
            if csum[k] > cthres:
                if entropy(mat_eff[:, k]) < np.log(ethres):
                    for p in range(self.p):
                        if np.sum(mat_eff[matidx[0:p, k], k]) > cthres * csum[k] or mat_eff[matidx[p, k], k] < 0.01:
                            eff[k] = matidx[0:p, k]
                            eff_gamma[k] = mat_eff[eff[k], k].round(4)
                            eff_mu[k] = self.beta_mu[eff[k], k].round(4)
                            break
        return eff, eff_gamma, eff_mu

    def train(self, XX, ytX, XtX, maxite=100, eps=0.01, verbose=False, loss=0.0):
        for ite in range(maxite):
            self.infer_q_beta(ytX, XtX)
            ll, mkl, elbo = self.get_elbo(XX, ytX, XtX)
            if verbose:
                print('*' * 70)
                print('Iteration-->{} . Likelihood: {:.2f} . KL: {:.2f} . ELBO: {:.2f}'.format(ite, ll, mkl, elbo))
            if abs(elbo - loss) < eps:
                break
            if ite == (maxite - 1):
                print("Algorithm not converged. Please make sure matched summary statistics and LD were provided!")
            loss = elbo


def get_XX_XtX_ytX_z(LD, Z, N):
    XX = np.ones(len(Z)) * N
    XtX = LD * N
    ytX = Z * np.sqrt(N)
    return XX, XtX, ytX


def get_HESS_h2_z(LD, Z, N, ptLD=0.2, ptp=1e-5):
    """calculate local heritabilities"""
    zsquare = Z ** 2
    pvalmin = chi2.sf(max(zsquare), 1)
    idx_retain = []
    idx_exclude = [i for i in range(len(Z))]
    while len(idx_exclude) > 0:  # P + T
        maxid = idx_exclude[np.argmax(zsquare[idx_exclude])]  # find the idx with smallest p-value
        pval = chi2.sf(zsquare[maxid], 1)
        idx_retain.append(maxid)
        idx_exclude = [i for i in idx_exclude if i not in np.where(np.abs(LD[maxid, :]) > ptLD)[0]]
        if pval > ptp:
            break
    Indidx = np.sort(idx_retain)  # obtain independent signals
    P = len(Indidx)
    LD_id = LD[np.ix_(Indidx, Indidx)]
    R_inv = np.linalg.inv(LD_id)
    vec_id = Z[Indidx]
    h2_hess = (np.dot(np.dot(vec_id.transpose(), R_inv), vec_id) - P) / (N - P)
    var_b = np.max(zsquare[Indidx]) / N
    if h2_hess < 0.0001:
        h2_hess = 0.0001
    if h2_hess > 0.2:
        h2_hess = 0.2
        print("Heritability estimates is exceeding 0.2, which is rare for a locus, please check!!")
    return h2_hess, var_b, pvalmin, P


def get_ld_ukb(ldlists, ite, z, args):
    """Get match GWAS summary statistics for each block"""
    ld, start, end = ldlists.iloc[ite, 0:3]
    ldfile = ld.replace('.npz', '')
    df_R, df_ld_snps = load_ld_npz(os.path.join(args.LDdir, ldfile))
    idx = df_R.index.intersection(z.index)
    if len(idx) < 10:
        print("Less than 10 variants matched in the range of {} to {}, skipping".format(start, end))
        return [], [], [], [], [], [], [], [], []
    pos = [int(i.split('.')[1]) for i in idx]
    effidx = [i for i in range(len(idx)) if ((pos[i] >= start) & (pos[i] < end))]
    effnum = len(effidx)
    if effnum <= 10:
        print('Less than 10 effective variants in the range of {} to {}, skipping'.format(start, end))
        return [], [], [], [], [], [], [], [], []
    print('{} variants in the range of {} to {}'.format(effnum, start, end))
    LD = df_R.loc[idx, idx]
    if args.h2:
        h2_hess, var_b, p_chi, K = ldlists.loc[ite, ['h2', 'varb', 'pval', 'K']]
        if np.isnan(K):
            return [], [], [], [], [], [], [], [], []
    else:
        h2_hess, var_b, p_chi, K = get_HESS_h2_z(LD.values, z.loc[idx, 1].values, args.N, ptp=args.gwp)
    if args.K is not None:
        K = args.K
    if p_chi > args.gwp:
        print('The smallest p-value in this region is larger than defined threshold, continue to next block\n')
        return [], [], [], [], [], [], [], [], []
    Z = z.loc[idx, 1].values
    LD = df_R.loc[idx, idx]
    XX, XtX, ytX = get_XX_XtX_ytX_z(LD, Z, args.N)
    return idx, effidx, XX, ytX, XtX, h2_hess, var_b, int(K), p_chi


def ukb(args):
    print("Using genome-wide fine-mapping mode with --ukb")
    if args.anno is None:
        print("Statistical fine-mapping...")
    elif args.aW is not None:
        print("Annotated fine-mapping...")
        anno = pd.read_csv(args.anno, sep='\t', index_col=0)
        W_sig = pd.read_csv(args.aW, sep='\t')
        sigidx = W_sig['sigidx'].tolist()
        W_new = W_sig['W'].values
        new_pi_vec = []
    else:
        print("Please choose either statistical fine-mapping or annotated fine-mapping")
        sys.exit()
    z = pd.read_csv(args.zdir, sep="\t", header=None, index_col=0)  # no header, first column idx(chr.pos.ref.alt)
    print("Summary statistics loaded at {}\n".format(time.strftime("%Y-%m-%d %H:%M")))
    ldlists = pd.read_csv(args.ukb, sep='\t')  # header with 3 columns
    pip = []
    pip_name = []
    cs = []
    cs_pip = []
    cs_eff = []
    z_vec = []
    for ite in range(len(ldlists)):
        idx, effidx, XX, ytX, XtX, h2_hess, var_b, K, p_chi = get_ld_ukb(ldlists, ite, z, args)
        if len(effidx) == 0:
            continue
        model = SparsePro(len(idx), K, XX, h2_hess, var_b)
        if args.anno is not None and args.aW is not None:
            ianno = anno.loc[idx]
            ANN = ianno.values[:, sigidx]
            new_pi = softmax(np.dot(ANN, W_new))
            model.update_pi(new_pi)
            new_pi_vec.extend([new_pi[i] for i in effidx])
        model.train(XX, ytX, XtX, verbose=args.verbose)
        mcs, eff_gamma, eff_mu = model.get_effect(args.cthres, args.ethres)
        pip_vec = model.get_PIP().round(4)
        pip.extend([pip_vec[i] for i in effidx])
        z_vec.extend([z.loc[idx, 1].values[i] for i in effidx])
        pip_name.extend([idx[i] for i in effidx])
        if not args.h2:
            ldlists.at[ite, 'h2'] = '{:.2e}'.format(h2_hess)
            ldlists.at[ite, 'pval'] = '{:.2e}'.format(p_chi)
            ldlists.at[ite, 'varb'] = '{:.2e}'.format(var_b)
            ldlists.at[ite, 'K'] = '{:.2e}'.format(K)
        if len(mcs) == 0:
            print("No effect detected")
            print()
            continue
        for e in mcs:
            if mcs[e][0] in effidx:
                mcs_idx = [idx[j] for j in mcs[e]]
                print('The {}-th effect contains effective variants:'.format(e))
                print('causal variants: {}'.format(mcs_idx))
                print('posterior inclusion probabilities: {}'.format(eff_gamma[e]))
                print('posterior causal effect size: {}'.format(eff_mu[e]))
                print()
                cs.append(mcs_idx)
                cs_pip.append(eff_gamma[e])
                cs_eff.append(eff_mu[e])
    if args.anno is not None and args.aW is not None:
        allPIP = pd.DataFrame({"idx": pip_name,
                               "z": z_vec,
                               "pip": pip,
                               "prior": new_pi_vec})
    else:
        allPIP = pd.DataFrame({"idx": pip_name,
                               "z": z_vec,
                               "pip": pip})
    allPIP.to_csv(os.path.join(args.save, "{}.pip".format(args.prefix)), sep='\t', header=False, index=False)
    allcs = pd.DataFrame({"cs": ['/'.join(i) for i in cs],
                          "pip": ['/'.join([str(j) for j in i]) for i in cs_pip],
                          "beta": ['/'.join([str(j) for j in i]) for i in cs_eff]})
    allcs.to_csv(os.path.join(args.save, "{}.cs".format(args.prefix)), sep='\t', header=True, index=False)
    ldlists.to_csv(os.path.join(args.save, "{}.h2".format(args.prefix)), sep='\t', header=True, index=False)


parser = argparse.ArgumentParser(description='SparsePro Commands:')
parser.add_argument('--ukb', type=str, default=None, help='genome-wide finemapping mode: path to LD lists')
parser.add_argument('--zdir', type=str, default=None, help='path to zscores files', required=True)
parser.add_argument('--LDdir', type=str, default=None, help='path to ld files', required=True)
parser.add_argument('--N', type=int, default=None, help='GWAS sample size', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', action="store_true", help='options for displaying more information')
parser.add_argument('--anno', type=str, default=None, help='path to annotation file')
parser.add_argument('--K', type=int, default=None, help='largest number of effect')
parser.add_argument('--gwp', type=float, default=1e-5, help='p-value threshold for filtering blocks to be fine-mapped')
parser.add_argument('--cthres', type=float, default=0.95, help='coverage level for effect group')
parser.add_argument('--ethres', type=float, default=20.0, help='entropy level for effect group')
parser.add_argument('--aW', type=str, default=None, help='enrichment file')
parser.add_argument('--h2', action="store_true", help='use previous h2 file as zld file')

args = parser.parse_args()
title()
if not os.path.exists(args.save):
    os.makedirs(args.save)

ukb(args)
