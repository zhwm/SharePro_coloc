import pandas as pd
import argparse
import os
import numpy as np
from scipy.special import softmax, expit
from scipy.stats import chi2
from scipy.stats import entropy

np.set_printoptions(precision=4, linewidth=200)


def title():
    print('**********************************************************************')
    print('* SharePro for accurate and efficient colocalization                 *')
    print('* Version 1.0.0                                                      *')
    print('* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *')
    print('**********************************************************************')


def get_HESS_h2_z(LD, Z, N, ptLD=0.2, ptp=1e-5):
    """calculate local heritabilities"""
    zsquare = Z ** 2
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
    if h2_hess > 0.9:
        h2_hess = 0.9
    return h2_hess, var_b


class SparseReg(object):
    def __init__(self, P, K, XX, h2, var_b, sigma):
        """initialize and set parameters"""
        self.p = P
        self.k = K
        self.beta_mu = np.zeros([self.p, self.k])
        self.beta_prior_tau = np.tile((1.0 / var_b * np.array([k + 1 for k in range(self.k)])), (self.p, 1))
        self.y_tau = 1.0 / (1 - h2)
        self.sigma = sigma
        self.beta_post_tau = np.tile(XX.reshape(-1, 1), (1, self.k)) * self.y_tau + self.beta_prior_tau
        self.delta = np.zeros((self.p, self.k))
        self.u1 = 0.5 * np.log(self.beta_post_tau / self.beta_prior_tau) + np.log(sigma / (1 - sigma))

    def infer_q_beta(self, ytX, XtX, GAM, k):
        """perform variational updates for the k-th effect with beta and c"""
        idxall = [x for x in range(self.k)]
        idxall.remove(k)
        beta_all_k = (GAM[:, idxall] * self.beta_mu[:, idxall] * self.delta[:, idxall]).sum(axis=1)
        self.beta_mu[:, k] = (ytX - np.dot(beta_all_k, XtX)) / self.beta_post_tau[:, k] * self.y_tau  # Update beta
        u = self.u1[:, k] + 0.5 * self.beta_mu[:, k] ** 2 * self.beta_post_tau[:, k]
        self.delta[:, k] = expit(u)  # Update c
        return u

    def get_elbo(self, XX, ytX, XtX, GAM):
        """get elbo with beta and c"""
        beta_all = (GAM * self.delta * self.beta_mu).sum(axis=1)
        ll1 = self.y_tau * np.dot(beta_all, ytX)
        ll2 = - 0.5 * self.y_tau * (
            (((GAM * self.delta * (self.beta_mu ** 2 + 1 / self.beta_post_tau)).sum(axis=1) * XX).sum()))
        W = GAM * self.delta * self.beta_mu
        WtRW = np.dot(np.dot(W.transpose(), XtX), W)
        ll3 = - 0.5 * self.y_tau * (WtRW.sum() - np.diag(WtRW).sum())
        ll = ll1 + ll2 + ll3
        mklbeta = - 0.5 * (GAM * self.delta * (
                self.beta_prior_tau * (self.beta_mu ** 2 + 1 / self.beta_post_tau) +
                np.log(self.beta_post_tau / self.beta_prior_tau) - 1)).sum()
        nozerodelta = self.delta[self.delta != 0]
        nozeroGAM = GAM[self.delta != 0]
        noonedelta = self.delta[self.delta != 1]
        nooneGAM = GAM[self.delta != 1]
        mkldelta = (nozeroGAM * nozerodelta * np.log(self.sigma / nozerodelta)).sum() - (
                nooneGAM * (1 - noonedelta) * np.log((1 - self.sigma) / (1 - noonedelta))).sum()
        return ll, mklbeta, mkldelta


class SharePro(object):
    def __init__(self, P, K, XX, h2, varb, sigma):
        """initialize and set parameters"""
        self.num = XX.shape[1]
        self.SR = [SparseReg(P, K, XX[:, i], h2[i], varb[i], sigma) for i in range(self.num)]
        self.p = P
        self.k = K
        self.gamma = np.zeros((self.p, self.k))
        self.sigma = sigma
        self.prior_pi = np.ones((self.p,)) * (1 / self.p)
        self.usum1 = self.num * np.log(1 - self.sigma) + np.log(self.prior_pi.transpose())

    def infer_q_s(self, ytX, XtX):
        """perform variational updates for s"""
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            u12 = np.array(
                [self.SR[i].infer_q_beta(ytX[:, i], XtX[i], self.gamma, k) for i in range(len(self.SR))])
            usum = np.log(1 + np.exp(-u12))
            uall = u12.sum(axis=0) + usum.sum(axis=0) + self.usum1
            self.gamma[:, k] = softmax(uall)

    def get_elbo(self, XX, ytX, XtX):
        """get elbo with s"""
        llsum, mklbetasum, mkldeltasum = zip(
            *[self.SR[i].get_elbo(XX[:, i], ytX[:, i], XtX[i], self.gamma) for i in range(len(self.SR))])
        gammaterm1 = (self.gamma * np.tile(self.prior_pi.reshape(-1, 1), (1, self.k))).sum()
        gammaterm2 = (self.gamma[self.gamma != 0] * np.log(self.gamma[self.gamma != 0])).sum()
        mklgamma = gammaterm1 - gammaterm2
        ll = sum(llsum)
        mklbeta = sum(mklbetasum)
        mkldelta = sum(mkldeltasum)
        elbo = ll + mklbeta + mkldelta + mklgamma
        return ll, mklbeta, mkldelta, mklgamma, elbo

    def train(self, XX, ytX, XtX, maxite=100, eps=0.01, verbose=True, loss=0.0):
        for ite in range(maxite):
            self.infer_q_s(ytX, XtX)
            ll, mklbeta, mkldelta, mklgamma, elbo = self.get_elbo(XX, ytX, XtX)
            if verbose:
                print('*' * 70)
                print('Iteration-->{} . Likelihood: {:.1f} . KL_b: {:.1f} . KL_c: {:.1f} . KL_s: {:.1f} . ELBO: {:.1f}'.
                      format(ite, ll, mklbeta, mkldelta, mklgamma, elbo))
            if abs(elbo - loss) < eps:
                break
            if ite == (maxite - 1):
                print("Algorithm not converged. Please make sure matched summary statistics and LD were provided!")
            loss = elbo

    def multiply_specific(self, i):
        """calculate c1(1-c2)*s"""
        allidx = [x for x in range(self.num)]
        allidx.remove(i)
        matspecific = self.SR[i].delta * (np.prod(np.array([(1 - self.SR[x].delta) for x in allidx]), axis=0))
        return matspecific

    def multiply_delta(self):
        """calculate c1*c2*s"""
        return np.prod(np.array([i.delta for i in self.SR]), axis=0)

    def get_effect(self, cthres=0.95, ethres=50):
        vidx = np.argsort(-self.gamma, axis=1)
        matidx = np.argsort(-self.gamma, axis=0)
        mat_eff = np.zeros((self.p, self.k))  # effective gamma
        for p in range(self.p):
            mat_eff[p, vidx[p, 0]] = self.gamma[p, vidx[p, 0]]
        matdelta = self.multiply_delta()
        csum = mat_eff.sum(axis=0).round(2)
        print("Attainable coverage for effect groups: {}".format(csum))  # statistical evidence
        eff = {}
        eff_gamma = {}
        eff_share = {}
        for k in range(self.k):
            if csum[k] > cthres:
                if entropy(mat_eff[:, k]) < np.log(ethres):
                    for p in range(self.p):
                        if np.sum(mat_eff[matidx[0:p, k], k]) > cthres * csum[k] or mat_eff[matidx[p, k], k] < 0.01:
                            eff[k] = matidx[0:p, k].tolist()
                            eff_gamma[k] = mat_eff[eff[k], k].round(4).tolist()
                            effgamma_n = eff_gamma[k] / csum[k]
                            eff_share[k] = sum(np.multiply(matdelta[eff[k], k], effgamma_n)).round(4)
                            break
        return eff, eff_gamma, eff_share


def zld(args):
    ldlists = pd.read_csv(args.zld, sep='\s+')  # header with 2 columns: z\tld
    print("LD list with {} LD blocks loaded\n".format(len(ldlists)))
    for ite in range(len(ldlists)):
        zfile, ldfile = ldlists.iloc[ite, 0:2]
        print("processing {}".format(zfile))
        zs = zfile.split(',')
        lds = ldfile.split(',')
        nums = len(zs)
        z = pd.concat([pd.read_csv(os.path.join(args.zdir, zs[i]), sep='\t', header=None, index_col=0)
                       for i in range(nums)], axis=1, join='inner')
        ldmat = [pd.read_csv(os.path.join(args.zdir, lds[i]), sep='\s+', header=None).values for i in range(nums)]
        assert len(z) == len(ldmat[0])
        assert all(len(i) == len(ldmat[0]) for i in ldmat)
        Z = z.values
        XX = np.ones(Z.shape) * args.N
        ytX = Z * np.sqrt(args.N)
        XtX = [i*j for i, j in zip(ldmat, args.N)]
        hess, varb = zip(*[get_HESS_h2_z(ldmat[i], Z[:, i], args.N[i], ptLD=args.ptLD, ptp=args.ptp)
                           for i in range(nums)])
        if args.hess is not None:
            hess = args.hess
        if args.varb is not None:
            varb = args.varb
        model = SharePro(Z.shape[0], args.K, XX, hess, varb, sigma=args.sigma)
        model.train(XX, ytX, XtX, verbose=args.verbose)
        eff, eff_gamma, eff_share = model.get_effect(cthres=args.cthres, ethres=args.ethres)
        for e in eff:
            mcs_idx = [z.index[j] for j in eff[e]]
            print('The {}-th effect group contains effective variants:'.format(e))
            print('causal variants: {}'.format(mcs_idx))
            print('variant probabilities for this effect group: {}'.format(eff_gamma[e]))
            print('shared probability for this effect group: {}'.format(eff_share[e]))
            print()
        ldlists.at[ite, 'h2'] = ','.join(['{:.2e}'.format(i) for i in hess])
        ldlists.at[ite, 'varb'] = ','.join(['{:.2e}'.format(i) for i in varb])
        allcs = pd.DataFrame.from_dict([eff, eff_share, eff_gamma]).transpose()
        allcs['cs'] = ['/'.join([z.index[j] for j in i]) for i in allcs[0]]
        allcs['share'] = allcs[1]
        allcs['variantProb'] = ['/'.join([str(j) for j in i]) for i in allcs[2]]
        allcs[['cs', 'share', 'variantProb']].to_csv(os.path.join(args.save, "{}.cs".format(zfile.replace(',', '_'))),
                     sep='\t', header=True, index=False)
    ldlists.to_csv(os.path.join(args.save, "{}.h2".format(args.prefix)), sep='\t', header=True, index=False)


parser = argparse.ArgumentParser(description='SharePro Commands:')
parser.add_argument('--zld', type=str, default=None,
                    help='index file contains path to matched zscore and ld lists', required=True)
parser.add_argument('--zdir', type=str, default=None, help='path to zscores files', required=True)
parser.add_argument('--N', type=int, default=None, nargs='+', help='sample sizes', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save results', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', action="store_true", help='options for displaying more information')
parser.add_argument('--K', type=int, default=None, help='largest number of causal signals', required=True)
parser.add_argument('--sigma', type=float, default=1e-5, help='prior colocalization probabilities')
parser.add_argument('--hess', type=float, default=None, nargs='+', help='heritability estimates, HESS estimates used as default')
parser.add_argument('--varb', type=float, default=None, nargs='+',
                    help='effect size variance estimated, estimates from HESS used as default')
parser.add_argument('--ptLD', type=float, default=0.2, help='P+T LD cutoff')
parser.add_argument('--ptp', type=float, default=1e-5, help='P+T p value cutoff')
parser.add_argument('--cthres', type=float, default=0.95, help='coverage level for credible sets')
parser.add_argument('--ethres', type=float, default=50.0, help='entropy level for credible sets')


args = parser.parse_args()
title()
if not os.path.exists(args.save):
    os.makedirs(args.save)

zld(args)