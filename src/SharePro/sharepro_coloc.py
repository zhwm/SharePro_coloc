import argparse
import logging
import pandas as pd
import numpy as np
from scipy.special import softmax, expit

np.set_printoptions(precision=4, linewidth=200)

def title():
    print('**********************************************************************')
    print('* SharePro for accurate and efficient colocalization                 *')
    print('* Version 5.0.0                                                      *')
    print('* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *')
    print('**********************************************************************')


class SparseReg(object):
    def __init__(self, P, K, N, sigma, varb):
        """initialize and set parameters"""
        self.p = P
        self.k = K
        self.n = N
        self.sigma = sigma
        self.y_tau = 1.0 / (1 - varb)
        self.priorb = 1.0 / varb
        self.postn = self.n * self.y_tau + self.priorb
        self.lamb = self.priorb / self.postn
        self.beta_mu = np.zeros([self.p, self.k])
        self.delta = np.zeros((self.p, self.k))
        self.u1 = np.log(self.sigma / (1 - self.sigma)) + 0.5 * np.log(self.lamb)

    def infer_q_beta(self, b_hat_s, R_s, GAM, k):
        """perform variational updates for the k-th effect with beta and c"""
        idxall = list(range(self.k))
        idxall.remove(k)
        beta_all_k = (GAM[:, idxall] * self.beta_mu[:, idxall] * self.delta[:, idxall]).sum(axis=1)
        res_beta = b_hat_s - np.dot(R_s, beta_all_k)
        self.beta_mu[:, k] = (1 - self.lamb) * res_beta
        u = self.u1 + 0.5 * self.beta_mu[:, k] ** 2 * self.postn
        self.delta[:, k] = expit(u)
        return u

    def get_elbo(self, b_hat_s, R_s, GAM):
        """get elbo with beta and c"""
        beta_all = (GAM * self.delta * self.beta_mu).sum(axis=1)
        ll1 = self.y_tau * np.dot(beta_all, self.n * b_hat_s)
        ll2 = - 0.5 * self.y_tau * (
            (((GAM * self.delta * (self.beta_mu ** 2 + 1 / self.postn)).sum(axis=1) * self.n).sum()))
        W = GAM * self.delta * self.beta_mu
        WtRW = np.dot(np.dot(W.transpose(), self.n * R_s), W)
        ll3 = - 0.5 * self.y_tau * (WtRW.sum() - np.diag(WtRW).sum())
        ll = ll1 + ll2 + ll3
        mklbeta = - 0.5 * (GAM * self.delta * (
                self.priorb * (self.beta_mu ** 2 + 1 / self.postn) +
                np.log(self.postn / self.priorb) - 1)).sum()
        nozerodelta = self.delta[self.delta != 0]
        nozeroGAM = GAM[self.delta != 0]
        noonedelta = self.delta[self.delta != 1]
        nooneGAM = GAM[self.delta != 1]
        mkldelta = (nozeroGAM * nozerodelta * np.log(self.sigma / nozerodelta)).sum() - (
                nooneGAM * (1 - noonedelta) * np.log((1 - self.sigma) / (1 - noonedelta))).sum()
        return ll, mklbeta, mkldelta


class SharePro(object):
    def __init__(self, P, K, num, sigma, Nlst, varblst):
        """initialize and set parameters"""
        self.p = P
        self.k = K
        self.num = num
        self.sigma = sigma
        self.prior_pi = np.ones((self.p,)) * (1 / self.p)
        self.gamma = np.zeros((self.p, self.k))
        self.usum1 = self.num * np.log(1 - self.sigma**(1/num)) + np.log(self.prior_pi)
        self.SR = [SparseReg(P, K, Nlst[i], self.sigma**(1/num), varblst[i]) for i in range(self.num)]

    def infer_q_s(self, bhatlst, ld):
        """perform variational updates for s"""
        for k in range(self.k):
            idxall = list(range(self.k))
            idxall.remove(k)
            unum = np.array([j.infer_q_beta(bhatlst[i], ld, self.gamma, k) for (i,j) in enumerate(self.SR)])
            unumexp = np.log(1 + np.exp(-unum))
            uall = unum.sum(axis=0) + unumexp.sum(axis=0) + self.usum1
            self.gamma[:, k] = softmax(uall)

    def get_effect(self, ld, cthres, pthres):
        vidx = np.argsort(-self.gamma, axis=1)
        matidx = np.argsort(-self.gamma, axis=0)
        mat_eff = np.zeros((self.p, self.k))
        for p in range(self.p):
            mat_eff[p, vidx[p, 0]] = self.gamma[p, vidx[p, 0]]
        matdelta = np.prod(np.array([i.delta for i in self.SR]), axis=0)
        csum = mat_eff.sum(axis=0).round(2)
        logging.info(f"Attainable coverage for effect groups: {csum}")
        eff = {}
        eff_gamma = {}
        eff_share = {}
        for k in range(self.k):
            if csum[k] > cthres:
                p = 0
                while np.sum(mat_eff[matidx[0:p, k], k]) < cthres * csum[k]:
                    p = p + 1
                cidx = matidx[0:p, k].tolist()
                purity = abs(ld[np.ix_(cidx, cidx)]).min()
                if purity > pthres:
                    eff[k] = cidx
                    eff_gamma[k] = mat_eff[eff[k], k].round(4).tolist()
                    effgamma_n = [i/sum(eff_gamma[k]) for i in eff_gamma[k]]
                    eff_share[k] = np.dot(matdelta[eff[k], k], effgamma_n).round(4)
        return eff, eff_gamma, eff_share

    def get_elbo(self, bhatlst, ld):
        """get elbo with s"""
        llsum, mklbetasum, mkldeltasum = zip(*[j.get_elbo(bhatlst[i], ld, self.gamma) for (i,j) in enumerate(self.SR)])
        gammaterm1 = (self.gamma * np.tile(self.prior_pi.reshape(-1, 1), (1, self.k))).sum()
        gammaterm2 = (self.gamma[self.gamma != 0] * np.log(self.gamma[self.gamma != 0])).sum()
        mklgamma = gammaterm1 - gammaterm2
        ll = sum(llsum)
        mklbeta = sum(mklbetasum)
        mkldelta = sum(mkldeltasum)
        elbo = ll + mklbeta + mkldelta + mklgamma
        return ll, mklbeta, mkldelta, mklgamma, elbo

    def train(self, bhatlst, ld, maxite, eps, ubound):
        loss = 100000.0
        for ite in range(maxite):
            self.infer_q_s(bhatlst, ld)
            ll, mklbeta, mkldelta, mklgamma, elbo = self.get_elbo(bhatlst, ld)
            logging.info('*' * 70)
            logging.info('Iteration-->{} . Likelihood: {:.1f} . KL_b: {:.1f} . KL_c: {:.1f} . KL_s: {:.1f} . ELBO: {:.1f}'.format(ite, ll, mklbeta, mkldelta, mklgamma, elbo))
            if abs(elbo - loss) < eps:
                converged = True
                break
            if ite == (maxite - 1) or elbo > ubound:
                print("Detected mismatch between summary statistics and LD matrix. Set K to 1...")
                converged = False
                break
            loss = elbo
        return converged

def get_bhat(beta, se, n):
    z = beta / se
    bhat = z / np.sqrt(z**2 + n)
    return bhat

def parse_args():
    parser = argparse.ArgumentParser(description='SharePro Commands:')
    parser.add_argument('--z', type=str, default=None, nargs='+', help='file to matched summary statistics', required=True)
    parser.add_argument('--ld', type=str, default=None, help='path to matched ld matrix', required=True)
    parser.add_argument('--save', type=str, default=None, help='path to save results', required=True)
    parser.add_argument('--sigma', type=float, default=1e-5, help='prior colocalization prob')
    parser.add_argument('--K', type=int, default=10, help='largest number of causal signals')
    parser.add_argument('--maxite', type=int, default=100, help='max number of iterations')
    parser.add_argument('--eps', type=float, default=0.1, help='convergence criterion')
    parser.add_argument('--ubound', type=int, default=1000000, help='upper bound for inconvergence')
    parser.add_argument('--cthres', type=float, default=0.95, help='attainable coverage threshold for effect groups')
    parser.add_argument('--pthres', type=float, default=0.5, help='purity threshold for effect groups')
    args = parser.parse_args()
    return args

def print_args(args):
    for arg in vars(args):
        logging.info(f"{arg}: {getattr(args, arg)}")

def main(args):
    zfile = [pd.read_csv(i, sep='\t') for i in args.z]
    ld = pd.read_csv(args.ld, sep='\s+', header=None).values
    bhatlst = [get_bhat(i['BETA'].values, i['SE'].values, i['N'].values) for i in zfile]
    Nlst = [max(i['N']) for i in zfile]
    varb = [(i**2).max() for i in bhatlst]
    model = SharePro(zfile[0].shape[0], args.K, len(zfile), args.sigma, Nlst, varb)
    mc = model.train(bhatlst, ld, args.maxite, args.eps, args.ubound)
    if mc:
        eff, eff_gamma, eff_share = model.get_effect(ld, args.cthres, args.pthres)
    else:
        model = SharePro(zfile[0].shape[0], 1, len(zfile), args.sigma, Nlst, varb)
        eff, eff_gamma, eff_share = model.get_effect(ld, args.cthres, args.pthres)
    for e in eff:
        mcs_idx = [zfile[0]['SNP'][j] for j in eff[e]]
        logging.info(f'The {e}-th effect group contains effective variants:')
        logging.info(f'causal variants: {mcs_idx}')
        logging.info(f'variant probabilities for this effect group: {eff_gamma[e]}')
        logging.info(f'shared probability for this effect group: {eff_share[e]}\n')
    allcs = pd.DataFrame()
    allcs['cs'] = ['/'.join(zfile[0]['SNP'][j]) for j in eff.values()]
    allcs['share'] = ['{:.4}'.format(j) for j in eff_share.values()]
    allcs['variantProb'] = ['/'.join([str(j) for j in i]) for i in eff_gamma.values()]
    allcs.to_csv('{}.sharepro.txt'.format(args.save), sep='\t', header=True, index=False)

if __name__ == '__main__':
    args = parse_args()
    logging.basicConfig(filename='{}.sharepro.log'.format(args.save), level=logging.INFO, filemode='w')
    print_args(args)
    main(args)