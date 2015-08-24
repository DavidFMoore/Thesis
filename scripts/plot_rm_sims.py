#! /usr/bin/env python

import numpy as np
from pylab import *
from capo import pspec
from glob import glob
import re
from scipy.stats import chi2
from scipy.special import erf

fig = figure()
ax = fig.add_subplot(111)

def Nsig(N):
    return 0.5*(1 + erf(N/np.sqrt(2.)))

def get_stats(P):
    N = float(P.shape[0])
    var = np.sum(P)/(N*(4.-2.))
    return chi2(4.,scale=var)

sims = ['6C', 'RM1k', 'RM2k', 'RM5k', 'RM10k']
zs = {}
Pk = {}
for sim in sims:
    zs[sim] = []
    Pk[sim] = []
    for f in glob('data/%s*.npz'%sim):
        z = float(re.compile('\d+\.\d+').findall(f)[0])
        _pk = np.load(f)
        _k = _pk['k'][0]
        _pk = _pk['I']
        i = np.argmin(np.abs(_k - 0.25))
        zs[sim].append(z)
        Pk[sim].append(get_stats(_pk[:,i]))
    i = np.argsort(zs[sim])
    zs[sim] = np.array(zs[sim]).take(i)
    Pk[sim] = np.array(Pk[sim]).take(i)

colors = ['c','b','g','0.5','k']
for i,sim in enumerate(sims):
    print colors[i], sim
    p_med = np.array([_pk.isf(Nsig( 0.)) for _pk in Pk[sim]])
    p_hi  = np.array([_pk.isf(Nsig(-1.)) for _pk in Pk[sim]])
    p_lo  = np.array([_pk.isf(Nsig( 1.)) for _pk in Pk[sim]])

    ax.errorbar(zs[sim], p_med, yerr=[p_med-p_lo, p_hi-p_med],
            color=colors[i],
            ecolor=colors[i])

ax.set_yscale('log', nonposy='clip')
ax.grid(True)
ax.set_ylim([10**0,10**4.5])
ax.set_ylabel(r'$\Delta^2(k)\ [\rm{mK}^2]$', size=15)
ax.set_xlabel(r'$z$', size=15)

fig.savefig('Pspec/figures/DeltakGoodK_rm.eps', fmt='eps')
show()
