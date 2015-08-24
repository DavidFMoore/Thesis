#! /usr/bin/env python

import numpy as np
from pylab import *
from capo import pspec
from scipy.stats import chi2
from scipy.special import erf

zs = ['7.05','9.77','11.13']
sims = ['6C', 'corr']#,'VLSS','hiRM']
colors = ['b','k','c','0.5']

fig = figure()

def Nsig(N):
    return 0.5*(1 + erf(N/np.sqrt(2.)))

def get_stats(P):
    pdfPk = []
    for i in range(P.shape[1]):
        N = float(P.shape[0])
        var = np.sum(P[:,i])/(N*(4.-2.))
        pdfPk.append(chi2(4.,scale=var))
    return pdfPk

def k3pk(k,p):
    return np.abs(k)**3 * p / (2.*np.pi**2)

def fmt_axis(ax):
    ax.set_yscale('log', nonposy='clip')
    #ax.set_xscale('log')
    ax.set_xlim([0.01,0.5])
    ax.set_ylim([10**1,10**9])
    ax.grid(True)

def rm_xlabel(ax):
    for xl in ax.get_xticklabels():
        xl.set_visible(False)

def fmt_row(iax,qax,index):
    iax.set_ylabel(r'$\Delta^2(k)\ [\rm{mK}^2]$', size=15)
    simlabel = "ABCDEF"[index]
    iax.text(0.05, 100., simlabel, size=15)
    qax.text(0.05, 100., simlabel, size=15)
    if index == 0:
        iax.set_title(r'$\mathcal{V}_{I} = \mathcal{V}_{xx} + \mathcal{V}_{yy}$', size=15)
        qax.set_title(r'$\mathcal{V}_{Q} = \mathcal{V}_{xx} - \mathcal{V}_{yy}$', size=15)
    if index < len(sims)-1:
        rm_xlabel(qax)
        rm_xlabel(iax)
    else:
        iax.set_xticks([0,0.1,0.2,0.3,0.4])
        iax.set_xlabel(r'$k\ [h\rm{Mpc}^{-1}]$', size=15)
        qax.set_xlabel(r'$k\ [h\rm{Mpc}^{-1}]$', size=15)
    for yl in qax.get_yticklabels():
        yl.set_visible(False)

has_model = [False]*len(sims)
for zi,z in enumerate(zs):
    for simj,sim in enumerate(sims):
        try:
            P = np.load('data/%s_%s.npz'%(sim,z))
        except(IOError):
            continue
        PI = 2.35*P['I'] #2.35 to correct for power-squared beam.
        PQ = 2.35*P['Q']
        if sim == 'VLSS':
            flx_scale = (150./74.)**(2.*0.79)
            PI *= flx_scale
            PQ *= flx_scale
        k = P['k']
        Iax = fig.add_subplot(1,2,1)
        Qax = fig.add_subplot(1,2,2)
        fmt_row(Iax, Qax, simj)
        for pk,ax in zip([PI,PQ],[Iax,Qax]):
            pk = get_stats(pk)
            p_med = np.array([_pk.isf(Nsig( 0.)) for _pk in pk])[:len(pk)/2]
            p_hi  = np.array([_pk.isf(Nsig(-1.)) for _pk in pk])[:len(pk)/2]
            p_lo  = np.array([_pk.isf(Nsig( 1.)) for _pk in pk])[:len(pk)/2]
            _k = k[0][:len(pk)/2]
            #ax.plot(_k,p_lo,color=colors[zi])
            if sim == 'corr':
                ax.errorbar(_k,p_med,yerr=[p_med-p_lo, p_hi-p_med],
                        fmt='--',
                        color=colors[zi],
                        ecolor=colors[zi])
            else:
                ax.errorbar(_k,p_med,yerr=[p_med-p_lo, p_hi-p_med],
                        color=colors[zi],
                        ecolor=colors[zi])
            fmt_axis(ax)
        if not has_model[simj]:
            Iax.plot(_k,pspec.k3pk_21cm(_k)*0.7**3,'0.5')
            has_model[simj] = True
fig.subplots_adjust(hspace=0, wspace=0)

fig.savefig('Pspec/figures/SimulationResults.eps',fmt='eps')
show()
