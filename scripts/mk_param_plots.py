#! /usr/bin/env python

import numpy as np
from pylab import *

VLSS = np.load('data/params_sampled_VLSS.npz')
SixC = np.load('data/params_sampled_6C_hiRM.npz')

def Normal(x,mu,sig):
    G = -0.5*(x-mu)**2 / sig**2
    return np.exp(G)/np.sqrt(2.*np.pi*sig**2)

def LogNormal(x,mu,sig):
    Vx = np.log(1.+(sig/mu)**2)
    P = -0.5*np.log(x/mu)**2 / Vx
    return np.exp(P) / (x*np.sqrt(2.*np.pi*Vx))

def SourceCounts(S, dNdS=False, log=True):
    Nsrc = len(S)
    if log:
        Plogs,logs = np.histogram(np.log(S), 30.)
        s = np.exp(0.5*(logs[1:]+logs[:-1]))
        Ps = Plogs/s
    else:
        Ps,s = np.histogram(S,bins=len(S)/100.)
        s = 0.5*(s[1:] + s[:-1])
    if dNdS:
        return s, Ps/0.76
    cdf = []
    for i in range(len(s)):
        cdf.append(np.sum(Ps[:i]))
    return s, (cdf[-1] - cdf)/0.76

def dNdS_6C(S, log=True):
    def lo(x): return (4000.*(0.88)**-0.76) * x**-1.75
    def hi(x): return 4000. * x**-2.81
    dNdS = np.where(S >= 0.88, hi(S), lo(S))
    if log:
        dS = np.ones_like(S)*(np.log(S[1])-np.log(S[0]))
    else:
        dS = np.ones_like(S)*(S[1]-S[0])
    NgtS = []
    for i in range(len(dNdS)):
        NgtS.append(np.sum(dNdS[:i]*dS[:i]))
    return (NgtS[-1] - NgtS)

def dNdS_VLSS(S, log=True):
    dNdS = 4865.* (S*(150./94.)**0.76)**-2.3
    if log:
        dS = np.ones_like(S)*(np.log(S[1])-np.log(S[0]))
    else:
        dS = np.ones_like(S)*(S[1]-S[0])
    NgtS = []
    for i in range(len(dNdS)):
        NgtS.append(np.sum(dNdS[:i]*dS[:i]))
    return (NgtS[-1] - NgtS)

ParamFig = figure(figsize=(8,8))
#Plot the Rotation Measure Distribution
RMax = ParamFig.add_subplot(212)

Prm,rm = np.histogram(VLSS['RM'], bins=len(VLSS['RM'])/30., density=True)
rm = 0.5*(rm[1:]+rm[:-1])
RMax.step(rm,Prm,'k',where='mid')
RMax.set_xlabel(r'$\Phi\ [\rm{m}^{-2}]$', size=15)
RMax.set_xlim([-200,200])
RMax.grid(True)

#Polarized Fraction
Fax  = ParamFig.add_subplot(222)
F = SixC['Pflx']/SixC['flx']
Pf,f = np.histogram(np.log10(F), bins=30, density=True)
f = 0.5*(f[1:] + f[:-1])

Fax.step(f, Pf, 'k', where='mid')
Fax.plot(f, Normal(f,np.log10(0.0201),np.log10(1+(0.0304/0.0201)**2)), '0.5')

Fax.set_xlabel(r'$\log_{10}(\Pi)$', size=15)
Fax.grid(True)
Fax.yaxis.tick_right()
Fax.set_xticks([-4.,-3.,-2.,-1.])

#Source Counts
Sax =  ParamFig.add_subplot(221)

S,NgtS = SourceCounts(SixC['flx'])
Sax.step(S, S**1.5*NgtS,'b', where='mid')
Sax.plot(S, S**1.5*dNdS_6C(S), 'c')

S,NgtS = SourceCounts(VLSS['flx'])
Sax.step(S, S**1.5*NgtS,'k', where='mid')
Sax.plot(S, S**1.5*dNdS_VLSS(S), '0.5')

Sax.set_yscale('log')
Sax.set_xscale('log')
Sax.set_ylim([1e1,1e4])
Sax.set_xlim([1e-1,1e2])
Sax.set_ylabel(r'$S^{3/2} N(>S)\ [\rm{Jy}^{3/2}\ \rm{sr}^{-1}]$', size=15)
Sax.set_xlabel(r'$S\ [\rm{Jy}]$', size=15)
Sax.grid(True)

Pfig = figure()
ax = Pfig.add_subplot(111)
P,NgtP = SourceCounts(SixC['Pflx'],dNdS=True)
ax.step(P, P**2.5 * NgtP,'b', where='mid')
Pmean = np.sum(P*NgtP)/np.sum(NgtP)
Pvar = np.sum(P**2*NgtP)/np.sum(NgtP)
Pvar = np.sqrt(Pvar - Pmean**2)
print "6C counts: <P> = %.2f +/- %.2f mJy"%(1000.*Pmean, 1000.*Pvar)

P,NgtP = SourceCounts(VLSS['Pflx'],dNdS=True)
ax.step(P, P**2.5 * NgtP,'k', where='mid')
Pmean = np.sum(P*NgtP)/np.sum(NgtP)
Pvar = np.sum(P**2*NgtP)/np.sum(NgtP)
Pvar = np.sqrt(Pvar - Pmean**2)
print "NVSS counts: <P> = %.2f +/- %.2f mJy"%(1000.*Pmean, 1000.*Pvar)

ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.set_ylabel(r'$P^{5/2}\ dN/dP\ [\rm{Jy}^{3/2}\ \rm{sr}^{-1}]$', size=15)
ax.set_xlabel(r'$P\ [\rm{Jy}]$', size=15)

ParamFig.savefig('Pspec/figures/SimParams.eps', fmt='eps')
Pfig.savefig('Pspec/figures/SimPolFlux.eps', fmt='eps')
show()
