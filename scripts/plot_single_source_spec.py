#! /usr/bin/env python

import numpy as np
from capo import pspec,pfb
from aipy import dsp
from pylab import *

c = 0.3 #m/ns

window = 'blackman-harris'

def nu_lam_swap(x):return c/x

def RMvis(freqs,RM):
    return np.exp(2.j*RM*(0.3/freqs)**2).real

Nchan = 1000
freqs = np.linspace(0.12,0.18,Nchan)
RMs = [0,3,10,30]

def k3pk(T,k,fq,B):
    z = pspec.f2z(fq)
    bm = 0.76
    scalar = pspec.X2Y(z) * k**3 / (2*np.pi**2* B*1e9 * bm)
    return scalar * 1e6 *np.abs(T)**2

def rebin_log(x,y,bin=10):
    logx = np.log(np.abs(x))
    hist1,bins = np.histogram(logx,bins=bin,weights=y)
    hist2,bins = np.histogram(logx,bins=bin)
    logx = 0.5 * (bins[1:] + bins[:-1])
    return np.e**logx,hist1/np.where(hist2 == 0, 1., hist2)

H0 = 71.
c = 3e5
Om = 0.27

def H(z): return H0 * np.sqrt(Om) * (1.+z)**1.5
def nu2l2(f): return (0.3/f)**2
def k_bad(fq,RM):
    _z = pspec.f2z(fq)
    k = nu2l2(fq) * RM * 4.* H(_z) / ( (1+_z) * c)
    print k / 0.71
    return k / 0.71

fig = figure()
colors='kcmb'
for j,RM in enumerate(RMs):
    Nann = 3
    for i in range(Nann):
        if i == 0 or i+1 == Nann: continue
        #if i%4 != 0: continue
        nchan_i = int(Nchan/Nann)
        ntaps = 3
        ch1,ch2 = i*nchan_i,(i+1)*nchan_i
        ch_in,ch_out = ch1-nchan_i, ch2+nchan_i
        freqs_i = freqs[ch_in:ch_out]
        dnu = freqs_i[1] - freqs_i[0]
        eta = np.fft.fftshift(np.fft.fftfreq(nchan_i,dnu))
        z_mid = np.median(pspec.f2z(freqs_i))
        k = np.abs(eta*pspec.dk_deta(z_mid))

        V_nu = RMvis(freqs_i,RM)
        V_nu *= pspec.jy2T(freqs_i)
        Trms = pfb.pfb(V_nu,taps=ntaps,window=window,fft=np.fft.ifft)
        Trms = np.fft.fftshift(Trms)
        DELTA = k3pk(Trms,k=k,fq=np.mean(freqs_i),B=0.1/Nann)


        valid = k!=0.

        DELTA = np.extract(valid,DELTA)
        DELTA /= DELTA.sum()
        k = np.extract(valid,k)

        _k,_DELTA = rebin_log(k,DELTA)

        ax = fig.add_subplot(len(RMs),1,j+1)
        #ax.step(k,DELTA,where='mid',color='k',label='z = %2.1f'%z_mid)
        ax.step(_k,_DELTA,where='mid',color=colors[j],label='z = %2.1f'%z_mid)

        ax.vlines(k_bad(np.mean(freqs_i),RM),1e-14,1e0,color='0.5')

        ax.set_yscale('log')
        ax.set_xscale('log')

        ax.set_ylabel('RM=%2.1f $m^{-2}$'%RM)
        ax.set_ylim(1e-14,1e0)
        ax.set_yticks([1e-1,1e-4,1e-7,1e-10,1e-13])

        #if j == len(RMs)-1: legend(loc='lower right')

        if j+1 != len(RMs): ax.set_xticks([])
        else: ax.set_xlabel("$k\\ [h\mathrm{Mpc}^{-1}]$",fontsize=15)

fig.subplots_adjust(bottom=0.1)
fig.savefig('Pspec/figures/SingleSourceSpec.eps',fmt='eps')

show()
