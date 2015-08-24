#! /usr/bin/env python

import numpy as np
import aipy as a
from capo.pspec import jy2T
from PolSim import *
from healpy import *
from pylab import *

Niter = 100
Nside = 64
aa = a.cal.get_aa('psa_null', np.array([0.165]))

D = Distributions.Dist('scripts/6C.cf', 0.76, 'data/Oppermann.copy.npz')

jy2T = jy2T(0.165)
Cell = []
for i in range(Niter):
    print "Working on iteration %d/%d"%(i+1,Niter)
    print "\t with %d sources"%D.Nsrc
    prms = D.SimSkyParams()
    x,y = prms['L'], prms['M']
    z = np.sqrt(1.-x**2 - y**2)

    m = np.zeros(nside2npix(Nside))
    for pi,fi,xi,yi,zi in zip(prms['P'],prms['F'],x,y,z):
        m[vec2pix(Nside, xi, yi, zi)] += jy2T*fi*pi
    Cell.append(anafast(m, lmax=1000))
Cell = np.array(Cell)
dCell = np.std(Cell, axis=0)
Cell = np.mean(Cell, axis=0)

#plot my Cell
plot(Cell,'b')
fill_between(range(len(Cell)), Cell-2.*dCell, Cell+2.*dCell, color='c', alpha=0.5)

#plot Bernardi points.
B = np.loadtxt('data/Bernardi_fig20.txt', delimiter=',')
ell_B = B[:,0]
Cl_B  = B[:,1]
E_B   = np.abs(B[:,2] - B[:,1])
plot(ell_B, Cl_B, 'kx')
errorbar(ell_B, Cl_B, yerr=E_B, ecolor='k', fmt=None)


xlim([10**1.5, 10**3.5])
ylim([10**0, 10**2.5])
gca().set_xscale('log')
gca().set_yscale('log', nonposy='clip')
xlabel(r'$\ell$', fontsize=15)
ylabel(r'$C_\ell\ [{\rm mK}^2]$', fontsize=15)
grid(True)
savefig('Pspec/figures/SimCell.eps', fmt='eps')
show()
