#! /usr/bin/env python

import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm as N
from pylab import *

nu0 = 0.15
c = 0.3
Nmc = 300
Dnu = 10**np.linspace(-4,-1,100)

def inner_int(dnu,phi):
    _nu = np.random.uniform(nu0-dnu/2., nu0+dnu/2.,size=Nmc)
    return np.mean(np.exp(-2.j*phi*c**2/_nu**2))

pdf = np.load('data/Oppermann.copy.npz')
Phi = pdf['arr_0']
pdf = pdf['arr_1']

invF = interp1d(pdf,Phi)
#invF = N(scale=25.)

attrition = []
for dnu in Dnu:
    _att = []
    #Phi = sorted(invF.rvs(size=Nmc))
    Phi = [invF(x) for x in np.random.uniform(pdf.min(), pdf.max(), size=30*Nmc)]
    for phi in Phi:
        _att.append(inner_int(dnu, phi))
    attrition.append(np.mean(np.abs(_att)**2))
plot(1000.*Dnu, attrition,'k')

ylabel(r'${\rm Attenuation\ Factor}$', size=15)
xlabel(r'$\Delta\nu\ [MHz]$',size=15)
gca().set_xscale('log')
grid(True)

savefig('Pspec/figures/BandwidthDepolarization.eps', fmt='eps')
show()
