#! /usr/bin/env python

import aipy as a
import numpy as np
from pylab import *

N = 100

fq = np.array([0.127,0.164])
aa = a.cal.get_aa('psa_null', fq)

x = np.linspace(-1,1,N)
y = np.zeros_like(x)
z = np.sqrt(1-x**2-y**2)

aa[0].set_active_pol('x')
BMx = aa[0].bm_response((x,y,z))
alt = np.linspace(-90., 90., N)

fig = figure(figsize=(8,3))
c = ['k','c']
ax = fig.add_subplot(121)#, polar=True)
for i,fi in enumerate(fq):
    Prad = np.mean(BMx[i])
    dBi = 10.*np.log10(BMx[i]/Prad)
    ax.plot(alt, dBi, color=c[i], lw=2, label=r'$%0.3f\ {\rm GHz}$'%fi)

ax.set_xlabel(r'${\rm Altitude}\ [^\circ]$', size=16)
ax.set_ylabel(r'$A(\theta)\ [{\rm dBi}]$', size=16)
ax.set_xlim([-90,90])
ax.set_xticks([-90,-60,-30,0,30,60,90])
ax.grid(True)
ax.legend(loc='lower center')

y,x = np.indices((N,N)).astype(float)
x = 2.*x/float(N) - 1.
y = 2.*y/float(N) - 1.
z = np.sqrt(1.-x**2-y**2)
x = x.flatten()
y = y.flatten()
z = z.flatten()
fq = np.linspace(0.1,0.2,203)
aa = a.cal.get_aa('psa_null', fq)
aa[0].set_active_pol('x')
BMx = np.abs(aa[0].bm_response((x,y,z)))#**2
BMx /= np.abs(aa[0].bm_response((0,0,1)))#**2
BMx = np.mean(BMx, axis=1)

Aeff = (0.3/fq)**2 / BMx

ax = fig.add_subplot(122)
ax.plot(fq, Aeff, 'k', lw=2)
ax.grid(True)
ax.set_xlim([0.1,0.2])
ax.set_ylabel(r'$A_{eff}\ [{\rm m}^2]$', size=16)
ax.set_xlabel(r'$\nu\ [{\rm GHz}]$', size=16)

fig.subplots_adjust(left=0.125, right=0.9, wspace=0.3, bottom=0.2)
fig.savefig('introduction/figures/PAPERBeam.eps', fmt='eps')
show()
