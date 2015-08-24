#! /usr/bin/env python

import numpy as np
from pylab import *
from capo import pspec

N = 100
fig = figure()
#plot a fringe of constant kperp

i,j = np.indices((N,N))
kperp = np.linspace(0.0,0.35,N)
fq = np.linspace(0.12,0.18,N)

fringe = np.cos(2.*np.pi*10.*kperp[j])

ax = fig.add_subplot(111)
imshow(fringe,
        extent=[kperp[0],kperp[-1],fq[0],fq[-1]],
        origin=(0,0),
        aspect='auto',
        cmap = cm.Blues,
        alpha=0.5
        )

for b in [30.,100.,200.,350.,500.]:
    b *= fq/0.3 #in lambda.
    b *= pspec.dk_du(7.5)
    plot(b, fq, 'k', lw=2)

ax.set_xlabel(r'$k_{||}\ [h{\rm Mpc}^{-1}]$', size=16)
ax.set_ylabel(r'${\rm Frequency}\ [{\rm GHz}]$', size=16)


fig.savefig('introduction/figures/BaselineTrack.eps', fmt='eps')
show()
