#! /usr/bin/env python

import numpy as np
from capo import pspec
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axes_grid1 import host_subplot
from pylab import *

import sys

H0 = 100. #km/s/Mpc
c = 3e5 #km/s
Om = 0.27

def H(z): return H0 * np.sqrt(Om) * (1.+z)**1.5 #km/s/Mpc
def nu2l2(f): return (0.3/f)**2 #m^2

def RM_bad(fq,k):
    _z = pspec.f2z(fq)
    return k * c * (1.+_z) / ( 4. * H(_z) * nu2l2(fq))

fq = np.linspace(0.12,0.18,1000) #GHz

print RM_bad(0.15,0.16*0.71)

k0 = 0.20 * 0.71
k1 = 0.25 * 0.71 #Mpc^{-1}, not hMpc^{-1}
k2 = 0.30 * 0.71


fig = figure()

ax = host_subplot(111,axes_class=AA.Axes)

ax.fill_between(fq,RM_bad(fq,k0),RM_bad(fq,k2),facecolor='0.75')

line = ax.plot(fq,RM_bad(fq,k1))[0]
line.set_color('k')
line.set_linewidth(2.0)

ax.set_ylabel('$\mathrm{Infecting\\ Rotation\\ Measure}\\ [m^{-2}]$')
ax.set_xlabel('$\mathrm{Frequency\\ [GHz]}$')
ax.axis["bottom"].label.set(size=15)
ax.axis["left"].label.set(size=15)

zs = pspec.z2f(np.array([7.,7.5,8.,8.5,9.,9.5,10.,10.5]))
ax_z = ax.twin()
ax_z.set_xticks(zs)
ax_z.set_xticklabels(['7','7.5','8','8.5','9','9,5','10','z=10.5'])
ax_z.axis["right"].major_ticklabels.set_visible(False)

xlim([fq[0],fq[-1]])

fig.savefig('Pspec/figures/bad_rm.eps',fmt='eps')

show()
