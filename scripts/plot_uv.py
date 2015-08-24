#! /usr/bin/env python

import numpy as np
import aipy as a
from pylab import *

fq = np.linspace(0.1,0.2,203)[110:149]
aa = a.cal.get_aa('psa898_v003', fq)

fig = figure(figsize=(8,6))

#plot antenna layout
ax = fig.add_subplot(211)
rot_m = a.coord.eq2top_m(0., aa.lat)
for i in range(32):
    x,y,z = np.dot(rot_m, aa[i].pos - aa[0].pos)
    plot(x*0.3,y*0.3,'ks')
ax.set_xlim([-10,220])
ax.grid(True)
ax.set_xlabel(r'$\rm{Easting}\ [m]$', size=16)
ax.set_ylabel(r'$\rm{Northing}\ [m]$', size=16)

#plot uv plane.

def plot_bl(ax, i, j, cyan=False):
    global aa
    u,v,w = aa.gen_uvw(i,j)
    u = u.flatten()
    v = v.flatten()
    if not cyan:
        ax.plot(u,v,'k.')
    else:
        ax.plot(u,v,'c.')


ax = fig.add_subplot(212)
#plot concentric circles
for r in [2,5,10,20,50,100]:
    th = np.arange(0, 2.*np.pi+0.2, 0.1)
    x,y = r*np.cos(th), r*np.sin(th)
    ax.plot(x,y, color='0.5', linestyle='-')

for i in range(32):
    for j in range(32):
        if i == j:
            continue
        if i/4 == j/4:
            continue
        plot_bl(ax, i, j, cyan=False)

plot_bl(ax, 0, 16, cyan=True)
plot_bl(ax, 0, 17, cyan=True)
plot_bl(ax, 1, 16, cyan=True)
plot_bl(ax, 16, 0, cyan=True)
plot_bl(ax, 16, 1, cyan=True)
plot_bl(ax, 17, 0, cyan=True)

ax.set_ylim([-10,10])
ax.set_ylabel(r'$v\ [\lambda]$', size=16)
ax.set_xlabel(r'$u\ [\lambda]$', size=16)
ax.grid(True)


fig.subplots_adjust(hspace=0.3)
fig.savefig('Pspec/figures/UVCoverage.eps', fmt='eps')
show()
