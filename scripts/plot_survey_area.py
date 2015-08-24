#! /usr/bin/env python

import numpy as np
import aipy as a
from healpy import *
from aipy.coord import *
from pylab import *

freq = 0.164
NSIDE=64
lat = -31.*np.pi/180.

def rot_bm(hpm, ha, dec):
    _m = eq2top_m(ha, dec)
    _xyz = pix2vec(NSIDE, range(len(hpm)))
    ncrd = np.dot(_m, _xyz)
    n_px = vec2pix(get_nside(hpm), ncrd[0], ncrd[1], ncrd[2])
    return hpm[n_px]

from bm_prms import prms
bm = prms['beam'](np.array([freq]), nside=NSIDE, lmax=20, mmax=20, deg=7)
bm.set_params(prms['bm_prms'])
px = range(nside2npix(NSIDE))
xyz = pix2vec(NSIDE, px)
poly = np.array([h.map[px] for h in bm.hmap])
Axx = np.polyval(poly, freq)
Axx = np.where(xyz[-1] >= 0, Axx, 0)
Axx = Axx*Axx
Axx *= np.max(Axx)/np.sum(Axx)

T = np.zeros_like(Axx)
t_int = 30.*Axx
files = __import__('glob').glob('/home/damo/Qspec/MyLST/cnt.0[1234567]*.uv')
print files
for f in files:
    print f
    uv = a.miriad.UV(f)
    a.scripting.uv_selector(uv, '0_16', 'xx')
    for (uvw,lst,(i,j)),d in uv.all():
        if lst >= 8.*np.pi/12.:
            continue
        if lst <= 1.*np.pi/12.:
            continue
        d = d[110:149]
        w = a.dsp.gen_window(len(d), window='blackman-harris')
        d = np.sum(d.real*w)/np.sum(w)
        T += rot_bm(d*t_int, lst, lat)

sr_per_pix = 4.*np.pi/float(len(T))

T /= 24.*3600. #t is in days now.
print "Total days: %.2f"%(np.sum(T)*4.*np.pi)
T /= sr_per_pix

FOV = np.sum(T)**2 / np.sum(T**2)
print "Field of view (sr): %.2f"%(FOV * sr_per_pix)

mollview(T*24.,
        unit=r'$t_{eff}\ [\rm{hr}\cdot\rm{sr}^{-1}]$',
        min=0.,
        max=150.,
        title='',
        )
graticule(dpar=30, dmer=30., coord='E')
gca().axis('off')
savefig('Pspec/figures/FOV.eps', fmt='eps')
show()
