#! /usr/bin/env python

import numpy as np
import aipy as a
import PolSim
from pylab import *
from healpy import *
import optparse

o = optparse.OptionParser()
o.add_option('-C','--cal',default='psa_null', help='Cal file')
o.add_option('-c','--chans',default='110_149', help='Channels to use')
o.add_option('-r','--remove',default=0, help='Number of sources to remove.')
opts, args = o.parse_args(__import__('sys').argv[1:])

Niter = 30

#freqs = np.linspace(0.1,0.2,203)[37:70]
#freqs = np.linspace(0.1,0.2,203)[85:110]
#freqs = np.linspace(0.1,0.2,203)[110:149]
cdown,cup = map(int, opts.chans.split('_'))
chans = slice(cdown,cup,1)
freqs = np.linspace(0.1,0.2,203)[chans]
bl_EW = 16.
N_zbins = 1
Nside = 64

#CONFIG_FILE = 'scripts/6C.cf'
CONFIG_FILE=args[0]
CAL_FILE = opts.cal
RM_DIST_FILE = '/home/damo/PolSim/input_data/Oppermann.copy.npz'
#RM_DIST_FILE = 'data/Oppermannx2.npz'

aa = a.cal.get_aa(CAL_FILE, freqs)

bm = PolSim.Beam(fromaipy=True, aa=aa).beam
D = PolSim.Dist(CONFIG_FILE, 0.76, RM_DIST_FILE)

def Trms_vs_fq(fq, jy_spec):
    #convert to mK
    Tspec = jy_spec * PolSim._pspec.jy2T(fq)
    fq0 = np.median(fq)
    z = PolSim._pspec.f2z(fq0)
    umag_fq0 = bl_EW * (fq0 / 0.150)
    k_perp = PolSim._pspec.dk_du(z) * umag_fq0
    etas = PolSim._pspec.f2eta(fq)
    k_par = PolSim._pspec.dk_deta(z) * etas
    k = np.sqrt(k_par**2 + k_perp**2)
    w = a.dsp.gen_window(len(Tspec), window='blackman-harris')
    Trms = np.fft.ifft(Tspec*w)
    return Trms, (k, k_par, k_perp)

r = Rotator(coord=['E','G'])
def hpm_gal2eq(N,_r):
    _nside = npix2nside(N)
    _px = _r(pix2vec(_nside, range(N)))
    return vec2pix(_nside,_px[0],_px[1],_px[2])
RMmap = read_map('data/faraday.fits')
RMmap = alm2map(map2alm(RMmap, lmax=int(2.*np.pi*16)), get_nside(RMmap))
RMmap = RMmap[hpm_gal2eq(len(RMmap), r)]
RMmap = ud_grade(RMmap, Nside)

I,Q = [],[]
for i in range(Niter):
    print "working on iteration %d/%d"%(i+1,Niter)
    SV = PolSim.SimVis(D, bm, bl_EW, freqs)
    #visXX = SV.SimVis(pol='xx', remove=int(opts.remove))
    #visYY = SV.SimVis(pol='yy', remove=int(opts.remove))
    visXX = 0
    visYY = 0
    for F,L,M,G,P in zip(SV.prms['F'], SV.prms['L'], SV.prms['M'], SV.prms['G'], SV.prms['P']):
        __I = SV.srcspec(F,L,G)
        z = np.sqrt(1.-L**2-M**2)
        RM = RMmap[vec2pix(Nside, L, M, z)]
        __Q = P*SV.RMphs(RM,0.)
        visXX += SV.beam.Response(L, M, SV.fqs, 'xx') * (1.+__Q)*__I
        visYY += SV.beam.Response(L, M, SV.fqs, 'yy') * (1.-__Q)*__I

    _I = visXX + visYY
    _Q = visXX - visYY

    _I, k = Trms_vs_fq(freqs, _I)
    _Q, k = Trms_vs_fq(freqs, _Q)

    w = a.dsp.gen_window(len(freqs), window='blackman-harris')
    neq = np.sum(w)**2/np.sum(w**2)
    BW = neq*(freqs[1]-freqs[0])
    fq0 = np.median(freqs)
    _I = PolSim._pspec.k3pk_from_Trms([_I],[1.], k=k[0], fq=fq0, B=BW)[0]
    _Q = PolSim._pspec.k3pk_from_Trms([_Q],[1.], k=k[0], fq=fq0, B=BW)[0]

    I.append(_I)
    Q.append(_Q)

I = np.array(I)
Q = np.array(Q)

#I = np.median(I, axis=0)
#Q = np.median(Q, axis=0)

np.savez('data/corr_%.2f.npz'%(PolSim._pspec.f2z(np.median(freqs))), I=I, Q=Q, k=k)
