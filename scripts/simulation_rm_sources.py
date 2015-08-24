#! /usr/bin/env python

import numpy as np
import PolSim

Niter = 1
freqs = np.linspace(0.1,0.2,203)[110:149]
bl_EW = 16.
RM_DIST_FILE = '/home/damo/PolSim/input_data/Opperman.copy.npz'

aa = a.cal.get_aa('psa_null', freqs)
bm = PolSim.Beam(fromaipy=True, aa=aa).beam
D = PolSim.Dist(__import__('sys').argv[1], 0.76, RM_DIST_FILE)

I,Q = [],[]
for i in range(Niter):
    print "Working on iteration %d/%d"%(i+1,Niter)
    SV = PolSim.SimVis(D, bm, bl_EW, freqs)
    visXX = SV.SimVis(pol='xx')
    visYY = SV.SimVis(pol='yy')
