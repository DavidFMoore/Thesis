#! /usr/bin/env python

import numpy as np
import pyfits as pf
from healpy import *
from pylab import *

m = read_map('/home/damo/Haslam/data/faraday.fits', hdu=3)
#m = alm2map(map2alm(m, lmax=int(2.*15.*np.pi)), get_nside(m))

mollview(m,
        title="",
        unit=r'${\rm Faraday\ Depth}\ [{\rm m}^{-2}]$',
        min=-140.,
        max=140.,
        )
graticule(dpar=30., dmer=30., coord='G')
savefig('Pspec/figures/Oppermann.eps', fmt='eps')
show()
