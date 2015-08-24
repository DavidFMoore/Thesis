#! /usr/bin/env python

import numpy as np
import aipy as a
from pylab import *

fq = np.linspace(0.09,0.21,203)
aa = a.cal.get_aa('psa898_v003', fq)
aa[0].set_active_pol('x')
pb = aa[0].passband().real

pb /= np.max(pb)

plot(fq, 10.*np.log10(pb),'k',lw=2)
#plot(fq, pb,'k',lw=2)
ylim([-5.5,0.5])
xlim([0.095,0.205])

ylabel(r'$g(\nu)\ [{\rm dB}]$', size=16)
xlabel(r'$\nu\ [{\rm GHz}]$', size=16)

grid(True)

savefig('introduction/figures/BandPass.eps', fmt='eps')
show()
