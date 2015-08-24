#! /usr/bin/env python

import numpy as np
from pylab import *


psa = np.load('data/flags_psa.npz')['x']
fq = np.linspace(0.1,0.2,psa.shape[-1])
fill_between(fq, np.zeros_like(psa), psa, color='0.5')

pgb = np.load('data/flags_pgb.npz')['x']
sfreq = 0.121143578125
sdf = 7.32421875e-05
fq = np.linspace(sfreq,sfreq + 1024*sdf,1024)
fill_between(fq, np.zeros_like(pgb), pgb, color='0.75')

ylabel(r'$\rm{Latency}$', size=16)
xlabel(r'$\rm{Frequency}\ [\rm{GHz}]$', size=16)
savefig('introduction/figures/RFI.eps', fmt='eps')
show()
