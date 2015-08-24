#! /usr/bin/env python

import numpy as np
from pylab import *

fq = np.linspace(0.1,0.2,1000)
colors = 'kcmb'

def gen_rm_spec(fq,rm):
    return 0.4*np.exp(2.j*rm*(0.3/fq)**2)

for i,rm in enumerate([3,10,30,100]):
    spec = gen_rm_spec(fq,rm).real
    plot(fq,spec.real + 3-i,colors[i])

legend()
yticks([])
ylabel(r'$Q\ [arbitrary\ units]$')
xlabel(r'$\rm{Frequency}\ [\rm{GHz}]$')
#savefig('figs/figure2.eps',fmt='eps')

show()
