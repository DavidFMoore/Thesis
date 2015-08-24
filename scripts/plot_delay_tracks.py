#! /usr/bin/env python

import numpy as np
import aipy as a
from pylab import *

fq = np.linspace(0.1,0.2,203)
d = np.load('data/psa919_0,31.npz')
times = d['times']
d = d['xx']
ker = np.where(d == 0, 0, 1)

delay = []
for i,_d in enumerate(d):
    w = a.dsp.gen_window(_d.shape[-1], window='blackman-harris')
    _d = np.fft.ifft(_d*w)
    _k = np.fft.ifft(ker[i]*w)
    _d, info = a.deconv.clean(_d, _k, tol=1e-4)
    _d += info['res']/a.img.beam_gain(_k)
    delay.append(np.fft.fftshift(_d))

delay = np.array(delay)

imshow(np.log10(np.abs(delay)**2),
        aspect='auto',
        origin=(0,0),
        interpolation='nearest'
        )
show()
