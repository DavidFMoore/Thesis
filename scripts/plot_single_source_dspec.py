#! /usr/bin/env python

import numpy as np
from aipy import dsp
from pylab import *

fq = np.linspace(0.1,0.2,1000)
colors = 'kcmb'

def gen_rm_spec(fq,rm):
    return np.exp(2.j*rm*(0.3/fq)**2)

x = [60., 128.,  300.,  865.]
y = [0.42, 0.27, 0.17, 0.107]
figure(figsize=(8,4))
for i,rm in enumerate([3,10,30,100]):

    subplot(121)
    spec = gen_rm_spec(fq,rm)
    plot(fq,0.4*spec.real + 3-i, colors[i])

    subplot(122)
    dspec = np.fft.fftshift(np.abs(np.fft.ifft(spec)))
    dlys = np.fft.fftshift(np.fft.fftfreq(1000,fq[1]-fq[0]))
    plot(dlys,dspec,colors[i],label='%d $m^{-2}$'%rm)
    text(x[i],y[i],"$%d\\ m^{-2}$"%rm,size='large')

subplot(121)
ylabel(r'$Q\ [\rm{arbitrary\ units}]$', fontsize=15)
xlabel(r'$\rm{Frequency}\ [\rm{GHz}]$', fontsize=15)
subplot(122)
yticks([])
ylabel(r'$|\widetilde{V}(u,v,\eta)|\ [\mathrm{arbitrary\ units}]$',fontsize=15)
xlabel('$\mathrm{Delay\ [ns]}$',fontsize=15)
xlim([-200,2000])

subplots_adjust(bottom=0.15)
savefig('Pspec/figures/SingleSource.eps',fmt='eps')
show()
