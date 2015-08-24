#! /usr/bin/env python

import aipy as a
import numpy as np
from pylab import *

fq = np.array(np.linspace(0.1,0.2,203))
aa = a.cal.get_aa('psa_null', fq)

im = a.img.Img(size=100, res=0.5)
N = int(im.shape[0])
x,y,z = im.get_top(center=(N/2,N/2))
x = x.flatten()
y = y.flatten()
z = z.flatten()
valid = np.where(z<=0.,1,0.)

aa[0].set_active_pol('x')
BMx = np.abs(aa[0].bm_response((x,y,z)))**2
BMx = BMx.reshape((203,N,N))
aa[0].set_active_pol('y')
BMy = np.abs(aa[0].bm_response((x,y,z)))**2
BMy = BMy.reshape((203,N,N))

valid.reshape((N,N))
del im

magicfq = np.argmin(np.abs(fq-0.164))
def get_metric(x,y):
    Aplus  = np.abs(x + y)**2
    Aminus = np.abs(x - y)**2
    return np.sum(Aminus)/np.sum(Aplus)

metric = []
for i in range(203):
    if i == magicfq:
        plot_bm = figure(figsize=(8,4))
        subplot(121)
        imshow(np.ma.array(BMx[i,:,:] + BMy[i,:,:], mask=valid),
                aspect='auto',
                interpolation='nearest')
        colorbar()
        gca().axis('off')
        subplot(122)
        imshow(np.ma.array(BMx[i,:,:] - BMy[i,:,:], mask=valid),
                aspect='auto',
                interpolation='nearest')
        colorbar()
        gca().axis('off')
    metric.append(get_metric(BMx[i],BMy[i]))

metric_fig = figure()
plot(fq,metric,'k')
xlim([0.120,0.180])
ylim([0,6e-3])
grid(True)
ylabel(r'$\mathcal{A}_-/\mathcal{A}_+$', size=15)
xlabel(r'$\rm{Frequency}\ [\rm{GHz}]$', size=15)
print metric[magicfq]
print fq[np.argmin(metric)]

plot_bm.savefig('Pspec/figures/BeamPlusMinus.eps', fmt='eps')
metric_fig.savefig('Pspec/figures/BeamLeakage.eps', fmt='eps')

show()
