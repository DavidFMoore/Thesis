#! /usr/bin/env python

import numpy as np
from pylab import *

dat = np.loadtxt('data/Bernardi_fig20.txt',comments='#',delimiter=',')
B = {}
B['l'] = dat[:,0]
B['Cl'] = dat[:,1]
B['dC'] = np.abs(dat[:,2] - dat[:,1]) * 0.5 * B['l'] * (B['l']+1.)/np.pi
B['Cl2'] = 0.5*B['Cl']*B['l']*(B['l']+1.)/np.pi

dat = np.loadtxt('data/Pen_fig4.txt',comments='#',delimiter=',')
P = {}
P['l'] = dat[:,0]
P['Cl2'] = 1e6 * dat[:,1]
P['Cl'] = P['Cl2'] / (0.5*P['l']*(P['l']+1.)/np.pi)


dat = np.load('data/Haslam.npz')
l,cl = dat['l'],dat['cl']
cl = np.extract(np.where(l <= 200),cl) * 1e6
l = np.extract(np.where(l <= 200),l)
cl2 = cl * l*(l+1.)/(2.*np.pi)
Pfrac = 2e-2
cl2 *= 10**-5#Pfrac**2

H={}
H['l']=np.linspace(10**2,10**4,1000)
H['Cl2']=np.exp(np.polyval(np.polyfit(np.log(l),np.log(cl2),2),np.log(H['l'])))
H['Cl'] = 2.*np.pi*H['Cl2'] / (H['l']*(H['l']+1.))

PAPERl = 4.*np.pi*16.
print PAPERl

fig = figure()
ax = fig.add_subplot(111)

ax.plot(B['l'],B['Cl2'],'k.')
ax.errorbar(B['l'],B['Cl2'],yerr=[[min(B['Cl2'][i]-1e-9,B['dC'][i]) for i in range(len(B['Cl2']))],B['dC']],ecolor='k',fmt=None)

ax.plot(P['l'],P['Cl2'], 'm')
ax.plot(H['l'],H['Cl2'], 'b:')
ax.plot(l,     cl2     , 'b:')

ax.set_xlim(10**1.5,10**3.5)
ax.set_ylim(1e4,1e8)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylabel('$\ell(\ell+1)\\ C_\ell / 2\\pi\\ [mK^2]$', fontsize=15)
ax.set_xlabel('$\ell$',fontsize=15)
grid(True)

fig.savefig('Pspec/figures/CellCompare.eps')
show()
