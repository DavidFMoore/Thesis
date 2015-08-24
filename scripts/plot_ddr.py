#! /usr/bin/env python

import numpy as np
from pylab import *
import aipy as a

lat=-31.*np.pi/180.

SIZE=500
m,l = np.indices((SIZE,SIZE))
m = 2.*(m/float(SIZE)-0.5)
l = 2.*(l/float(SIZE)-0.5)
z = np.sqrt(1. - m**2 - l**2)
mask = np.where(z > 0, 0, 1)

xyz = np.array([l,m,z])
top_m = a.coord.top2eq_m(0.,lat)
radec = np.array([[np.dot(top_m,xyz[:,i,j]) for i in range(SIZE)] for j in range(SIZE)])
ra,dec = a.coord.eq2radec(radec.transpose([2,1,0]))
delay = (32./0.3)*np.cos(dec)*np.sin(ra)
fringe = 1000.*(32./0.3)*(1./a.const.sidereal_day)*np.cos(ra)*np.cos(dec)

def plot_coords():
    o = Circle((0.,0.), radius=1., fill=False)
    gca().add_patch(o)
    top_m = a.coord.eq2top_m(0.,lat)
    for _ra in np.linspace(0,2.*np.pi,9)[:-1]:
        _dec = np.linspace(-np.pi/2.,np.pi/2.,100)
        _ra = np.ones_like(_dec)*_ra
        top = np.dot(top_m,a.coord.radec2eq([_ra,_dec]))
        is_up = np.where(top[2] > 0, 0, 1)
        x = np.ma.array(top[0], mask=is_up)
        y = np.ma.array(top[1], mask=is_up)
        #x = (SIZE/2.)*(x + 1.)
        #y = (SIZE/2.)*(y + 1.)
        plot(x,y,'k:')
    for _dec in np.linspace(-np.pi/2.,np.pi/2.,7)[:-1]:
        _ra = np.linspace(0.,2.*np.pi, 100)
        _dec = np.ones_like(_ra)*_dec
        top = np.dot(top_m, a.coord.radec2eq([_ra,_dec]))
        is_up = np.where(top[2] > 0, 0, 1)
        x = np.ma.array(top[0], mask=is_up)
        y = np.ma.array(top[1], mask=is_up)
        #x = (SIZE/2.)*(x + 1.)
        #y = (SIZE/2.)*(y + 1.)
        plot(x,y,'k:')

def fmt():
    plot_coords()
    colorbar()
    xticks([])
    yticks([])
    fig.patch.set_visible(False)
    gca().axis('off')

fig = figure(figsize=(8,4))
subplot(121)
imshow(delay,
        aspect='auto',
        extent=[-1,1,-1,1],
        origin=(0,0))
fmt()
subplot(122)
imshow(fringe,
        aspect='auto',
        extent=[-1,1,-1,1],
        origin=(0,0))
fmt()

savefig('introduction/Interferometry/figures/DelayDelayRate.eps',fmt='eps')
show()
