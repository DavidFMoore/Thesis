#! /usr/bin/env python

import numpy as np
import aipy as a
from pylab import *

lat = -31. * (np.pi/180.)

ras = np.linspace(0,2*np.pi,20)[:-1]
decs = np.array([60,30,0,-30,-60])*np.pi/180.
ra_dec = np.array([[[ra,dec] for ra in ras] for dec in decs])
ra_dec = ra_dec.transpose([2,0,1])
XYZ = a.coord.radec2eq(ra_dec)
XYZ.shape = (3,XYZ.shape[1]*XYZ.shape[2])
ra_dec.shape = (2, ra_dec.shape[1]*ra_dec.shape[2])
m = a.coord.eq2top_m(0,lat)
xyz = np.array([np.dot(m,XYZ[:,i]) for i in range(XYZ.shape[-1])])
is_up = np.where(xyz[:,2] > 0, 1, 0)

x = xyz[:,0].copy()#.take(is_up)
y = xyz[:,1].copy()#.take(is_up)
ra = ra_dec[0,:].copy()
dec = ra_dec[1,:].copy()

x = x.compress(is_up)
y = y.compress(is_up)
ra = ra.compress(is_up)
dec = dec.compress(is_up)

def parang(H,d):
    global lat
    up = np.cos(lat)*np.sin(H)
    down = np.sin(lat)*np.cos(d) - np.cos(lat)*np.sin(d)*np.cos(H)
    return np.arctan2(up,down)

def parrot(H,d,unit=np.array([1,0])):
    unit /= np.sqrt(np.dot(unit.T,unit))
    psi = -1.*parang(H,d)
    m = np.array([[np.cos(psi), -np.sin(psi)],[ np.sin(psi),np.cos(psi)]])
    return np.array([np.dot(m[:,:,i],unit) for i in range(m.shape[-1])]).T

o = Circle((0.,0.),
        radius=1,
        fill=False)
#print lines of constant RA,dec
top_m = a.coord.eq2top_m(0.,lat)
for _ra in np.linspace(0,2.*np.pi,9)[:-1]:
    _dec = np.linspace(-np.pi/2.,np.pi/2.,100)
    _ra = np.ones_like(_dec)*_ra
    top = np.dot(top_m,a.coord.radec2eq([_ra,_dec]))
    is_up = np.where(top[2] > 0, 0, 1)
    plot(np.ma.array(top[0], mask=is_up),np.ma.array(top[1], mask=is_up),'k:')
for _dec in np.linspace(-np.pi/2.,np.pi/2.,7)[:-1]:
    _ra = np.linspace(0.,2.*np.pi, 100)
    _dec = np.ones_like(_ra)*_dec
    top = np.dot(top_m, a.coord.radec2eq([_ra,_dec]))
    is_up = np.where(top[2] > 0, 0, 1)
    plot(np.ma.array(top[0], mask=is_up),np.ma.array(top[1], mask=is_up),'k:')

vec1 = parrot(ra,dec,unit=np.array([0.,1.]))
vec2 = parrot(ra,dec,unit=np.array([1.,0.]))
for _x,_y,_v1,_v2 in zip(x,y,vec1.T,vec2.T):
    _v1 *= 1e-1
    _v2 *= 1e-1
    plot([_x-_v1[0]/2., _x+_v1[0]/2.], [_y-_v1[1]/2., _y+_v1[1]/2.],
        color='k',
        linewidth=2)
    plot([_x-_v2[0]/2., _x+_v2[0]/2.], [_y-_v2[1]/2., _y+_v2[1]/2.],
        color='k',
        linewidth=2)
gca().add_patch(o)
xlim([-1.1,1.1])
ylim([-1.1,1.1])
xticks([])
yticks([])
xlabel(r'$\rm{East}$', size=18)
ylabel(r'$\rm{North}$', size=18)
gca().axis('off')

savefig('introduction/Interferometry/figures/ParallacticAngle.eps', fmt='eps')
show()
