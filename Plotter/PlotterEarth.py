#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


"""
Plot XY of XYZ data as line a 2D plot with Earth
args:
    xyz - (ndarray) a 3xn array containing x, y, and z position data
    num - int the figure number
"""
def plotXY_withEarth(xyz, num, outpath, fname):
    fig = plt.figure()
    ax = plt.gca()
    plt.axis('equal')
    plt.xlabel('x in (km)')
    plt.ylabel('y in (km)')
    nus = np.linspace(start=0,stop=2.*np.pi,num=100,endpoint=True)
    m = Basemap(projection='ortho', resolution='l', lat_0=90, lon_0=0, ax=ax)
    m.bluemarble(scale=1);
    m.fillcontinents(color='k',lake_color='aqua')
    ax.plot(xyz[0]*1000+6371*1000,xyz[1]*1000+6371*1000,color='purple')
    ax.set_xlim([-(6371+1000)*1000+6371*1000,(6371+1000)*1000+6371*1000])
    ax.set_ylim([-(6371+1000)*1000+6371*1000,(6371+1000)*1000+6371*1000])

    plt.show(block=False)
    plt.savefig(os.path.join(outpath,fname + '.png'), format='png', dpi=300, pad_inches=0.1, bbox_inches="tight")
    plt.close('all')


"""
Plot XZ of XYZ data as line a 2D plot with Earth
args:
    xyz - (ndarray) a 3xn array containing x, y, and z position data
    num - int the figure number
"""
def plotXZ_withEarth(xyz, num, outpath, fname):
    fig = plt.figure()
    ax = plt.gca()
    plt.axis('equal')
    plt.xlabel('x in (km)')
    plt.ylabel('y in (km)')
    nus = np.linspace(start=0,stop=2.*np.pi,num=100,endpoint=True)
    m = Basemap(projection='ortho', resolution='l', lat_0=0, lon_0=-90, ax=ax)
    m.bluemarble(scale=1);
    m.fillcontinents(color='k',lake_color='aqua')
    ax.plot(xyz[0]*1000+6371*1000,xyz[2]*1000+6371*1000,color='purple')
    ax.set_xlim([-(6371+1000)*1000+6371*1000,(6371+1000)*1000+6371*1000])
    ax.set_ylim([-(6371+1000)*1000+6371*1000,(6371+1000)*1000+6371*1000])

    plt.show(block=False)
    plt.savefig(os.path.join(outpath,fname + '.png'), format='png', dpi=300, pad_inches=0.1, bbox_inches="tight")
    plt.close('all')


