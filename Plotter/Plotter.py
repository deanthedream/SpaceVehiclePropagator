#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import os
import matplotlib.pyplot as plt


"""
Plot XY of XYZ data as line a 2D plot
args:
    xyz - (ndarray) a 3xn array containing x, y, and z position data
    num - int the figure number
"""
def plotXY(xyz, num, outpath, fname):
    fig = plt.figure(num=num)
    ax=fig.add_axes([.1,.1,.8,.8]) # This is the background axis
    plt.axis('equal')
    plt.plot(xyz[0],xyz[1],color='purple')

    plt.xlabel('x in (km)')
    plt.ylabel('y in (km)')

    plt.savefig(os.path.join(outpath,fname + '.png'), format='png', dpi=300, pad_inches=0.5, bbox_inches="tight")
    plt.close('all')

"""
Plot XZ of XYZ data as line a 2D plot
args:
    xyz - (ndarray) a 3xn array containing x, y, and z position data
    num - int the figure number
"""
def plotXZ(xyz, num, outpath, fname):
    plt.figure(num=num)
    plt.axis('equal')
    plt.plot(xyz[0],xyz[2],color='purple')

    plt.xlabel('x in (km)')
    plt.ylabel('z in (km)')

    plt.savefig(os.path.join(outpath,fname + '.png'), format='png', dpi=300, pad_inches=0.3, bbox_inches="tight")
    plt.close('all')



