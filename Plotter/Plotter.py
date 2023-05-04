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
    plt.figure(num=num)
    plt.axis('equal')
    plt.plot(xyz[0],xyz[1],color='purple')

    plt.xlabel('x in (km)')
    plt.ylabel('y in (km)')

    plt.savefig(os.path.join(outpath,fname), format='png', dpi=300)
    plt.close()

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

    plt.savefig(os.path.join(outpath,fname), format='png', dpi=300)
    plt.close()



