#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np

"""
Keplerian Orbital Elements to XYZ Converstion
args:
    a - (ndarray) semi-major axis shape n
    e - (ndarray) eccentricity shape n
    i - (ndarray) inclination shape n
    W - (ndarray) longitude of the ascending node shape n
    w - (ndarray) argument of periapsis shape n
    v - (ndarray) true anomaly shape n
returns:
    xyz - (ndarray) 3xn cartesian coordinates
"""
def XYZfromKOE(a,e,i,W,w,v):
    xyz = np.zeros((3))
    r = a*(1-e**2)/(1+e*np.cos(v))
    xyz[0] = r*(np.cos(W)*np.cos(w+v)-np.sin(W)*np.sin(w+v)*np.cos(i))
    xyz[1] = r*(np.sin(W)*np.cos(w+v)+np.cos(W)*np.sin(w+v)*np.cos(i))
    xyz[2] = r*np.sin(i)*np.sin(w+v)
    return xyz

def XdYdZdfromKOE(a,e,i,W,w,v,mu):
    # [v_PQW] [-np.sqrt(mu/p)*np.sin(v),np.sqrt(mu/p)*(e+np.cos(v)),0.]

    #[IJK/PQW]
    # [[np.cos(W)*np.cos(w)-np.sin(W)*np.sin(w)*np.cos(i),-np.cos(W)*np.sin(w)-np.sin(W)*np.cos(w)*np.cos(i), np.sin(W)*np.sin(i)],
    # [np.sin(W)*np.cos(w)+np.cos(W)*np.sin(w)*np.cos(i), -np.sin(W)*np.sin(w)+np.cos(W)*np.cos(w)*np.cos(i), -np.cos(W)*np.sin(i)],
    # [np.sin(w)*np.sin(i),                               np.cos(w)*np.sin(i),                                np.cos(i)]]

    xdydzd = np.zeros((3))
    xdydzd[0] = (np.cos(W)*np.cos(w)-np.sin(W)*np.sin(w)*np.cos(i))*-np.sqrt(mu/a)*np.sin(v)\
        + (-np.cos(W)*np.sin(w)-np.sin(W)*np.cos(w)*np.cos(i))*np.sqrt(mu/a)*(e+np.cos(v))
    xdydzd[1] = (np.sin(W)*np.cos(w)+np.cos(W)*np.sin(w)*np.cos(i))*-np.sqrt(mu/a)*np.sin(v)\
        + (-np.sin(W)*np.sin(w)+np.cos(W)*np.cos(w)*np.cos(i))*np.sqrt(mu/a)*(e+np.cos(v))
    xdydzd[2] = np.sin(w)*np.sin(i)*-np.sqrt(mu/a)*np.sin(v) + np.cos(w)*np.sin(i)*np.cos(w)*np.sin(i)
    return xdydzd

