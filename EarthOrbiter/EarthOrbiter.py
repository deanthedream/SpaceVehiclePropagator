#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
from KOE.KOE import XYZfromKOE
from KOE.KOE import XdYdZdfromKOE
from Plotter.Plotter import plotXY
from Plotter.Plotter import plotXZ
from Plotter.PlotterEarth import plotXY_withEarth
from Plotter.PlotterEarth import plotXZ_withEarth
import argparse


"""
Args:
    a semi-major axis in km
    e eccentricity
    i inclination in rad
    W longitude of the ascending node in rad
    w argument of periapsis in rad
    v0 true anomaly in rad
    mu 398600.4415 km^3 s^-2 from  Barycentric Coordinate Time 
        https://iau-a3.gitlab.io/NSFA/NSFA_cbe.html#GME2009
"""
class EarthOrbiter:
    def __init__(self, a, e, i, W, w, v0, mu=398600.4415, epsmult=4.0):
        self.a = a
        self.e = e
        self.i = i
        self.W = W
        self.w = w
        self.v = v0
        self.mu = mu
        self.epsmult = epsmult
        self.T = self.orbitalPeriod(self.a,self.mu)

        #  Compute State Space
        self.x = np.concatenate((XYZfromKOE(a,e,i,W,w,v0),\
            XdYdZdfromKOE(a,e,i,W,w,v0,mu)))
        r0 = self.x[:3]
        v0 = self.x[3:]

        # constants
        self.r0norm = np.sqrt(np.sum(r0 ** 2.0, 0))  # ||r0||
        self.v0norm2 = np.sum(v0 * v0, 0)  # ||v0||^2
        self.nu0 = np.sum(r0 * v0, 0)  # r0 \cdot v0
        self.beta = 2 * self.mu / self.r0norm - self.v0norm2  # -2E
        self.alpha = self.beta / self.mu
        self.nu0osmu = self.nu0 / np.sqrt(self.mu)

    def orbitalPeriod(self,a,mu):
        return 2*np.pi*np.sqrt(a**3./mu)

    def takeStep(self, dt):
        Phi = self.calcSTM_vallado(dt)
        self.updateState(np.dot(Phi, self.x))


    def updateState(self, x):
        self.x = x
        r0 = self.x[:3]
        v0 = self.x[3:]

        # constants
        self.r0norm = np.sqrt(np.sum(r0 ** 2.0, 0))  # ||r0||
        self.v0norm2 = np.sum(v0 * v0, 0)  # ||v0||^2
        self.nu0 = np.sum(r0 * v0, 0)  # r0 \cdot v0
        self.beta = 2 * self.mu / self.r0norm - self.v0norm2  # -2E
        self.alpha = self.beta / self.mu
        self.nu0osmu = self.nu0 / np.sqrt(self.mu)


    """
    Algorithm 8 from Vallado Pg 93
    """
    def calcSTM_vallado(self, dt):
        # classify orbits
        epsval = 1e-12  #tolerance for classifciation

        eorbs = self.alpha >= epsval  # circle/ellipse
        porbs = np.abs(self.alpha) < epsval  # parabola
        horbs = self.alpha <= -epsval  #  Hyperbola

        #  initialize positions
        xi = 0.

        #  ellipse
        if self.alpha >= epsval:  # circle/ellipse
            atmp = self.alpha
            tmp = np.sqrt(self.mu) * dt * atmp
            circinds = np.abs(atmp - 1) > epsval
            if circinds:
                tmp *= 0.97

            xi = tmp

        #  parabola
        if np.abs(self.alpha) < epsval:  # parabola
            r = self.x[:3]
            v = self.x[3:]

            h = np.cross(r.T, v0.T).T
            p = np.sum(h * h, 0) / self.mu

            s = np.arctan2(1.0, (3.0 * np.sqrt(self.mu / p ** 3.0) * dt)) / 2.0
            w = np.arctan((np.tan(s)) ** (1.0 / 3.0))
            xi = np.sqrt(p) * 2.0 / np.tan(2 * w)
            self.alpha = 0

        #  hyperbola
        if self.alpha <= -epsval:  #  Hyperbola
            a = 1.0 / (self.alpha)
            xi = (
                np.sign(dt)
                * np.sqrt(-a)
                * np.log(
                    -2
                    * self.mu
                    * self.alpha
                    * dt
                    / (
                        self.nu0
                        + np.sign(dt)
                        * np.sqrt(-self.mu * self.alpha)
                        * (1.0 - self.r0norm * self.alpha)
                    )
                )
            )

        # loop
        counter = 0
        r = self.r0norm
        xiup = 10.0 * np.max(np.abs((xi, r)))
        while (
            np.max(np.abs(xiup))
            > self.epsmult * np.spacing(np.max(np.abs((xi, r))))
        ) and (counter < 1000):
            ps = xi ** 2.0 * self.alpha
            c2, c3 = self.psi2c2c3(ps)
            r = (
                xi ** 2.0 * c2
                + self.nu0osmu * xi * (1 - ps * c3)
                + self.r0norm * (1 - ps * c2)
            )
            xiup = (
                np.sqrt(self.mu) * dt
                - xi ** 3.0 * c3
                - self.nu0osmu * xi ** 2.0 * c2
                - self.r0norm * xi * (1 - ps * c3)
            ) / r
            xi += xiup
            counter += 1

        if counter == 1000:
            raise ValueError(
                "Failed to converge on xi: %e/%e"
                % (
                    np.max(np.abs(xiup)),
                    self.epsmult * np.spacing(np.max(np.abs((xi, r)))),
                )
            )

        # kepler solution
        f = 1.0 - xi ** 2.0 / self.r0norm * c2
        g = dt - xi ** 3.0 / np.sqrt(self.mu) * c3
        F = np.sqrt(self.mu) / r / self.r0norm * xi * (ps * c3 - 1.0)
        G = 1.0 - xi ** 2.0 / r * c2

        # Phi = np.zeros(([6] * 2))
        # st = 6
        # Phi[st : st + 6, st : st + 6] = np.vstack(
        #     (
        #         np.hstack((np.eye(3) * f, np.eye(3) * g)),
        #         np.hstack((np.eye(3) * F, np.eye(3) * G)),
        #     )
        # )

        Phi = np.vstack(
            (
                np.hstack((np.eye(3) * f, np.eye(3) * g)),
                np.hstack((np.eye(3) * F, np.eye(3) * G)),
            )
        )

        return Phi

    def psi2c2c3(self, psi0):

        c2 = 0.#np.zeros(len(psi0))
        c3 = 0.#np.zeros(len(psi0))

        psi12 = np.sqrt(np.abs(psi0))
        #pos = psi0 >= 0
        #neg = psi0 < 0
        if psi0 >= 0: #np.any(pos):
            c2 = (1 - np.cos(psi12)) / psi0
            c3 = (psi12 - np.sin(psi12)) / psi12 ** 3.0
        else: #psi0 < 0 if any(neg):
            c2 = (1 - np.cosh(psi12)) / psi0
            c3 = (np.sinh(psi12) - psi12) / psi12 ** 3.0

        tmp = c2 + c3 == 0
        if tmp == 0: #any(tmp):
            c2 = 1.0 / 2.0
            c3 = 1.0 / 6.0

        return c2, c3


def main(a=6371+400, e=0.02, i=0, W=0, w=0, v0=0, dt=60, mu=398600.4415):
    #  Contains a sample propagation
    EO = EarthOrbiter(a, e, i, W, w, v0, mu=398600.4415)

    n = np.ceil(EO.T/dt).astype(int)
    xyzs = np.zeros((n+1,6))
    xyzs[i] = EO.x
    for i in np.arange(n+1):
        EO.takeStep(dt)
        xyzs[i] = EO.x

    outpath = './Output'
    fname = 'XY'
    num=1
    plotXY(xyzs.T, num, outpath, fname)
    fname = 'YZ'
    num=2
    plotXZ(xyzs.T, num, outpath, fname)
    fname = 'XY_Earth'
    num=3
    plotXY_withEarth(xyzs.T, num, outpath, fname)
    fname = 'YZ_Earth'
    num=4
    plotXZ_withEarth(xyzs.T, num, outpath, fname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a ArcHydro schema')
    parser.add_argument('--a', metavar='path', required=True,
                        help='semi-major axis in km')
    parser.add_argument('--e', metavar='path', required=True,
                        help='eccentricity')
    parser.add_argument('--i', metavar='path', required=True,
                        help='inclination')
    parser.add_argument('--W', metavar='path', required=True,
                        help='longitude of ascending node')
    parser.add_argument('--w', metavar='path', required=True,
                        help='argument of periapsis')
    parser.add_argument('--v0', metavar='path', required=True,
                        help='True Anomaly in rad')
    parser.add_argument('--dt', metavar='path', required=False,
                        help='The Time Spacing in (s) to plot')
    parser.add_argument('--mu', metavar='path', required=False,
                        help='Gravitational Parameter in km^3s^-2')

    main(a, e, i, W, w, v0, dt, mu)

