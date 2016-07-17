#coding: utf8
#File mesh.py
#José Amoreira
#July 2016

import numpy as np

def mkmesh(n,a,b):
    """mkmesh(a,b,n). Generates a uniform centered mesh, defined by np.arrays
       xc[n] storing the center cell positions, and xf[n+1] storing the cell
       faces positions.
    """
    assert isinstance(n,int), "mkmesh: n is not an integer."
    xf = np.linspace(a, b, n+1)
    xc = (xf[1:] + xf[:-1])/2.0
    return xc, xf


