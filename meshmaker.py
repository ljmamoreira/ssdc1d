#coding: utf8
#File mesh.py
#José Amoreira
#July 2016

#uniform, centered, 1D mesh maker

import numpy as np

def ucmesh(n,a,b):
    """ucmesh(n, xmin, xmax). Generates a uniform centered mesh with n cells,
       defined by np.arrays xc[n] storing the center cell positions, and xf[n+1]
       storing the cell faces positions.
    """
    assert isinstance(n,int), "mkmesh: n is not an integer."
    xf = np.linspace(a, b, n+1)
    xc = (xf[1:] + xf[:-1])/2.0

    mesh = (xc, xf)
    return mesh


