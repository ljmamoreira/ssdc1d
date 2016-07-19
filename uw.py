#coding: utf8
#File uw.py
#José Amoreira
#July 2016

#Computation of the coefficients of the algebraic equations in standard form
#using the upwind scheme 

import numpy as np

def mkcoeffs(mesh, F, D, srcCoeffs, bdrVals):
    """mkcoeffs(mesh, F, D, scrCoeffs, bdrVals).
       Computes the standard form coefficients for diffusion with given
       convection 1D steady-state problem using the LDS scheme.
       F and D: (ndarrays) convection and diffusion strength coefficients,
       srcCoeffs: (tuple) M and N arrays representing the source
       bdrVals: (tuple) lower and upper boudary values.
    """

    xc, xf = mesh
    phib_W, phiB_E =bdrVals
    M, N = srcCoeffs 

    aE = np.zeros_like(xc)
    aW = np.zeros_like(xc)
    aP = np.zeros_like(xc)
    b = np.zeros_like(xc)
    S = np.zeros_like(xc)
