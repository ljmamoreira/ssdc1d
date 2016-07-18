#coding: utf8
#File lds.py
#José Amoreira
#July 2016

#Computation of the coefficients of the algebraic equations in standard form
#using the linear differencing scheme 

import numpy as np

def mkcoeffs(mesh, F, D, srcCoefs, bdrVals):
    xc, xf = mesh
    phib_W, phiB_E =bdrVals
    M, N = srcCoefs 

    aE = np.zeros_like(xc)
    aW = np.zeros_like(xc)
    aP = np.zeros_like(xc)
    b = np.zeros_like(xc)
    S = np.zeros_like(xc)

    aE[:-1] = D[1:-1] - F[1:-1] * (xf[1:-1] - xc[:-1]) / (xc[1:] - xc[:-1])
    aW[1:] =  D[1:-1] + F[1:-1] * (xc[1:]  - xf[1:-1]) / (xc[1:] - xc[:-1])

    S = M
    S[0]  -= D[0]  + F[0]
    S[-1] -= D[-1] - F[-1]

    b = N
    b[0]  += (D[0] + F[0]) * phib_W
    b[-1] += (D[-1] - F[-1]) * phib_E

    aP = aE + aW + F[1:] - F[-1:] - S
    return aP, aE, aW, b
