#coding: utf8
#File cddiscretize.py
#José Amoreira
#July 2016

#Computation of the coefficients of the algebraic equations in standard form
#using any of the lds, the upwind, the hybrid or the ais schemes

import numpy as np

#Different shemes implemented
def lds(mesh, FD, srcCoeffs, bdrVals):
    xc, xf = mesh
    F, D = FD
    M, N = srcCoeffs
    phib_W, phib_E = bdrVals
    aW =  D[1:-1] + F[1:-1] * (xc[1:]  - xf[1:-1]) / (xc[1:] - xc[:-1])
    aW = np.insert(aW, 0, 0.0)
    aE = D[1:-1] - F[1:-1] * (xf[1:-1] - xc[:-1]) / (xc[1:] - xc[:-1])
    aE = np.append(aE,0.0)
    
    S = M
    S[0]  -= D[0]  + F[0]
    S[-1] -= D[-1] - F[-1]

    b = N
    b[0]  += (D[0]  + F[0])  * phib_W
    b[-1] += (D[-1] - F[-1]) * phib_E
    
    aP = aE + aW + F[1:]-F[:-1] - S

    return (aP, aW, aE, b)


def uw(mesh, F, D, srcCoeffs, bdrVals, coeffs):
    pass
    #return (aP, aW, aE, b)


def hyb(mesh, F, D, srcCoeffs, bdrVals, coeffs):
    pass
    #return (aP, aW, aE, b)


def ais(mesh, F, D, srcCoeffs, bdrVals, coeffs):
    pass
    #return (aP, aW, aE, b)


schemeFunc = {"lds": lds, "uw": uw, "hyb": hyb, "ais": ais}


def stdEqCoeffs(scheme, mesh, FD, srcCoeffs, bdrVals):
    """mkcoeffs(mesh, FD, scrCoeffs, bdrVals).
       Computes the standard form coefficients for diffusion with given
       convection 1D steady-state problem.
       scheme: (string, one of "lds", "uw", "hyb", "ais") scheme used
       FD: (F, D)
       F and D: (ndarrays) convection and diffusion strength coefficients
       srcCoeffs: (tuple) M and N arrays representing the source
       bdrVals: (tuple) lower and upper boudary values.
    """
    
    aP, aW, aE, b = schemeFunc[scheme](mesh, FD, srcCoeffs, bdrVals)

    return aP, aW, aE, b
