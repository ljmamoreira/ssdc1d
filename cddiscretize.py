#coding: utf8
#File cddiscretize.py
#José Amoreira
#July 2016

#Computation of the coefficients of the algebraic equations in standard form
#using any of the lds, the upwind, the hybrid or the ais schemes

import numpy as np

#Linear interpolation coefficients, computed at the inner faces (N-1 elements)
#At the W boundary fwW=1; at the E boundary feE=1.
def lic(mesh):
    xc, xf = mesh
    fwW = (xc[1:] - xf[1:-1]) / (xc[1:] - xc[:-1])
    feE = (xf[1:-1] - xc[:-1]) / (xc[1:] - xc[:-1])
    return fwW, feE

#Convected interpolatios coefficients, computed at faces inner (N-1 elements)
#At the W boundary gwW=1; at the E boundary geE=1
def cic(P):
    xsi = 0.1
    f = P[1:-1] / (xsi + np.abs(P[1:-1]))
    hwW = 0.5 * (1.0 + f)  #needed for CVs 1:
    heE = 0.5 * (1.0 - f)  #needed for CVs :-1
    return hwW, heE


#Different schemes implemented
def lds(mesh, FD, srcCoeffs, bdrVals):
    F, D = FD
    M, N = srcCoeffs
    phib_W, phib_E = bdrVals
    fwW, feE = lic(mesh)

    aW = D[1:-1] + F[1:-1] * fwW 
    aW = np.insert(aW, 0, 0.0)
    aE = D[1:-1] - F[1:-1] * feE
    aE = np.append(aE, 0.0)
    
    S = M
    b = N
    S[0]  -= D[0]  + F[0]
    b[0]  += (D[0]  + F[0])  * phib_W

    S[-1] -= D[-1] - F[-1]
    b[-1] += (D[-1] - F[-1]) * phib_E
    
    aP = aE + aW + F[1:]-F[:-1] - S

    return (aP, aW, aE, b)


def uws(mesh, FD, srcCoeffs, bdrVals):
    F, D = FD
    M, N = srcCoeffs 
    phib_W, phib_E =bdrVals

    aW =  D[1:-1] + np.maximum(0.0,  F[1:-1])
    aW = np.insert(aW, 0, 0.0)
    aE = D[1:-1] + np.maximum(0.0, -F[1:-1])
    aE = np.append(aE, 0.0)
    
    S = M
    b = N
    S[0]  -= D[0]  + np.maximum(0.0,  F[0])
    b[0]  += (D[0]  + np.maximum(0.0, F[0] )) * phib_W

    S[-1] -= D[-1] + np.maximum(0.0, -F[-1])
    b[-1] += (D[-1] + np.maximum(0.0, F[-1])) * phib_E

    aP = aW + aE + F[1:] - F[:-1] - S

    return (aP, aW, aE, b)


def hyb(mesh, FD, srcCoeffs, bdrVals):
    F, D = FD
    ncvs = len(F)-1
    M, N = srcCoeffs 
    phib_W, phib_E = bdrVals
    fwW, feE = lic(mesh)

    aW = np.amax(
           np.column_stack((np.zeros(ncvs-1), D[1:-1] + F[1:-1]*fwW, F[1:-1])),
         axis=1)
    aW = np.insert(aW, 0, 0.0)
    aE = np.amax(
           np.column_stack((-F[1:-1], D[1:-1] - F[1:-1]*feE, np.zeros(ncvs-1))),
            axis=1)
    aE = np.append(aE, 0.0)
    S = M
    b = N
    Xw = D[0]
    if F[0] > -D[0]:
        Xw += F[0]
    S[0] -= Xw
    b[0] += Xw * phib_W

    Xe = D[-1]
    if F[-1] < D[-1]:
        Xe -= F[-1]
    S[-1] -= Xe
    b[-1] += Xe * phib_E

    aP = aW + aE + F[1:] - F[:-1] - S

    return (aP, aW, aE, b)


def cai(mesh, FD, srcCoeffs, bdrVals):
    F, D = FD
    M, N = srcCoeffs
    phib_W, phib_E = bdrVals
    fwW, feE = lic(mesh)
    hwW, heE = cic(F/D)

    geE = feE * heE / (1- feE - heE + 2*feE*heE)
    gwW = fwW * hwW / (1- fwW - hwW + 2*fwW*hwW)
    
    aW = D[1:-1] + F[1:-1] * gwW 
    aW = np.insert(aW, 0, 0.0)
    aE = D[1:-1] - F[1:-1] * geE
    aE = np.append(aE, 0.0)
    
    S = M
    b = N
    S[0]  -= D[0]  + F[0]
    b[0]  += (D[0]  + F[0])  * phib_W

    S[-1] -= D[-1] - F[-1]
    b[-1] += (D[-1] - F[-1]) * phib_E
    
    aP = aE + aW + F[1:]-F[:-1] - S

    return (aP, aW, aE, b)


schemeFunc = {"lds": lds, "uws": uws, "hyb": hyb, "cai": cai}


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


