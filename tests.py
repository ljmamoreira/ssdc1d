#coding: utf8
#File tests.py
#José Amoreira
#July 2016

#Solve Example 5.1 and following from Versteeg and Malalakesera, check the
#results from those presented there


import numpy as np
import meshmaker
import lds

#Define here a local version of init(), in order to keep the parameters

def init(mesh):
    """init(mesh). Defines the values for density, velocity, conduction
       coefficient, convection and diffusion strength parameters F and D
       source terms and boundary values at each CV.
    """
    xc, xf = mesh

    #Values below taken from Versteeg & Malalasekera Ex. 5.1b.
    #DO NOT CHANGE THE VALUES!

    #Density, velocity, convection strength parameter F:
    rho = 1.0
    v = 2.5
    F = np.ones_like(xf) * rho * v

    #Conduction coefficient and diffusion strength parameter D:
    gamma = 0.1 
    D = np.zeros_like(xf)
    D [1:-1] = gamma / (xc[1:] - xc[:-1])
    D[0] = gamma / (xc[0] - xf[0])
    D[-1] = gamma / (xf[-1] - xc[-1])

    #Source terms
    M = np.zeros_like(xc)
    N = np.zeros_like(xc)
    srcCoeffs = (M, N)

    #Boundary values
    phib_W = 1.0
    phib_E = 0.0
    bdrVals = (phib_W, phib_E)

    return F, D, srcCoeffs, bdrVals

#Small auxiliary function: returns True if arrays are very close
def isSmallDiff(a1,a2):
    d = a1 - a2
    return np.dot(d,d) < np.finfo(float).eps

#Small auxiliary function 
def okIfTrue(ok):
    if ok:
        return "OK."
    return "Wrong."


def test_mkmesh():
    mesh = meshmaker.mkmesh(5, 0.0, 1.0)
    xc, xf = mesh
    centerOK = isSmallDiff(xc, np.array([0.1, 0.3, 0.5, 0.7, 0.9]))
    print "Mesh center coords",okIfTrue(centerOK)
    if not centerOK:
        print "Expected: [0.1, 0.3, 0.5, 0.7, 0.9]"
        print "Computed: ",xc
    faceOK = isSmallDiff(xf, np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]))
    print "Mesh boundary coordinates",okIfTrue(faceOK)
    if not faceOK:
        print "Expected: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]"
        print "Computed: ",xf


def test_lds():
    mesh = meshmaker.mkmesh(5, 0.0, 1.0)
    F, D, srcCoeffs, bdrVals = init(mesh)
    stdFormCoefss = lds.mkcoeffs(mesh, F, D, srcCoeffs, bdrVals)
    aP, aW, aE, b = stdFormCoeffs
    aPOK = isSmallDiff(aP, np.array([ 2.75,  1.0,   1.00,  1.00,  0.25]))
    aWOK = isSmallDiff(aW, np.array([ 0.00,  1.75,  1.75,  1.75,  1.75]))
    aEOK = isSmallDiff(aE, np.array([-0.75, -0.75, -0.75, -0.75,  0.00]))
    bOK  = isSmallDiff(b,  np.array([ 3.50,  0.00,  0.00,  0.00,  0.00]))
    print "aP", okIfTrue(aPOK)
    if not aPOK:
        print "Expected: [ 2.75,  1.0,   1.00,  1.00,  0.25]"
        print "Computed: ",aP
    print "aW", okIfTrue(aWOK)
    if not aWOK:
        print "Expected: [ 0.00,  1.75,  1.75,  1.75,  1.75]"
        print "Computed: ",aW
    print "aE", okIfTrue(aEOK)
    if not aEOK:
        print "Expected: [-0.75, -0.75, -0.75, -0.75,  0.00]"
        print "Computed: ",aE
    print "b ", okIfTrue(bOK)
    if not bOK:
        print "Expected: [ 3.50,  0.00,  0.00,  0.00,  0.00]"
        print "Computed: ",b


def test_uw():
    pass

def test_hyb():
    pass

def test_ais():
    pass

