#coding: utf8
#File tests.py
#José Amoreira
#July 2016

#Solve Example 5.1 and following from Versteeg and Malalakesera, check the
#results from those presented there


import numpy as np
import meshmaker
import cddiscretize as discr
#import lds
#import uw
import solve
import analytic

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

    #mu: parameter in the analytic solution
    mu = rho * v / gamma

    return F, D, srcCoeffs, bdrVals, mu

#Small auxiliary functions for reporting results conformance
def isSmallDiff(a1,a2):
    d = a1 - a2
    return np.dot(d,d) < np.finfo(float).eps

def okIfTrue(ok):
    if ok:
        return "OK."
    return "wrong."

def printReport(name,expected,computed):
    isOK = isSmallDiff(expected, computed)
    print name,okIfTrue(isOK)
    if not isOK:
        print "   Expected: ", expected
        print "   Computed: ", computed


#Tests start here

def test_ucmesh():
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    xc, xf = mesh
    xcExpected = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    xfExpected = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    printReport("Mesh center coords", xcExpected, xc)
    printReport("Mesh face coords", xfExpected, xf)


def test_lds():
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    F, D, srcCoeffs, bdrVals, mu = init(mesh)
    stdFormCoeffs = discr.stdEqCoeffs("lds",mesh, (F, D), srcCoeffs, bdrVals)
    aP, aW, aE, b = stdFormCoeffs
    aPExpected = np.array([ 2.75,  1.0,   1.00,  1.00,  0.25])
    aWExpected = np.array([ 0.00,  1.75,  1.75,  1.75,  1.75])
    aEExpected = np.array([-0.75, -0.75, -0.75, -0.75,  0.00])
    bExpected =  np.array([ 3.50,  0.00,  0.00,  0.00,  0.00])
    printReport("aP", aPExpected, aP)
    printReport("aW", aWExpected, aW)
    printReport("aE", aEExpected, aE)
    printReport("b",  bExpected, b)
    ylds = np.around(solve.solve(stdFormCoeffs),decimals=4)
    yexpected = [1.0356, 0.8694, 1.2573, 0.3521, 2.4644]
    xc, xf = mesh
    yanalytic = analytic.solution(xc, 0.0, 1.0, bdrVals, mu)
    for i,(x,y,ye,ya) in enumerate(zip(xc,ylds,yexpected,yanalytic)):
        print ("{:2d} "+"{:4f} "*4).format(i+1, x, y,ye, ya)
    


def test_uw():
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    F, D, srcCoeffs, bdrVals, mu = init(mesh)
    stdFormCoeffs = uw.mkcoeffs(mesh, F, D, srcCoeffs, bdrVals)
    aP, aW, aE, b = stdFormCoeffs
    aPExpected = np.array([ 4.00,  3.50,  3.50,  3.50,  4.00])
    aWExpected = np.array([ 0.00,  3.00,  3.00,  3.00,  3.00])
    aEExpected = np.array([ 0.50,  0.50,  0.50,  0.50,  0.00])
    bExpected =  np.array([ 3.50,  0.00,  0.00,  0.00,  0.00])
    printReport("aP", aPExpected, aP)
    printReport("aW", aWExpected, aW)
    printReport("aE", aEExpected, aE)
    printReport("b",  bExpected, b)
    yuw = np.around(solve.solve(stdFormCoeffs),decimals=4)
    yexpected = [0.9998, 0.9987, 0.9921, 0.9524, 0.7143]
    xc, xf = mesh
    yanalytic = analytic.solution(xc, 0.0, 1.0, bdrVals, mu)
    for i,(x,y,ye,ya) in enumerate(zip(xc,yuw,yexpected,yanalytic)):
        print ("{:2d} "+"{:4f} "*4).format(i+1, x, y,ye, ya)

def test_hyb():
    pass

def test_ais():
    pass

