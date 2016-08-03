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


#pretty printer for arrays
def printArray(a,label=""):
    print label+" "+' '.join('{:7.4f}'.format(x) for x in a)

#Define here a local version of init(), in order to keep the parameters
def init(mesh,velocity):
    """init(mesh). Defines the values for density, velocity, conduction
       coefficient, convection and diffusion strength parameters F and D
       source terms and boundary values at each CV.
    """
    xc, xf = mesh

    #Values below taken from Versteeg & Malalasekera Ex. 5.1b.
    #DO NOT CHANGE THE VALUES! 

    #Density, velocity, convection strength parameter F:
    rho = 1.0
    v = velocity
    F = np.ones_like(xf) * rho * v

    #Conduction coefficient and diffusion strength parameter D:
    gamma = 0.1 
    D = np.zeros_like(xf)
    D [1:-1] = gamma / (xc[1:] - xc[:-1])
    D[0] = gamma / (xc[0] - xf[0])
    D[-1] = gamma / (xf[-1] - xc[-1])

    printArray(F/D, "P:")
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

    return (F, D), srcCoeffs, bdrVals, mu

#Small auxiliary functions for reporting results conformance
def isSmallDiff(a1,a2):
    d = a1 - a2
    return np.dot(d,d) < np.finfo(float).eps

def okIfTrue(ok):
    if ok:
        return "OK"
    return "Wrong"


def report_aPWEb(computed, versteeg):
    aPc, aWc, aEc, bPc = computed
    aPv, aWv, aEv, bPv = versteeg
    print ("{:2s}  "+4*"{:^7s} {:^7s}  ").format('i', 'aPc', 'aPv', 'aWc',
            'aWv', 'aEc', 'aEv', 'bPc', 'bPv')
    for (i, (apc, apv, awc, awv, aec, aev, bc, bv)) in enumerate(zip(
        aPc, aPv, aWc, aWv, aEc, aEv, bPc, bPv)):
        print ('{: 2d}  '+4*'{:7.4f} {:7.4f}  ').format(
            i+1, apc, apv, awc, awv, aec, aev, bc, bv)
    sumUp ='    ' + '  '.join('{:^15s}'.format(okIfTrue(isSmallDiff(ac,av)))
        for ac, av in zip(computed, versteeg))
    print sumUp
    return reduce(lambda u,v: u and v, [isSmallDiff(ac,av) for ac, av in
        zip(computed, versteeg)])


def report_y(xc, yc, yv, ya):
    print
    print ("{:^2s}  {:^3s}" + 3*"  {:^6s}").format('i', 'x', 'yc', 'yv', 'ya')
    for i, (x, fc, fv, fa) in enumerate(zip(xc, yc, yv, ya)):
        print ("{:2d}  {:3.1f}"+3*"  {:6.4f}").format(i,x,fc,fv,fa)
    result = isSmallDiff(yc,yv)
    print "Solution: ", okIfTrue(result)
    return result
    



#Tests start here

def test_ucmesh():
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    xc, xf = mesh
    xcExpected = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    xfExpected = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    printReport("Mesh center coords", xcExpected, xc)
    printReport("Mesh face coords", xfExpected, xf)



#5.1(i): v=0.1, N=5, LDS
def test_V51i():
    print "\n"+30*"="
    print "V&M, Example 5.1(i): LDS, N=5, v=0.1"
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    FD, srcCoeffs, bdrVals, mu = init(mesh, velocity=0.1)
    stdFormCoeffs = discr.stdEqCoeffs("lds",mesh, FD, srcCoeffs, bdrVals)
    aP, aW, aE, bP = stdFormCoeffs
    aP_v = [1.55, 1.00, 1.00, 1.00, 1.45]
    aW_v = [0.00, 0.55, 0.55, 0.55, 0.55]
    aE_v = [0.45, 0.45, 0.45, 0.45, 0.00]
    bP_v = [1.10, 0.00, 0.00, 0.00, 0.00]
    passed = report_aPWEb((aP,aW,aE,bP), (aP_v,aW_v,aE_v,bP_v))
    yc = np.around(solve.solve(stdFormCoeffs),decimals=4)
    xc, xf = mesh
    fa = analytic.AnalyticSolution(0.0, 1.0, bdrVals, mu)
    ya = fa(xc)
    yv = [0.9421, 0.8006, 0.6276, 0.4163, 0.1579]
    passed = passed and report_y(xc, yc, yv, ya)
    return passed



#5.1(ii): v=2.5, N=5, LDS
def test_V51ii():
    print "\n"+30*"="
    print "V&M, Example 5.1(ii): LDS, N=5, v=2.5"
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    FD, srcCoeffs, bdrVals, mu = init(mesh, velocity=2.5)
    stdFormCoeffs = discr.stdEqCoeffs("lds",mesh, FD, srcCoeffs, bdrVals)
    aP, aW, aE, bP = stdFormCoeffs
    aP_v = [2.75, 1.00, 1.00, 1.00, 0.25]
    aW_v = [0.00, 1.75, 1.75, 1.75, 1.75]
    aE_v = [-0.75, -0.75, -0.75, -0.75, 0.00]
    bP_v = [3.50, 0.00, 0.00, 0.00, 0.00]
    passed = report_aPWEb((aP,aW,aE,bP), (aP_v,aW_v,aE_v,bP_v))
    yc = np.around(solve.solve(stdFormCoeffs),decimals=4)
    xc, xf = mesh
    fa = analytic.AnalyticSolution(0.0, 1.0, bdrVals, mu)
    ya = fa(xc)
    yv = [1.0356, 0.8694, 1.2573, 0.3521, 2.4644]
    passed = report_y(xc, yc, yv, ya) and passed
    return passed


#5.2(i): v=0.1, N=5, UWS
def test_V52i():
    print "\n"+30*"="
    print "V&M, Example 5.2(i): UWS, N=5, v=0.1"
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    FD, srcCoeffs, bdrVals, mu = init(mesh, velocity=0.1)
    stdFormCoeffs = discr.stdEqCoeffs("uws", mesh, FD, srcCoeffs, bdrVals)
    aP, aW, aE, bP = stdFormCoeffs
    aP_v = np.array([ 1.60,  1.10,  1.10,  1.10,  1.60])
    aW_v = np.array([ 0.00,  0.60,  0.60,  0.60,  0.60])
    aE_v = np.array([ 0.50,  0.50,  0.50,  0.50,  0.00])
    bP_v = np.array([ 1.10,  0.00,  0.00,  0.00,  0.00])
    passed = report_aPWEb((aP,aW,aE,bP), (aP_v,aW_v,aE_v,bP_v))
    yc = np.around(solve.solve(stdFormCoeffs),decimals=4)
    xc, xf = mesh
    fa = analytic.AnalyticSolution(0.0, 1.0, bdrVals, mu)
    ya = fa(xc)
    yv = [0.9337, 0.7879, 0.6130, 0.4031, 0.1512]
    passed = report_y(xc, yc, yv, ya) and passed
    return passed



#5.2(ii): v=2.5, N=5, UWS
def test_V52i():
    print "\n"+30*"="
    print "V&M, Example 5.2(ii): UWS, N=5, v=2.5"
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    FD, srcCoeffs, bdrVals, mu = init(mesh, velocity=2.5)
    stdFormCoeffs = discr.stdEqCoeffs("uws", mesh, FD, srcCoeffs, bdrVals)
    aP, aW, aE, bP = stdFormCoeffs
    aP_v = np.array([ 4.00,  3.50,  3.50,  3.50,  4.00])
    aW_v = np.array([ 0.00,  3.00,  3.00,  3.00,  3.00])
    aE_v = np.array([ 0.50,  0.50,  0.50,  0.50,  0.00])
    bP_v = np.array([ 3.50,  0.00,  0.00,  0.00,  0.00])
    passed = report_aPWEb((aP,aW,aE,bP), (aP_v,aW_v,aE_v,bP_v))
    yc = np.around(solve.solve(stdFormCoeffs),decimals=4)
    xc, xf = mesh
    fa = analytic.AnalyticSolution(0.0, 1.0, bdrVals, mu)
    ya = fa(xc)
    yv = [0.9998, 0.9987, 0.9921, 0.9524, 0.7143]
    passed = report_y(xc, yc, yv, ya) and passed
    return passed



#5.3: v=2.5, N=5, HYB
def test_V53():
    print "\n"+30*"="
    print "V&M, Example 5.3: HYB, N=5, v=2.5"
    mesh = meshmaker.ucmesh(5, 0.0, 1.0)
    FD, srcCoeffs, bdrVals, mu = init(mesh, velocity=2.5)
    stdFormCoeffs = discr.stdEqCoeffs("hyb", mesh, FD, srcCoeffs, bdrVals)
    aP, aW, aE, bP = stdFormCoeffs
    aP_v = np.array([ 3.50,  2.50,  2.50,  2.50,  3.50])
    aW_v = np.array([ 0.00,  2.50,  2.50,  2.50,  2.50])
    aE_v = np.array([ 0.00,  0.00,  0.00,  0.00,  0.00])
    bP_v = np.array([ 3.50,  0.00,  0.00,  0.00,  0.00])
    passed = report_aPWEb((aP,aW,aE,bP), (aP_v,aW_v,aE_v,bP_v))
    yc = np.around(solve.solve(stdFormCoeffs),decimals=4)
    xc, xf = mesh
    fa = analytic.AnalyticSolution(0.0, 1.0, bdrVals, mu)
    ya = fa(xc)
    yv = [1.0000, 1.0000, 1.0000, 1.0000, 0.7143]
    passed = report_y(xc, yc, yv, ya) and passed
    return passed




def test_cai():
    mesh = meshmaker.ucmesh(25, 0.0, 1.0)
    FD, srcCoeffs, bdrVals, mu = init(mesh)
    stdFormCoeffs = discr.stdEqCoeffs("cai", mesh, FD, srcCoeffs, bdrVals)
    aP, aW, aE, b = stdFormCoeffs

    printArray(aW, "aW:")
    printArray(aE, "aE:")
    printArray(aP, "aP:")
    printArray(b,  "bP:")
    printArray((np.abs(aE) + np.abs(aW))/np.abs(aP), "Sc:")
    ycai = np.around(solve.solve(stdFormCoeffs),decimals=4)
    xc, xf = mesh
    f_ana = analytic.AnalyticSolution(0.0, 1.0, bdrVals, mu)
    yanalytic = f_ana(xc)
    print "  #      x        yc        ya"
    for i,(x,y,ya) in enumerate(zip(xc,ycai,yanalytic)):
        print ("{:3d}  "+"{:4f}  "*3).format(i+1, x, y, ya)
    print "Error: ",f_ana.error(xc,ycai)




if __name__ == "__main__":
    passed = test_V51i()
    print "Global result:", okIfTrue(passed)
    passed = test_V51ii() and passed
    print "Global result:", okIfTrue(passed)
    passed = test_V52i() and passed
    print "Global result:", okIfTrue(passed)
    passed = test_V53() and passed
    print "Global result:", okIfTrue(passed)
    
