#coding: utf8
#File ss1dcd.py
#José Amoreira
#August 2016

#Compute solutions to steady-state 1D difusion with given convection



import numpy as np
import sys
import meshmaker
import discretize as discr
import solve
import aux


if __name__ == "__main__":

    #Default parameter values
    n = 5
    v = 2.5; rho=1.0; gamma = 0.1;
    phib_W=1.0; phib_E=0.0
    xsi = 0.5
    meth = "cai"
    #Update from assignments (syntax: var=val) in the command line 
    for assignment in sys.argv[1:]:
        var,val = (x.strip() for x in assignment.split("="))
        locals()[var] = type(locals()[var])(val)

    print "Parameters:"
    print "n =  ", n
    print "meth=", meth
    print "v =  ", v
    print "rho =", rho
    print "gam =", gamma
    print "phia=", phib_W
    print "phib=", phib_E
    print "xsi= ", xsi

    mu = rho * v / gamma

    mesh = meshmaker.ucmesh(n, 0.0, 1.0)
    xc, xf = mesh

    F = np.ones_like(xf) * rho * v
    D = np.zeros_like(xf)
    D[1:-1] = gamma / (xc[1:] - xc[:-1])
    D[0] = gamma / (xc[0] - xf[0])
    D[-1] = gamma / (xf[-1] - xc[-1])
    FD = (F,D)

    M = np.zeros_like(xc)
    N = np.zeros_like(xc)
    srcCoeffs = (M,N)

    bdrVals = (phib_W, phib_E)

    stdFormCoeffs = discr.mkCoeffs(meth, mesh, FD, srcCoeffs, bdrVals, xsi=xsi)
    aP, aW, aE, b = stdFormCoeffs
    
    aux.printArray(aW, "aW:")
    aux.printArray(aE, "aE:")
    aux.printArray(aP, "aP:")
    aux.printArray(b,  "bP:")
    aux.printArray((np.abs(aE) + np.abs(aW))/np.abs(aP), "Sc:")
    ycai = np.around(solve.solve(stdFormCoeffs),decimals=4)
    xc, xf = mesh
    f_ana = aux.AnalyticSolution(0.0, 1.0, bdrVals, mu)
    yanalytic = f_ana(xc)
    print "  #      x        yc        ya"
    for i,(x,y,ya) in enumerate(zip(xc,ycai,yanalytic)):
        print ("{:3d}  "+"{:4f}  "*3).format(i+1, x, y, ya)
    print "Error: ",f_ana.error(xc,ycai)


