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
from scipy.optimize import minimize




#def fvL2Error(meth, mesh, FD, srcCoeffs, bdrVals, xsi):
#    stdFormCoeffs = discr.mkCoeffs(meth, mesh, FD, srcCoeffs, bdrVals, xsi=xsi)
#    yfv = solve.solve(stdFormCoeffs)
#    xc = mesh[0]
#    v, rho, gamma = physPars
#    mu = rho * v / gamma
#    f_ana = aux.AnalyticSolution(0.0, 1.0, bdrVals, mu)
#    return f_ana.error(xc,yfv)

class FVSolution(object):
    def __init__(self,mesh,meth,FD,srcCoeffs, bdrVals, xsi):
        self.mesh = mesh
        self.meth = meth
        self.FD = FD
        self.srcCoeffs = srcCoeffs
        self.bdrVals = bdrVals
        self.xsi = xsi
        F,D = FD
        self.xc, self.xf = mesh
        dx0 = xc[0] - xf[0]
        mu = F[0] / (D[0] * dx0)
        self.f_ana = aux.AnalyticSolution(0.0, 1.0, bdrVals, mu)

    def compute(self):
        self.stdFormCoeffs = discr.mkCoeffs(self.meth, self.mesh, self.FD,
                self.srcCoeffs, self.bdrVals, xsi=self.xsi)
        self.yfv = solve.solve(self.stdFormCoeffs)


    def eofxsi(self,xsi):
        self.xsi = xsi
        self.compute()
        return self.f_ana.error(self.xc,self.yfv)
        


if __name__ == "__main__":

    #Default parameter values
    n = 5
    v = 2.5; rho=1.0; gamma = 0.1
    phib_W=1.0; phib_E=0.0
    xsi = 0.5
    meth = "cai"
    #Update from assignments (syntax: var=val) in the command line 
    for assignment in sys.argv[1:]:
        var,val = (x.strip() for x in assignment.split("="))
        locals()[var] = type(locals()[var])(val)

    bdrVals = (phib_W, phib_E)

    M = np.zeros(n)
    N = np.zeros(n)
    srcCoeffs = (M,N)


    mesh = meshmaker.ucmesh(n, 0.0, 1.0)
    xc, xf = mesh

    bxsi = []
    print ("{:^7s}"+4*"{:^12s}"+"{:^6s}").format("v","lds","uws","hyb","cai", "bxsi")
    for v in (0.1*i for i in range(1,41)):
        line = "{:4.2f}   ".format(v)
        physPars = (v, rho, gamma)
        F = np.ones_like(xf) * rho * v
        D = np.zeros_like(xf)
        D[1:-1] = gamma / (xc[1:] - xc[:-1])
        D[0] = gamma / (xc[0] - xf[0])
        D[-1] = gamma / (xf[-1] - xc[-1])
        FD = (F,D)
        for meth in ("lds", "uws", "hyb", "cai"):
            FVsys = FVSolution(mesh,meth,FD,srcCoeffs, bdrVals, xsi)

            if meth == "cai":
                optimal = minimize(FVsys.eofxsi, 0.5)
                err = optimal.fun
                best_xsi = optimal.x[0]
                bxsi.append(best_xsi)
                line += "{:6.4e}  {:6.4f}".format(err,best_xsi)
            else:
                err = FVsys.eofxsi(0.5)
                line += "{:6.4e}  ".format(err)
        print line

    print
    bxsi = np.array(bxsi)
    best_xsi = np.average(bxsi)
    print "Average xsi:", best_xsi
    print "Std xsi:    ", np.std(bxsi)

    with open('errtab.dat','w') as ofile:
        print ("{:^7s}"+4*"{:^12s}").format("v","lds","uws","hyb","cai")
        for v in (0.1*i for i in range(1,41)):
            line = "{:4.2f}   ".format(v)
            physPars = (v, rho, gamma)
            F = np.ones_like(xf) * rho * v
            D = np.zeros_like(xf)
            D[1:-1] = gamma / (xc[1:] - xc[:-1])
            D[0] = gamma / (xc[0] - xf[0])
            D[-1] = gamma / (xf[-1] - xc[-1])
            FD = (F,D)
            for meth in ("lds", "uws", "hyb", "cai"):
                FVsys = FVSolution(mesh,meth,FD,srcCoeffs, bdrVals, xsi)
                err = FVsys.eofxsi(best_xsi)
                line += "{:6.4e}  ".format(err)
            ofile.write(line+'\n')
            print line


