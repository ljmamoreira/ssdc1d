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




#GFVSolution: General finite volume solution for LDS, UWS or HYB schemes
class GFVSolution(object):
    def __init__(self,mesh,scheme,FD,srcCoeffs, bdrVals):
        self.mesh = mesh
        self.scheme = scheme
        self.FD = FD
        self.srcCoeffs = srcCoeffs
        self.bdrVals = bdrVals
        F,D = FD
        xc, xf = mesh
        dx0 = xc[0] - xf[0]
        mu = F[0] / (D[0] * dx0)
        self.f_ana = aux.AnalyticSolution(0.0, 1.0, bdrVals, mu)
        #self.xsi is needed for conformance with the interface of
        #discretize.mkCoeffs (in compute())
        self.xsi = 0

    def compute(self):
        self.stdFormCoeffs = discr.mkCoeffs(self.scheme, self.mesh, self.FD,
                self.srcCoeffs, self.bdrVals, xsi=self.xsi)
        self.yfv = solve.solve(self.stdFormCoeffs)

    def error(self):
        self.compute()
        xc = self.mesh[0]
        return self.f_ana.error(xc, self.yfv)

#CASSolution: CAS solution, iherits from GVFSolution
class CASSolution(GFVSolution):
    def __init__(self,mesh,scheme,FD,srcCoeffs,bdrVals,xsi):
        super(CASSolution, self).__init__(mesh, scheme, FD, srcCoeffs, bdrVals)
        self.xsi = xsi

    def eofxsi(self,xsi):
        save_xsi = self.xsi
        self.xsi = xsi
        result = self.error()
        self.xsi = save_xsi
        return result


#Common interface to create FVSolution objects
def fvsolution(mesh,scheme,FD,srcCoeffs,bdrVals,*args):
    if scheme == "cas":
        xsi = args[0]
        return CASSolution(mesh,scheme,FD,srcCoeffs,bdrVals,xsi)
    else:
        return GFVSolution(mesh,scheme,FD,srcCoeffs, bdrVals)



def best_xsi(mesh, physPars, bdrVals, srcCoeffs):
    v, rho, gamma = physPars
    best_xsis = []
    for v in np.linspace(0.1, 5.0, 50):
        physPars = v, rho, gamma
        F = np.ones_like(xf) * rho * v
        D = np.zeros_like(xf)
        D[1:-1] = gamma / (xc[1:] - xc[:-1])
        D[0] = gamma / (xc[0] - xf[0])
        D[-1] = gamma / (xf[-1] - xc[-1])
        FD = F,D
        FVSys = fvsolution(mesh, "cas", FD, srcCoeffs, bdrVals, 0.5)
        optimal = minimize(FVSys.eofxsi, 0.5)
        best_xsis.append(optimal.x[0])
        #print "{:4.1f}  {:7.4f}  {:11.4e}".format(v,optimal.x[0],optimal.fun)
    bxsi = np.average(best_xsis)
    return bxsi



def convergence_rate(physPars, bdrVals, xsi):
    v, rho, gamma = physPars
    with open("conv_rate.dat", "w") as erfile:
        for n in [5, 11, 21, 41, 81, 161, 321]:
            dx = 1.0 / n
            mesh = meshmaker.ucmesh(n, 0.0, 1.0)
            xc, xf = mesh
            F = np.ones_like(xf) * rho * v
            D = np.zeros_like(xf)
            D[1:-1] = gamma / (xc[1:] - xc[:-1])
            D[0] = gamma / (xc[0] - xf[0])
            D[-1] = gamma / (xf[-1] - xc[-1])
            FD = (F,D)
            M = np.zeros(n)
            N = np.zeros(n)
            srcCoeffs = (M,N)
            line = "{:6.4f}  ".format(dx)
            for scheme in ("lds", "uws", "hyb", "cas"):
                FVsys = fvsolution(mesh,scheme,FD,srcCoeffs, bdrVals, xsi)
                err = FVsys.error()
                line += "{:6.4e}  ".format(err)
            print line
            erfile.write(line+'\n')


def standardSol(n,scheme,xsi=0.5):
    mesh = meshmaker.ucmesh(n, 0.0, 1.0)
    v = 2.5; rho = 1.0; gamma = 0.1; 
    xc,xf = mesh
    F = np.ones(n+1) * v * rho
    D = np.zeros(n+1); 
    D[0] = gamma/(xc[0]-xf[0])
    D[1:-1] = gamma/(xc[1:]-xc[:-1])
    D[-1] = gamma/(xf[-1] - xc[-1])
    FD = (F,D)
    srcCoeffs = (np.zeros(n), np.zeros(n))
    bdrVals = (1.0, 0.0)
    return fvsolution(mesh, scheme, FD, srcCoeffs, bdrVals, xsi)



def mkFD(physPars, mesh):
    v, rho, gamma = physPars
    xc, xf = mesh
    F = np.ones_like(xf) * rho * v
    D = np.zeros_like(xf)
    D[1:-1] = gamma / (xc[1:]- xc[:-1])
    D[0] = gamma / (xc[0] - xf[0])
    D[-1] = gamma / (xf[-1] - xc[-1])
    return (F,D)


if __name__ == "__main__":

    #Default parameter values
    n = 5
    v = 2.5; rho=1.0; gamma = 0.1
    phib_W=1.0; phib_E=0.0
    xsi = 0.5
    scheme = "cas"
    #Update from assignments (syntax: var=val) in the command line 
    for assignment in sys.argv[1:]:
        var,val = (x.strip() for x in assignment.split("="))
        locals()[var] = type(locals()[var])(val)

    physPars = (v, rho, gamma)
    bdrVals = (phib_W, phib_E)

    M = np.zeros(n)
    N = np.zeros(n)
    srcCoeffs = (M,N)

    mesh = meshmaker.ucmesh(n, 0.0, 1.0)
    xc, xf = mesh
    
    #Determination of best xsi in the velocity range [0:5]
    bxsi = best_xsi(mesh, physPars, bdrVals, srcCoeffs)
    print "best_xsi:", bxsi

    #Compare solutions
    physPars = (v, rho, gamma)
    FD = mkFD(physPars, mesh)
    lds = fvsolution(mesh, "lds", FD, srcCoeffs, bdrVals);
    uws = fvsolution(mesh, "uws", FD, srcCoeffs, bdrVals)
    hyb = fvsolution(mesh, "hyb", FD, srcCoeffs, bdrVals)
    cas = fvsolution(mesh, "cas", FD, srcCoeffs, bdrVals, bxsi)
    lds.compute()
    uws.compute()
    hyb.compute()
    cas.compute()
    print cas.f_ana(mesh[0])
    with open('phi.dat','w') as pfile:
        print (5*"{:^8s}").format("x", "lds", "uws", "hyb", "cas")
        for x,y1,y2,y3,y4 in zip(mesh[0], lds.yfv, uws.yfv, hyb.yfv, cas.yfv):
            line = (5*"{:8.4f}").format(x,y1,y2,y3,y4)
            print line
            pfile.write(line+'\n')



    #Compare errors
    with open('errtab.dat','w') as ofile:
        print ("{:^7s}"+4*"{:^12s}").format("v","lds","uws","hyb","cas")
        for v in (0.1*i for i in range(1,41)):
            line = "{:4.2f}   ".format(v)
            physPars = (v, rho, gamma)
            F = np.ones_like(xf) * rho * v
            D = np.zeros_like(xf)
            D[1:-1] = gamma / (xc[1:] - xc[:-1])
            D[0] = gamma / (xc[0] - xf[0])
            D[-1] = gamma / (xf[-1] - xc[-1])
            FD = (F,D)
            for scheme in ("lds", "uws", "hyb", "cas"):
                FVsys = fvsolution(mesh,scheme,FD,srcCoeffs, bdrVals, bxsi)
                err = FVsys.error()
                line += "{:6.4e}  ".format(err)
            ofile.write(line+'\n')
            print line

    v = 2.5; rho = 1.0; gamma = 0.1; physPars = (v, rho, gamma)
    convergence_rate(physPars, bdrVals, bxsi)
