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

def test_mkmesh():
    mesh = meshmaker.mkmesh(5, 0.0, 1.0)
    xc, xf = mesh
    dc = xc - np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    if np.dot(dc,dc)<1.e-20:
        print "Mesh center coordinates OK."
    else:
        print "Mesh center coordinates wrong."
        print "Expected: [0.1, 0.3, 0.5, 0.7, 0.9]"
        print "Computed: ",xc
    df = xf - np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    if np.dot(df,df)<1.e-20:
        print "Mesh boundary coordinates OK."
    else:
        print "Mesh boundary coordinates wrong."
        print "Expected: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]"
        print "Computed: ",xf
    

def test_uw():
    pass

def test_hyb():
    pass

def test_ais():
    pass

