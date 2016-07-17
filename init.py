#coding: utf8
#File init.py
#José Amoreira
#July 2016

#Function init(xf, xc) sets the initial values of all physical parameters for
#the 1D, steady-state, diffusion-convection (fixed convection) problem

import numpy as np



def init(xc, xf):
    """init(). Defines the values for density, velocity, conduction coefficient,
    convection and diffusion strength parameters F and D
       source terms and boundary values at each CV.
    """
    #Values below taken from Versteeg & Malalasekera Ex. 5.1b.
    #UPDATE COMMENT IF YOU CHANGE THE VALUES!

    #Density, velocity, convection strength parameter F:
    rho = 1.0
    v = 2.5
    F = rho * v

    #Conduction coefficient and diffusion strength parameter D:
    gamma = 0.1 
    D = np.zeros_like(xf)
    D [1:-1] = gamma / (xc[1:] - xc[:-1])
    D[0] = gamma / (xc[0] - xf[0])
    D[-1] = gamma / (xf[-1] - xc[-1])

    #Boundary values
    phib_W = 1.0
    phib_E = 0.0

    return F, D, phib_W, phib_E
