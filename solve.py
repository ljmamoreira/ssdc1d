#coding: utf8
#File solve.py
#José Amoreira
#July 2016

#Form system of equations from standard form coefficiens; compute and return the
#solution

import numpy as np

def solve(stdFormCoeffs):
    aP, aW, aE, b = stdFormCoeffs
    a = np.diag(aP) - np.diag(aE[:-1], k=1) - np.diag(aW[1:], k=-1)
    y = np.linalg.solve(a,b)
    return y
