#coding: utf8
#File tests.py
#José Amoreira
#July 2016

#Calculation of analytic solution

import numpy as np

def solution(x, x0, x1, bdrVals, mu):
    y0,y1 = bdrVals
    ex = np.exp(mu*x)
    e0 = np.exp(mu*x0); e1 = np.exp(mu*x1)
    y = y0 + (y1-y0)/(e1-e0)*(ex-e0)
    return y


