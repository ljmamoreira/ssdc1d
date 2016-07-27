#coding: utf8
#File analytic.py
#José Amoreira
#July 2016

#Calculation of analytic solution and L2 error
#The analytic function is created with the desired parameters via
#myFunc = AnalyticSolution(x0,x1, (y0,y1), mu)
#Then it can be called simply by: myFunc(x)
#To find the error of a solution x,y: myFunc.error(x,y)

import numpy as np

def _solution(x, x0, x1, y0, y1, mu):
    ex = np.exp(mu*x)
    e0 = np.exp(mu*x0); e1 = np.exp(mu*x1)
    y = y0 + (y1-y0)/(e1-e0)*(ex-e0)
    return y

class AnalyticSolution(object):
    def __init__(self, x0, x1, bdrVals, mu):
        self.x0 = x0
        self.x1 = x1
        self.phiA, self.phiB = bdrVals
        self.mu = mu

    def __call__(self, x):
        return _solution(x, self.x0, self.x1, self.phiA, self.phiB, self.mu)

    def error(self,x,yc):
        delta = yc - self(x)
        return np.sqrt(np.dot(delta,delta)/len(delta))


