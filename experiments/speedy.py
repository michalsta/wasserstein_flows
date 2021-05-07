import networkx
import IsoSpecPy
import math
import sys
from pprint import pprint



def empiric_gradient_iso(exp, the_l, point, exp_ab_cost, th_ab_cost, dval = 0.0001):
    res = []
    the_l = [i.copy() for i in the_l]
    for i, sf in zip(the_l, point):
        i.scale(sf)

    the_s = IsoSpecPy.IsoDistribution.LinearCombination(the_l, point)
    base = exp.abyssalWassersteinDistance(the_s)
    print("base:", base)

    for idx in range(len(the_l)):
        p = point[:]
        p[idx] += dval
        the_s = IsoSpecPy.IsoDistribution.LinearCombination(the_l, p)
        grad_part = (exp.abyssalWassersteinDistance(the_s, exp_ab_cost + the_ab_cost) - base) / dval
        res.append(grad_part)

    return res

import random
import numpy as np
N = 10
def get_probs():
    res = [random.randint(1,10) for x in range(N)]
    #res = probs + [N*5 - sum(probs)]
    #print("Probs:", res)
    return res
def get_masses():
    res = list(np.random.choice(range(N*3), size=N, replace=False))
#        res = [3,3]
    while len(res) != len(set(res)):
        res = [random.randint(0,3+N) for x in range(N)]
    #print("Masses:", res)
    return res

def random_spectrum():
    return IsoSpecPy.IsoDistribution(masses=get_masses(), probs = get_probs())

import scipy.optimize

def rescaled(iso, scale):
    ret = iso.copy()
    ret.scale(scale)
    return ret

def deconvolve(EXP, THEs, exp_ab_cost, the_ab_cost):
    def minimized_f(point):
        if any(x < 0.0 for x in point):
            return 10000000000000000.0
        scsp = [rescaled(THE, x_i) for THE, x_i in zip(THEs, point)]
        res = awsd(EXP, scsp, exp_ab_cost, the_ab_cost)
        print(point, res)
        return res
    def grad(point):
        return empiric_gradient(EXP, THEs, point, exp_ab_cost, the_ab_cost)
    print("susp:", minimized_f([1.0, 0.8]))

    print(scipy.optimize.minimize(minimized_f, [10.1] * len(THEs), jac=grad))

if __name__ == '__main__':
    EXP = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [1.0, 1.0])
    THE1 = IsoSpecPy.IsoDistribution(masses=[1.001, 2.0], probs = [0.2, 0.1])

    THE2 = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [1.0, 0.0])
    THEs = [THE1, THE2]

    print(awsd(EXP, THEs, 0.5, 0.5))
    print(empiric_gradient(EXP, THEs, [1.0, 1.0], 0.5, 0.5))
    sys.exit(0)
    deconvolve(EXP, THEs, 0.5, 0.5)

#    EXP = random_spectrum()
#    THEs = [random_spectrum() for _ in range(10)]

    print(empiric_gradient(EXP, THEs, [0.1]*len(THEs), 0.5, 0.5))
    lkklnjk

    i = 1
    n = 1
    while i == n:
        EXP = IsoSpecPy.IsoDistribution(masses=get_masses(), probs = get_probs())
        THE1 = IsoSpecPy.IsoDistribution(masses=get_masses(), probs = get_probs())
        #EXP = IsoSpecPy.IsoDistribution(masses=[1.0, 2.1], probs = [1.0, 0.3])
        #THE1 = IsoSpecPy.IsoDistribution(masses=[2.0, 3.0], probs = [1.1, 2.0])

        n = awsd(EXP, [THE1], 0.5, 0.5)
        i = EXP.abyssalWassersteinDistance(THE1, 1.0, 1.0)
        print(i, n)
#    print(awsd(EXP, [THE1], 0.5, 0.5))
#    print(EXP.abyssalWassersteinDistance(THE1, 1.0, 1.0))
