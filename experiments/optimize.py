import networkx
import IsoSpecPy
import math
import sys
from pprint import pprint

from parameters import int_fact, integerize
from distances import awsd_checked as awsd
from gradients import checked_grad

def empiric_gradient(exp, the_l, point, exp_ab_cost, th_ab_cost):
    dval = 0.0001
    res = []
    assert dval * int_fact > 1.0
    the_l = [i.copy() for i in the_l]
    for i, sf in zip(the_l, point):
        i.scale(sf)

    base = awsd(exp, the_l, exp_ab_cost, th_ab_cost)
    print("base:", base)

    for idx in range(len(the_l)):
        c = the_l[idx].copy()
        c.scale(dval + 1.0) #point[idx])
        n_the_l = the_l[:idx] + [c] + the_l[idx+1:]
        print("Aaa")
        grpart = (awsd(exp, n_the_l, exp_ab_cost, th_ab_cost) - base) / dval
        print("grpart", awsd(exp, n_the_l, exp_ab_cost, th_ab_cost))
        res.append(grpart)

    return res

def empiric_gradient_iso(exp, the_l, point, exp_ab_cost, th_ab_cost):
    dval = 0.0001
    res = []
    the_l = [i.copy() for i in the_l]
    for i, sf in zip(the_l, point):
        i.scale(sf)

    base = awsd(exp, the_l, exp_ab_cost, th_ab_cost)
    print("base:", base)

    for idx in range(len(the_l)):
        c = the_l[idx].copy()
        c.scale(dval + 1.0) #point[idx])
        n_the_l = the_l[:idx] + [c] + the_l[idx+1:]
        print("Aaa")
        grpart = (awsd(exp, n_the_l, exp_ab_cost, th_ab_cost) - base) / dval
        print("grpart", awsd(exp, n_the_l, exp_ab_cost, th_ab_cost))
        res.append(grpart)

    return res

import random
import numpy as np

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
        return checked_grad(EXP, THEs, point, exp_ab_cost, the_ab_cost)
    print("susp:", minimized_f([1.0, 0.8]))

    print(scipy.optimize.minimize(minimized_f, [15340.1] * len(THEs), jac=grad))

if __name__ == '__main__':
    EXP = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [1.0, 1.0])
    THE1 = IsoSpecPy.IsoDistribution(masses=[1.001, 2.0], probs = [0.2, 0.1])

    THE2 = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [1.0, 0.0])
    THEs = [THE1, THE2]

    #print(awsd(EXP, THEs, 0.5, 0.5))
    #print(empiric_gradient(EXP, THEs, [1.0, 1.0], 0.5, 0.5))
    #sys.exit(0)
    deconvolve(EXP, THEs, 0.5, 0.5)

#    EXP = random_spectrum()
#    THEs = [random_spectrum() for _ in range(10)]

    print(checked_grad(EXP, THEs, [0.1]*len(THEs), 0.5, 0.5))
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
