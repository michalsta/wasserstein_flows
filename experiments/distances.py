import math
import numpy as np

from IsoSpecPy import IsoDistribution

from flows import awsd
from flows_cl import AbWSDistCalc


def awsd_checked(exp, the_l, exp_ab_cost, th_ab_cost):
    val1 = awsd(exp, the_l, exp_ab_cost, th_ab_cost)

    DC = AbWSDistCalc(exp, the_l, exp_ab_cost, th_ab_cost, lambda x, y: abs(x-y))
    val2 = DC.value_at([1.0] * len(the_l))
    #val2 = exp.abyssalWassersteinDistance(IsoDistribution.LinearCombination(the_l, [1.0]*len(the_l)), (exp_ab_cost + th_ab_cost))
    #val3 = exp.abyssalWassersteinDistanceGrad(the_l, [1.0] * (1+len(the_l)), [1.0] * (1+len(the_l)), exp_ab_cost, th_ab_cost)

    if not math.isclose(val1, val2, rel_tol=1e-05):# and math.isclose(val1, val3, rel_tol=1e-05)):
        print("mismatching distances:", val1, val2)
        raise Exception()
    return val1

if __name__ == "__main__":
    from test_spectra import *
    awsd_checked(EXP, [THE1, THE2], 1.0, 1.0)
