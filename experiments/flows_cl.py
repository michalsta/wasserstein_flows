import copy
import math
import networkx
import numpy as np
from pprint import pprint
from IsoSpecPy import IsoDistribution
from parameters import integerize, integerize_single, int_fact
from nx_print import nx_print

negative_penalty = 10000.0
negative_const_penalty = 1e10

def integerize_spectrum(spectrum):
    return [(label, integerize_single(prob)) for label, prob in spectrum]

def deiso(spectrum):
    if isinstance(spectrum, IsoDistribution):
        return list(zip(spectrum.masses, spectrum.probs))
    return spectrum

class AbWSDistCalc:
    def __init__(self, experimental, theoreticals, exp_ab_cost, the_ab_cost, distance_fun):
        self.G = networkx.DiGraph()

        self.experimental = integerize_spectrum(deiso(experimental))
        self.theoreticals = [copy.deepcopy(deiso(s)) for s in theoreticals]
        self.theoretical_p_to_probs = [{label:prob for label,prob in the_s} for the_s in self.theoreticals]

        exp_sum = sum(prob for _, prob in self.experimental)

        self.exp_ab_cost = integerize_single(exp_ab_cost)
        self.the_ab_cost = integerize_single(the_ab_cost)
        self.tot_ab_cost = self.exp_ab_cost + self.the_ab_cost

        for label, prob in self.experimental:
            self.G.add_edge('source', (label, -1), capacity = prob, weight = 0)
            self.G.add_edge((label, -1), 'middle', weight = self.exp_ab_cost)
            
        for idx, the_s in enumerate(self.theoreticals):
            for label, _ in the_s:
                self.G.add_edge((label, idx), 'sink', capacity = 0, weight = 0)
                self.G.add_edge('middle', (label, idx), weight = self.the_ab_cost)


        for elabel, _ in self.experimental:
            for idx, the_s in enumerate(self.theoreticals):
                for tlabel, _ in the_s:
                    dist = integerize_single(distance_fun(elabel, tlabel))
                    if dist < self.tot_ab_cost:
                        self.G.add_edge((elabel, -1), (tlabel, idx), weight = integerize_single(distance_fun(elabel, tlabel)))

        self.G.add_edge('source', 'middle', capacity = 0, weight = 0)
        self.G.add_edge('middle', 'sink', capacity = exp_sum, weight = 0)

        self.point = None
        self.realpoint = None

    def set_point(self, point):
        assert len(point) == len(self.theoreticals)
        print("Point:", point)

        point = tuple(point)
        if self.realpoint == point:
            return

        self.realpoint = point
        self.point = tuple(max(x, 0.0) for x in point)

        print(self.point)

        edges_view = self.G.edges
        tot_the_f = 0
        for idx, (coord, the_s) in enumerate(zip(self.point, self.theoreticals)):
            for label, prob in the_s:
                ip = integerize_single(coord * prob)
                edges_view[((label, idx), 'sink')]['capacity'] = ip
                tot_the_f += ip

        self.G.edges[("source", "middle")]["capacity"] = tot_the_f

        #nx_print(self.G)

        self.flow = networkx.max_flow_min_cost(self.G, 'source', 'sink')
        #pprint(self.flow)


    def check_point(self):
        try:
            self.flow
        except AttributeError:
            raise RuntimeError("You must set point before calling value() or gradient() method")
        

    def value(self):
        self.check_point()

        tot_cost = 0.0
        for start_n, neigh in self.flow.items():
            for end_n, fl_v in neigh.items():
                tot_cost += fl_v * self.G.edges[start_n, end_n]['weight']

        neg_penalty = sum(0.0 if x >= 0.0 else -negative_penalty*x for x in self.realpoint)

        val = neg_penalty + tot_cost/(int_fact*int_fact)
        print("Val:", val)
        return val


    def gradient(self):
        self.check_point()

        grad = np.zeros(len(self.theoreticals), dtype=float)
        for th_gr_lab in self.G.predecessors('sink'):
            if th_gr_lab != 'middle':
                label, idx = th_gr_lab
                if self.flow['middle'][th_gr_lab] > 0.0:
                    grad[idx] += self.the_ab_cost * self.theoretical_p_to_probs[idx][label]
                else:
                    delta = min(
                                    self.G[exp_p][th_gr_lab]['weight'] if th_gr_lab != 'middle'
                                    else math.inf
                                    for exp_p in self.G.predecessors(th_gr_lab)
                                )
                    grad[idx] += (delta - self.exp_ab_cost) * self.theoretical_p_to_probs[idx][label]

        neg_penalty = np.array([0.0 if x >= 0.0 else -negative_penalty for x in self.realpoint])
        grad = neg_penalty + grad/int_fact
        print("Grad:", grad)
        return grad

    def value_at(self, point):
        self.set_point(point)
        return self.value()

    def gradient_at(self, point):
        self.set_point(point)
        return self.gradient()

if __name__ == '__main__':

    EXP = [ 
            ((0.0, 0.0), 2.0),
            ((1.0, 2.0), 4.0),
          ]

    THE1 = [((0.0, 0.0), 1.0)]
    THE2 = [((1.0, 2.0), 1.0)]

    def euclidean_distance(point1, point2):
        return math.sqrt(sum((x_i - y_i)*(x_i - y_i) for x_i, y_i in zip(point1, point2)))

    wasserstein_distancer = AbWSDistCalc(EXP, [THE1, THE2], 500.0, 500.0, euclidean_distance)

    wasserstein_distancer.set_point([2.0, 5.0])
    print(wasserstein_distancer.value())
    print(wasserstein_distancer.gradient())

    print("------------------------")
    import scipy.optimize

    print(scipy.optimize.minimize(wasserstein_distancer.value_at, [733.0, 2467.0], jac=wasserstein_distancer.gradient_at, method='BFGS'))
