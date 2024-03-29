import copy
import math
import networkx
import numpy as np
from pprint import pprint
from IsoSpecPy import IsoDistribution
from wasserstein_flows.parameters import integerize, integerize_single, int_fact

negative_penalty = 10000.0
negative_const_penalty = 1e10
do_prints = False

def integerize_spectrum(spectrum):
    return [(label, integerize_single(prob)) for label, prob in spectrum]

def deiso(spectrum):
    if isinstance(spectrum, IsoDistribution):
        return list(zip(spectrum.masses, spectrum.probs))
    return spectrum

class AbWSDistCalc:
    def __init__(self, experimental, theoreticals, exp_ab_cost, the_ab_cost, distance_fun = None):
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


        if distance_fun is None:
            peaks = [(label, prob, -1) for label, prob in self.experimental]
            for idx, the_s in enumerate(self.theoreticals):
                peaks.extend((label, prob, idx) for label, prob in the_s)
            peaks.sort()
            for idx in range(len(peaks)-1):
                self.G.add_edge((peaks[idx][0], peaks[idx][2]), (peaks[idx+1][0], peaks[idx+1][2]), weight = integerize_single(peaks[idx+1][0]-peaks[idx][0]))

        else:
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
        if do_prints: print("Point:", point)

        point = tuple(point)
        if self.realpoint == point:
            return

        self.realpoint = point
        self.point = tuple(max(x, 0.0) for x in point)

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

        if do_prints: print("Val:", val)

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

        if do_prints: print("Grad:", grad)

        return grad

    def value_at(self, point):
        self.set_point(point)
        return self.value()

    def gradient_at(self, point):
        self.set_point(point)
        return self.gradient()


class MatrixDistance:
    def __init__(self, dist_matrix):

        assert dist_matrix.dtype in [np.uint32, np.uint64], "Input matrin must be in one of numpy's unsigned integer types (np.uint32, np.uint64)"
        self.G = networkx.DiGraph()

        x_dim, y_dim = dist_matrix.shape
        self.x_dim = x_dim
        self.y_dim = y_dim

        for x in range(x_dim):
            for y in range(y_dim):
                self.G.add_edge(f"src{x}", f"tgt{y}", capacity=1, weight=dist_matrix[x,y])

        for x in range(x_dim):
            self.G.add_edge("source", f"src{x}", capacity=1, weight=0)

        for y in range(y_dim):
            self.G.add_edge(f"tgt{y}", "sink", capacity=1, weight=0)

        self.mfmc = networkx.max_flow_min_cost(self.G, 'source', 'sink')

    def cost(self):
        tot_cost = 0.0
        for start_n, neigh in self.mfmc.items():
            for end_n, fl_v in neigh.items():
                tot_cost += fl_v * self.G.edges[start_n, end_n]['weight']
        return tot_cost

    def flows(self):
        ret = np.full(self.x_dim, -1)
        for start_n, neigh in self.mfmc.items():
            if start_n[:3] != 'src':
                continue
            start_idx = int(start_n[3:])
            for end_n, fl_v in neigh.items():
                if end_n[:3] != 'tgt':
                    continue
                end_idx = int(end_n[3:])
                if fl_v > 0:
                    assert ret[start_idx] == -1 and fl_v == 1, "Something is wrong, this shouldn't happen"
                    ret[start_idx] = end_idx
        return ret

#if __name__ == '__main__':
if False:
    EXP = [ 
            ((0.0, 0.0), 2.0),
            ((1.0, 2.0), 4.0),
          ]

    THE1 = [((0.0, 0.0), 1.0)]
    THE2 = [((1.0, 2.0), 1.0)]

    def euclidean_distance(point1, point2):
        return math.sqrt(sum((x_i - y_i)*(x_i - y_i) for x_i, y_i in zip(point1, point2)))

    wasserstein_distancer = AbWSDistCalc(EXP, [THE1, THE2], 5.0, 5.0, euclidean_distance)

    wasserstein_distancer.set_point([2.0, 5.0])
    print(wasserstein_distancer.value())
    print(wasserstein_distancer.gradient())

    print("------------------------")
    import scipy.optimize

    print(scipy.optimize.minimize(wasserstein_distancer.value_at, [12.0, 15.0], jac=wasserstein_distancer.gradient_at, method='SLSQP'))


if __name__ == '__main__':
    in_m = np.matrix([[1,2],[3,4], [1,2]], dtype=np.uint32)
    WS = MatrixDistance(in_m)
    print(WS.cost())
    print(WS.flows())