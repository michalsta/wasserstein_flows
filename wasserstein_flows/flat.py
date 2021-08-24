import pathlib
from math import inf
from tqdm import tqdm
import numpy as np

def samesign(x, y):
    if x<0 and y<0:
        return True
    if y>=0 and x>= 0:
        return True
    return False


epsilon = 1e-8
debug = False
#debug = True

import cppyy

here = pathlib.Path(__file__).parent.resolve()
cppyy.include(here/'cmirror.cpp')

class FlatGraph:
    def __init__(self, exp, the_l, exp_ab_cost, the_ab_cost):

        self.exp_ab_cost = exp_ab_cost
        self.the_ab_cost = the_ab_cost
        self.tot_ab_cost = exp_ab_cost + the_ab_cost

        L = [(mass, prob, -1) for mass, prob in zip(exp.masses, exp.probs)]
        for idx, the in enumerate(the_l):
            L.extend((mass, prob, idx) for mass, prob in zip(the.masses, the.probs))

        L.sort()

        masses, probs, idxes = zip(*L)
        '''
        self.masses = cppyy.gbl.std.vector[float](self.masses)
        self.probs = cppyy.gbl.std.vector[float](self.probs)
        self.idxes = cppyy.gbl.std.vector[int](self.idxes)
        
        nconfs = len(self.masses)

        self.inline_flows = cppyy.gbl.std.vector[float](nconfs-1) #[0.0] * (nconfs-1) # Positive is to right, neg. to left
        self.costs = cppyy.gbl.std.vector[float]([self.masses[i+1] - self.masses[i] for i in range(nconfs-1)])

        self.directed_probs = cppyy.gbl.std.vector([prob if idx >= 0 else -prob for _, prob, idx in L])
        self.from_abyss_flows = cppyy.gbl.std.vector([-x for x in self.directed_probs])
'''
#        self.cmirror = cppyy.gbl.CMirror(cppyy.gbl.std.move(self.masses), cppyy.gbl.std.move(self.probs), cppyy.gbl.std.move(self.idxes))
        self.cmirror = cppyy.gbl.CMirror(masses, probs, idxes, exp_ab_cost, the_ab_cost)

        self.masses = self.cmirror.masses
        self.probs = self.cmirror.probs
        self.idxes = self.cmirror.idxes
        self.inline_flows = self.cmirror.G.inline_flows
        self.costs = self.cmirror.G.costs
        self.directed_probs = self.cmirror.G.directed_probs
        self.from_abyss_flows = self.cmirror.G.from_abyss_flows

#        for i in range(nconfs):
#            print("AAA", self.acceptable_modification_range(i))

        #self.self_verify()

    def __str__(self):
        ret = []
        return '\n'.join(self.peak_summary(i) for i in range(len(self)))
        for mass, prob, idx in zip(self.masses, self.directed_probs, self.idxes):
            ret.append(f"{idx}\t{mass}\t{prob}\t")
        return '\n'.join(ret)

    def __len__(self):
        return len(self.probs)

    def inflow_from_left(self, idx):
        if idx > 0:
            return self.inline_flows[idx-1]
        return 0.0

    def inflow_from_right(self, idx):
        if idx < len(self)-1:
            return -self.inline_flows[idx]
        return 0.0

    def self_verify(self):
        if not debug:
            return
        for i in range(len(self)):
            print(self.peak_summary(i))
            assert abs(self.inflow_from_left(i) + self.inflow_from_right(i) + self.from_abyss_flows[i] + self.directed_probs[i]) < epsilon
            if self.directed_probs[i] == 0.0:
                assert self.from_abyss_flows[i] == 0.0
            else:
                assert self.from_abyss_flows[i] * self.directed_probs[i] <= 0.0
                assert 0.0 <= abs(self.from_abyss_flows[i]) <= self.probs[i]
                assert -epsilon <= self.yankable_flow(i) + self.stuffable_flow(i) - abs(self.directed_probs[i]) <= epsilon
                assert self.stuffable_flow(i) >= 0.0
                assert self.yankable_flow(i) >= 0.0
                #a,b = self.acceptable_modification_range(i)
                #assert a <= 0.0 <= b
                #print(a, b, -self.yankable_flow(i))
                #assert a == -self.yankable_flow(i)
                

    def acceptable_modification_range(self, idx):
        return tuple(sorted((-self.from_abyss_flows[idx] - self.directed_probs[idx], -self.from_abyss_flows[idx])))

    def yankable_flow(self, idx):
        return self.cmirror.yankable_flow(idx)
        '''
        if self.directed_probs[idx] > 0.0:
            return -self.from_abyss_flows[idx]
        return -self.directed_probs[idx] - self.from_abyss_flows[idx]
        '''

    def stuffable_flow(self, idx):
        return self.cmirror.stuffable_flow(idx)
        '''
        if self.directed_probs[idx] > 0.0:
            return self.directed_probs[idx] + self.from_abyss_flows[idx]
        return self.from_abyss_flows[idx]
        '''

    def can_modify_into_point_flow(self, idx, delta):
        if delta >= 0.0:
            return delta <= self.stuffable_flow(idx)
        return -delta <= self.yankable_flow(idx)
#        x,y = self.acceptable_modification_range(idx)
#        return x - epsilon <= delta <= y + epsilon

    def can_send(self, src, tgt, howmuch):
        #print("SRC:", src, self.can_modify_into_point_flow(src, -howmuch))
        #print("TGT:", tgt, self.can_modify_into_point_flow(tgt, howmuch))
        #print("ACCMRTGT:", self.acceptable_modification_range(src))
        return src == tgt or (self.can_modify_into_point_flow(src, -howmuch) and self.can_modify_into_point_flow(tgt, howmuch))

    def send(self, src, tgt, howmuch):
        self.cmirror.G.send(src, tgt, howmuch)
        '''
        if src == tgt:
            return
        if tgt < src:
            return self.send(tgt, src, -howmuch)

        #print("src:", self.peak_summary(src))
        #print("tgt:", self.peak_summary(tgt))
        #print("Amount:", howmuch)
        if debug:
            assert self.can_send(src, tgt, howmuch)

        #print(self.from_abyss_flows[src], howmuch)
        self.from_abyss_flows[src] += howmuch
        #print(self.from_abyss_flows[src])
        self.from_abyss_flows[tgt] -= howmuch

        for i in range(src, tgt):
            self.inline_flows[i] += howmuch
'''
        if debug:
            self.self_verify()

    def total_cost(self):
        cost = sum(c * abs(f) for c, f in zip(self.costs, self.inline_flows))
        cost += sum(self.get_abyss_cost(idx)*abs(self.from_abyss_flows[idx]) for idx in range(len(self)))
        return cost

    def get_abyss_cost(self, idx):
        return self.cmirror.get_abyss_cost(idx)
        '''
        return self.exp_ab_cost if self.idxes[idx] == -1 else self.the_ab_cost
        '''

    def stuffing_into_abyss_cost(self, idx):
        return self.cmirror.stuffing_into_abyss_cost(idx)
        '''
        #print("SIA", self.peak_summary(idx))
        if self.from_abyss_flows[idx] > 0.0:
            return (-self.get_abyss_cost(idx), self.from_abyss_flows[idx])
        #print("CSDC", self.get_abyss_cost(idx), self.probs[idx] + self.from_abyss_flows[idx])
        return (self.get_abyss_cost(idx), self.probs[idx] + self.from_abyss_flows[idx])
        '''

    def yanking_out_of_abyss_cost(self, idx):
        return self.cmirror.yanking_out_of_abyss_cost(idx)
        '''
        if self.from_abyss_flows[idx] < 0.0:
            return (-self.get_abyss_cost(idx), -self.from_abyss_flows[idx])
        return (self.get_abyss_cost(idx), self.probs[idx] - self.from_abyss_flows[idx])
        '''

    def peak_summary(self, idx):
        return f"{idx}: spectrum no: {self.idxes[idx]}\tmass: {self.masses[idx]}\tdirprob: {self.directed_probs[idx]}\tfrom_ab_flow: {self.from_abyss_flows[idx]},\tstuffable: {self.stuffable_flow(idx)}\tyankable: {self.yankable_flow(idx)}"

    def delta_cost(self, src, tgt):
        '''The change in cost that will happen if we send flow from src to tgt, and how many units of flow can we send at that cost'''
        if debug:
            assert 0 <= src < len(self)
            assert 0 <= tgt < len(self)
        return self.cmirror.delta_cost(src, tgt)
        '''
        if src == tgt:
            return (0.0, 0.0)
        cost, flow = self.yanking_out_of_abyss_cost(src)
        c, f = self.stuffing_into_abyss_cost(tgt)
        cost += c
        flow = min(f, flow)
        #print("Flow:", flow)
        #assert self.can_send(src, tgt, flow)

        if src < tgt:
            for i in range(src, tgt):
                if self.inline_flows[i] >= 0.0:
                    cost += self.costs[i]
                else:
                    cost -= self.costs[i]
                    flow = min(flow, -self.inline_flows[i])
        else:
            for i in range(tgt, src):
                if self.inline_flows[i] <= 0.0:
                    cost += self.costs[i]
                else:
                    cost -= self.costs[i]
                    flow = min(flow, self.inline_flows[i])

        return (cost, flow)
        '''

    def is_optimal(self):
        for src in range(len(self)):
            if self.yankable_flow(src) > 0.0:
                for tgt in range(len(self)):
                    if self.stuffable_flow(tgt) > 0.0 and self.delta_cost(src, tgt)[0] < 0.0:
                        if debug:
                            print(self.peak_summary(src))
                            print(self.peak_summary(tgt))
                            print("not optimal:", self.delta_cost(src, tgt)[0])
                        return False
        return True

    def naive_optimize_once(self):
        for src in range(len(self)):
            for tgt in range(len(self)):
                while True: # in lieu of a do: ... until loop
                    dc, df = self.delta_cost(src, tgt)
                    #print(str(self))
                    #print(src, tgt, dc, df)
                    if dc < 0.0:
                        self.send(src, tgt, df)
                    else:
                        break

    def pushout_once(self, idx):
        return self.cmirror.pushout_once(idx)
        '''
        if self.yankable_flow(idx) == 0.0:
            return False
        dcost_left = 0.0
        fl_left = 0.0
        left_idx = idx-1
        while left_idx >= 0:
            if self.stuffable_flow(left_idx) > 0.0:
                dcost_left, fl_left = self.delta_cost(idx, left_idx)
                if (dcost_left < 0.0 and fl_left > 0.0) or dcost_left > self.tot_ab_cost:
                    break
            left_idx -= 1
        dcost_right = 0.0
        fl_right = 0.0
        right_idx = idx+1
        while right_idx < len(self):
            if self.stuffable_flow(right_idx) > 0.0:
                dcost_right, fl_right = self.delta_cost(idx, right_idx)
                if (dcost_right < 0.0 and fl_right > 0.0) or dcost_right > self.tot_ab_cost:
                    break
            right_idx += 1
#        print("Dcosts:", dcost_left, dcost_right)
        if dcost_left >= 0.0 and dcost_right >= 0.0:
            return False
        if dcost_left < dcost_right:
            #print(dcost_left, dcost_right, left_idx, idx, fl_left)
            self.send(idx, left_idx, fl_left)
        else:
            self.send(idx, right_idx, fl_right)
        return True
        '''

    def pushout_optimize(self):
        pushed_out = False
        for i in range(len(self)):
            while self.pushout_once(i):
                pushed_out = True
                pass
        return pushed_out

    def optimize(self):
        while self.pushout_optimize(): #not self.is_optimal():
            print("DSCA")
#            self.pushout_optimize()


if __name__ == "__main__":
    from test_spectra import *
    from flows_cl import AbWSDistCalc
    import time
    import sys
    i = 0
    while True:
        i += 1
        rerandomize(i)
        LEXP, LTHEs = rerandomize()
        print()
        print()
        #print(str(FG))
        #print(FG.peak_summary(0))
        #print(FG.peak_summary(1))
        #print(FG.delta_cost(1, 0))
    #    print("Cost:", FG.total_cost(), FG.is_optimal())
    #   FG.send(1, 0, 1.0)
        if debug:
            AW_time = time.time()
            AWDC = AbWSDistCalc(LEXP, LTHEs, 1.0, 1.0, lambda x,y: abs(x-y))
            constr_time = time.time() - AW_time
            print("Flow graph construction time:", constr_time)
            AW_time = time.time()
            rwsd = AWDC.value_at([1.0]*len(LTHEs))
            AW_time = time.time() - AW_time
            print("Real AWSD:", rwsd, "elapsed time:", AW_time)

        FG_time = time.time()
        FG = FlatGraph(LEXP, LTHEs, 1.0, 1.0)
        FG.optimize()
        if debug:
            assert FG.is_optimal()
        FG_time = time.time() - FG_time
        print("FG time:", FG_time)
        print("FG cost:", FG.total_cost())
        if debug:
            print(FG.total_cost(), round(rwsd))
            assert FG.total_cost() == round(rwsd)
        if not debug:
            sys.exit(0)
        FG.pushout_optimize()
        FG.pushout_optimize()
        FG.pushout_optimize()
    #    FG.naive_optimize_once()
        print("Cost:", FG.total_cost(), FG.is_optimal())
        assert FG.is_optimal()
    #    FG.naive_optimize_once()
    #    print("Cost:", FG.total_cost(), FG.is_optimal())
        FG.self_verify()
