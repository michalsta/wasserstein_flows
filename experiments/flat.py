from math import inf

def samesign(x, y):
    if x<0 and y<0:
        return True
    if y>=0 and x>= 0:
        return True
    return False


epsilon = 1e-8

class FlatGraph:
    def __init__(self, exp, the_l, exp_ab_cost, the_ab_cost):

        self.exp_ab_cost = exp_ab_cost
        self.the_ab_cost = the_ab_cost

        L = [(mass, prob, -1) for mass, prob in zip(exp.masses, exp.probs)]
        for idx, the in enumerate(the_l):
            print(L)
            L.extend((mass, prob, idx) for mass, prob in zip(the.masses, the.probs))

        L.sort()

        self.masses, self.probs, self.idxes = zip(*L)
        
        nconfs = len(self.masses)

        self.inline_flows = [0.0] * (nconfs-1) # Positive is to right, neg. to left
        self.costs = [self.masses[i+1] - self.masses[i] for i in range(nconfs-1)]

        self.directed_probs = [prob if idx >= 0 else -prob for _, prob, idx in L]
        self.from_abyss_flows = [-x for x in self.directed_probs]

#        for i in range(nconfs):
#            print("AAA", self.acceptable_modification_range(i))

        self.self_verify()

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
        for i in range(len(self)):
            print(self.peak_summary(i))
            assert abs(self.inflow_from_left(i) + self.inflow_from_right(i) + self.from_abyss_flows[i] + self.directed_probs[i]) < epsilon
            if self.directed_probs == 0.0:
                assert self.from_abyss_flows[i] == 0.0
            else:
                assert self.from_abyss_flows[i] * self.directed_probs[i] <= 0.0
                assert 0.0 <= abs(self.from_abyss_flows[i]) <= self.probs[i]
                a,b = self.acceptable_modification_range(i)
                print(a,b)
                assert a <= 0.0 <= b
                

    def acceptable_modification_range(self, idx):
        return tuple(sorted((-self.from_abyss_flows[idx] - self.directed_probs[idx], -self.from_abyss_flows[idx])))

    def can_modify_abyss_flow(self, idx, delta):
        if delta == 0.0:
            return True
        x,y = self.acceptable_modification_range(idx)
        return x <= delta <= y

    def can_send(self, src, tgt, howmuch):
        return src == tgt or (self.can_modify_abyss_flow(src, howmuch) and self.can_modify_abyss_flow(tgt, -howmuch))

    def send(self, src, tgt, howmuch):
        if src == tgt:
            return
        if tgt < src:
            return self.send(tgt, src, -howmuch)

        assert self.can_send(src, tgt, howmuch)

        self.from_abyss_flows[src] += howmuch
        self.from_abyss_flows[tgt] -= howmuch

        for i in range(src, tgt):
            self.inline_flows[i] += howmuch

        self.self_verify()

    def total_cost(self):
        cost = sum(c * abs(f) for c, f in zip(self.costs, self.inline_flows))
        cost += sum(self.get_abyss_cost(idx)*abs(self.from_abyss_flows[idx]) for idx in range(len(self)))
        return cost

    def get_abyss_cost(self, idx):
        return self.exp_ab_cost if self.idxes[idx] == -1 else self.the_ab_cost

    def stuffing_into_abyss_cost(self, idx):
        if self.from_abyss_flows[idx] > 0.0:
            return (-self.get_abyss_cost(idx), self.from_abyss_flows[idx])
        return (self.get_abyss_cost(idx), self.probs[idx] + self.from_abyss_flows[idx])

    def yanking_out_of_abyss_cost(self, idx):
        if self.from_abyss_flows[idx] < 0.0:
            return (-self.get_abyss_cost(idx), -self.from_abyss_flows[idx])
        return (self.get_abyss_cost(idx), self.probs[idx] - self.from_abyss_flows[idx])

    def peak_summary(self, idx):
        return f"{idx}: {self.idxes[idx]}\t{self.masses[idx]}\tdirprob: {self.directed_probs[idx]}\tfrom_ab_flow: {self.from_abyss_flows[idx]}"

    def delta_cost(self, src, tgt):
        '''The change in cost that will happen if we send flow from src to tgt, and how many units of flow can we send at that cost'''
        assert 0 <= src < len(self)
        assert 0 <= tgt < len(self)
        if src == tgt:
            return (0.0, 0.0)
        cost, flow = self.yanking_out_of_abyss_cost(src)
        c, f = self.stuffing_into_abyss_cost(tgt)
        cost += c
        flow = min(f, flow)

        if src < tgt:
            for i in range(src, tgt):
                if self.inline_flows[i] >= 0:
                    cost += self.costs[i]
                else:
                    cost -= self.costs[i]
                    flow = min(flow, -self.inline_flows[i])
        else:
            for i in range(tgt, src):
                if self.inline_flows[i] <= 0:
                    cost += self.costs[i]
                else:
                    cost -= self.costs[i]
                    flow = min(flow, -self.inline_flows[i])

        return (cost, flow)

    def is_optimal(self):
        for src in range(len(self)):
            for tgt in range(len(self)):
                if self.delta_cost(src, tgt)[0] < 0.0:
                    return False
        return True

    def naive_optimize_once(self):
        for src in range(len(self)):
            for tgt in range(len(self)):
                while True: # in lieu of a do: ... until loop
                    dc, df = self.delta_cost(src, tgt)
                    print(str(self))
                    print(src, tgt, dc, df)
                    if dc < 0.0:
                        self.send(src, tgt, df)
                    else:
                        break

    def pushout_once(self, idx):
        if self.from_abyss_flows[idx] >= 0.0:
            return False
        dcost_left = 0.0
        fl_left = 0.0
        left_idx = idx-1
        while left_idx >= 0:
            dcost_left, fl_left = self.delta_cost(idx, left_idx)
            if dcost_left < 0.0:
                break
            left_idx -= 1
        dcost_right = 0.0
        fl_right = 0.0
        right_idx = idx+1
        while right_idx < len(self):
            dcost_right, fl_right = self.delta_cost(idx, right_idx)
            if dcost_right < 0.0:
                break
            right_idx += 1
        if dcost_left >= 0.0 and dcost_right >= 0.0:
            return False
        if dcost_left < dcost_right:
            self.send(idx, left_idx, fl_left)
        else:
            self.send(idx, right_idx, fl_right)
        return True


    def pushout_optimize(self):
        for i in range(len(self)):
            while self.pushout_once(i):
                pass


if __name__ == "__main__":
    from test_spectra import *
    i = 0
    while True:
        i += 1
        rerandomize(i)
        LEXP, LTHEs = rerandomize()
        FG = FlatGraph(LEXP, LTHEs, 1.0, 1.0)
        print(str(FG))
        print(FG.peak_summary(0))
        print(FG.peak_summary(1))
        print(FG.delta_cost(1, 0))
    #    print("Cost:", FG.total_cost(), FG.is_optimal())
    #    FG.send(1, 0, 1.0)
        print("Cost:", FG.total_cost(), FG.is_optimal())
        FG.pushout_optimize()
    #    FG.naive_optimize_once()
        print("Cost:", FG.total_cost(), FG.is_optimal())
        assert FG.is_optimal()
    #    FG.naive_optimize_once()
    #    print("Cost:", FG.total_cost(), FG.is_optimal())
        FG.self_verify()
