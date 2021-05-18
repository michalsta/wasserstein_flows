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
            assert abs(self.inflow_from_left(i) + self.inflow_from_right(i) + self.from_abyss_flows[i] + self.directed_probs[i]) < epsilon
            if self.directed_probs == 0:
                assert self.from_abyss_flows[i] == 0
            else:
                assert self.from_abyss_flows[i] * self.directed_probs[i] <= 0
                assert 0.0 <= abs(self.from_abyss_flows[i]) <= self.probs[i]
                

    def acceptable_modification_range(self, idx):
        return tuple(sorted((self.from_abyss_flows[idx] + self.directed_probs[idx], -self.from_abyss_flows[idx])))

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
            self.send(tgt, src, -howmuch)

        assert self.can_send(src, tgt, howmuch)

        self.from_abyss_flows[src] += howmuch
        self.from_abyss_flows[tgt] -= howmuch

        if tgt < src:
            self.send(tgt, src, -howmuch)
        
        for i in range(src, tgt):
            self.inline_flows[i] += howmuch
        

    def delta_cost(self, src, tgt):
        '''The change in cost that will happen if we send flow from src to tgt, and how many units of flow can we send at that cost'''
        '''
        if src == tgt:
            return 0
        c = 0.0
        m = inf
        if self.from_abyss_flow[src] >= 0.0:
            c += self.
        m = min(from_abyss_flow
        if src < tgt:
            m = inf
            c = 0.0
            for i in range(src, tgt):
                if self.inline_flows[i] >= 0:
                    c += self.costs[i]
                else:
                    c -= self.costs[i]
                    m = min(
        '''

if __name__ == "__main__":
    from test_spectra import *
    FG = FlatGraph(EXP, THEs, 1.0, 1.0)
    print(str(FG))
    FG.send(0, 1, -1.0)
    FG.self_verify()
