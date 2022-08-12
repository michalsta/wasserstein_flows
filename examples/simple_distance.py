import numpy as np

from wasserstein_flows.flows_cl import MatrixDistance

in_m = np.matrix([[1,2],[3,4], [1,2]], dtype=np.uint32)
print("Input matrix:")
print(in_m)
WS = MatrixDistance(in_m)
print("Transport cost:")
print(WS.cost())
print("Transport scheme: [a_1, a_2, a_3, ..., a_n] : number a_i means that ith source is sending its output to a_i sink. If a_i == -1 then ith point is not sending its output anywhere (can happen with more sources than sinks)")
print(WS.flows())
