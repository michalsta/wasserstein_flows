from pprint import pprint
import networkx as nx


def nx_print(G):
    ret = "-------------------------------------\n"
    for node in sorted(G.nodes, key = str):
        ret += str(node) + ':\n'
        for succ in G.successors(node):
            ret += '\t' + str(succ) + '\t' + str(G.edges[node, succ]) + '\n'
    ret += "-------------------------------------\n"
    print(ret)
