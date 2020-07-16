from collections import defaultdict, Counter
import networkx as nx
import json
import sys
import logging
from tqdm import tqdm
from itertools import combinations
import pandas as pd
import math


def write_cytoscape(G, filename):
    cs = nx.cytoscape_data(G)
    with open(filename, 'w') as fo:
        fo.write(json.dumps(cs))


def giant_component(G):
    giant_component_nodes = max(nx.connected_components(G.to_undirected()), key=len)
    return nx.subgraph(G, giant_component_nodes)


def collapse(G, nodes, new_id, **kwargs):
    all_successors = set()
    all_predecessors = set()
    for node in nodes:
        if node in G.nodes:
            all_successors.update(G.successors(node))
            all_predecessors.update(G.predecessors(node))
    all_successors = all_successors - set(nodes)
    all_predecessors = all_predecessors - set(nodes)
    unfrozen = nx.DiGraph(G)
    if all_successors or all_predecessors:
        unfrozen.remove_nodes_from(nodes)
        unfrozen.add_node(new_id, **kwargs)
        for node in all_successors:
            unfrozen.add_edge(new_id, node)
        for node in all_predecessors:
            unfrozen.add_edge(node, new_id)
    return unfrozen


def is_directed_module(G, X):
    all_pred = set()
    all_succ = set()
    for x in X:
        all_pred.add(frozenset(G.predecessors(x)) - X)
        all_succ.add(frozenset(G.successors(x)) - X)
    return len(all_pred) <= 1 and len(all_succ) <= 1


def find_modules(G):
    all_sets = set()
    for u in G:
        all_sets.add(frozenset(G.predecessors(u)))
        all_sets.add(frozenset(G.successors(u)))
    modules = []
    for s in all_sets:
        if len(s) > 1 and is_directed_module(G, s):
            modules.append(s)
    return modules
