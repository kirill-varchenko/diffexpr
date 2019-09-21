import networkx as nx
from xml.etree import ElementTree
from os import listdir
import itertools
import matplotlib.pyplot as plt
import json
import sys

import logging

import KEGGgraph

min_genes = 5 # min number of genes in the pathway to cutoff


logging.basicConfig(level=logging.INFO,
                    filename='diffexprproject.log', filemode='w',
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Load kegg.gs data
with open('kegg.gs.json', 'r') as fi:
    kegg_gs = json.load(fi)
print(len(kegg_gs), 'pathways in kegg.gs')
with open('kegg.gs.dise.json', 'r') as fi:
    kegg_gs_dise = json.load(fi)
print(len(kegg_gs_dise), 'pathways in kegg.gs.dise')

# Common dict of kegg pathways
kegg_pathways = {}
kegg_pathways_titles = {}
kegg_total_genes = set()
for k in kegg_gs:
    hsa, title = k.split(' ', maxsplit=1)
    kegg_pathways[hsa] = set(kegg_gs[k])
    kegg_pathways_titles[hsa] = title
    kegg_total_genes.update(kegg_pathways[hsa])
for k in kegg_gs_dise:
    hsa, title = k.split(' ', maxsplit=1)
    kegg_pathways[hsa] = set(kegg_gs_dise[k])
    kegg_pathways_titles[hsa] = title
    kegg_total_genes.update(kegg_pathways[hsa])
print('...with total {} genes\n'.format(len(kegg_total_genes)))

# Load diff exprenssed genes
diff_expressed = set()
with open('gse_16357_001.csv', 'r') as fi:
    for line in fi:
        if '/' not in line:
            diff_expressed.add(line.strip())
print(len(diff_expressed), 'diff expressed genes loaded\n')

# intersection between genes in diff expr dataset and pathway genes from KEGG
intersections = {k: len(kegg_pathways[k] & diff_expressed) for k in kegg_pathways}
intersections = {k: v for k, v in intersections.items() if v >= min_genes}
print('{} pathways left after cutoff by min {} genes'.format(len(intersections), min_genes))
n_top = 5 # print top pathways
print('Top {}:'.format(n_top))
for key, value in sorted(intersections.items(), key=lambda x: x[1], reverse=True)[:n_top]:
    print('{:>3} : {:>10} : {}'.format(value, key, kegg_pathways_titles[key]))
print()

print('diff expr genes:')
mapped_diff_genes = diff_expressed & kegg_total_genes
print(len(mapped_diff_genes), 'mapped')
unmapped_diff_genes = diff_expressed - kegg_total_genes
print(len(unmapped_diff_genes), 'not mapped\n')

supergraph = nx.DiGraph(name='SuperGraph')
removed_gene_nodes = set()

for pathway in intersections:
    G = KEGGgraph.getgraph(pathway)
    #G = KEGGgraph.parseKGML2(xmlpath)

    # remove isolated nodes
    # unG = G.to_undirected()
    # for subG in nx.connected_components(unG):
    #     if len(subG) == 1:
    #         G.remove_nodes_from(subG)
    #         removed_nodes.update(subG)
    #
    # supergraph.update(G)

print(nx.info(supergraph))
# unG = supergraph.to_undirected()
# comp = []
# for subG in nx.connected_components(unG):
#     comp.append((len(subG), len(subG & diff_expressed)))
#
# comp.sort(reverse=True, key=lambda x: x[0])
# print(comp)
#
# print(len(removed_nodes & diff_expressed))
