from collections import defaultdict
import networkx as nx
import itertools
import matplotlib.pyplot as plt
import json
import sys
import logging

from KEGGgraph import KEGGparser

min_genes = 5  # min number of genes in the pathway to cutoff

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

# Load diff expressed genes
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
n_top_pathways = 5  # print top pathways
print('Top {}:'.format(n_top_pathways))
for key, value in sorted(intersections.items(), key=lambda x: x[1], reverse=True)[:n_top_pathways]:
    print('{:>3} : {:>10} : {}'.format(value, key, kegg_pathways_titles[key]))
print()

KP = KEGGparser(save_local=True, local_kgml_dir='/home/kirill/sources/R/Dif_expression_profiles_project/kegg_data/')

supergraph = nx.DiGraph()
removed_gene_nodes = set()

print('Looking in graphs separately')
diff_ascendants = defaultdict(set)
for pathway in intersections:
    G = KP[pathway]
    if G is None:
        continue

    supergraph.update(G)

    # nodes with diff expr genes in this graph
    diff_expressed_nodes = G.nodes & diff_expressed

    # for each node (from diff expr) making a set of  other diff expt gene nodes reachable from it
    for node in diff_expressed_nodes:
        diff_ascendants[node].update(nx.descendants(G, node) & diff_expressed_nodes)

freq = {k: len(v) for k, v in diff_ascendants.items()}
freq = sorted(freq.items(), key=lambda x: x[1], reverse=True)

n_top_genes = 20
print('Top {} genes:'.format(n_top_genes))
for gene, fr in freq[:n_top_genes]:
    print('{:>4} : {:>7} : {}'.format(fr, gene, supergraph.nodes[gene]['graphicalname']))
print()

print('Looking in the supergraph')
diff_ascendants = defaultdict(set)
diff_expressed_nodes = supergraph.nodes & diff_expressed

# for each node (from diff expr) making a set of  other diff expt gene nodes reachable from it
for node in diff_expressed_nodes:
    diff_ascendants[node].update(nx.descendants(supergraph, node) & diff_expressed_nodes)

freq = {k: len(v) for k, v in diff_ascendants.items()}
freq = sorted(freq.items(), key=lambda x: x[1], reverse=True)

print('Top {} genes:'.format(n_top_genes))
for gene, fr in freq[:n_top_genes]:
    print('{:>4} : {:>7} : {}'.format(fr, gene, supergraph.nodes[gene]['graphicalname']))
