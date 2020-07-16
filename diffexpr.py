from collections import defaultdict, Counter
import networkx as nx
import json
import sys
import logging

from KEGGgraph import KEGGparser
from FileDict import FileDict


def optimal_greedy(graph, diff_expressed):
    main_ancestorz = Counter()
    diff_expressed_nodes = graph.nodes & diff_expressed

    all_children = defaultdict(set)

    for node in diff_expressed_nodes:
        diff_expr_descendants = nx.descendants(graph, node) & diff_expressed_nodes
        all_children[node].update(diff_expr_descendants)
        all_children[node].add(node)

    while all_children:
        whomax = max(all_children, key=lambda x: len(all_children[x]))
        main_ancestorz[whomax] += 1
        all_children = {key: all_children[key] - all_children[whomax] for key in all_children}
        all_children = {key: value for key, value in all_children.items() if value}

    return main_ancestorz


with open('config.json', 'r') as fo:
    config = json.load(fo)

logging.basicConfig(level=logging.INFO,
                    filename='diffexprproject.log', filemode='w',
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Load kegg.gs data
with open('data/kegg.gs.json', 'r') as fi:
    kegg_gs = json.load(fi)
print(len(kegg_gs), 'pathways in kegg.gs')
with open('data/kegg.gs.dise.json', 'r') as fi:
    kegg_gs_dise = json.load(fi)
print(len(kegg_gs_dise), 'pathways in kegg.gs.dise')

# Common dict of kegg pathways
kegg_pathways = {}
kegg_total_genes = set()
for k in kegg_gs:
    hsa, title = k.split(' ', maxsplit=1)
    kegg_pathways[hsa] = set(kegg_gs[k])
    kegg_total_genes.update(kegg_pathways[hsa])
for k in kegg_gs_dise:
    hsa, title = k.split(' ', maxsplit=1)
    kegg_pathways[hsa] = set(kegg_gs_dise[k])
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
# filtering intersection
intersections = {k: v for k, v in intersections.items() if v >= config['min diff expr genes in pathways']}
print('{} pathways left after cutoff by min {} genes'.format(len(intersections),
                                                             config['min diff expr genes in pathways']))

KP = KEGGparser(save_local=config['save KGML local'],
                local_kgml_dir=config['local KGML dir'],
                genes_only=config['genes only'])

supergraph = nx.DiGraph()

print('Looking in graphs separately')
diff_ancestors = defaultdict(set)
all_mega_ancestorz = Counter()
for pathway in intersections:
    G = KP[pathway]
    if G is None:
        continue

    supergraph.update(G)

    # nodes with diff expr genes in this graph
    diff_expressed_nodes = G.nodes & diff_expressed

    # for each node (from diff expr) making a set of  other diff expt gene nodes reachable from it
    for node in diff_expressed_nodes:
        diff_expr_descendants = nx.descendants(G, node) & diff_expressed_nodes
        diff_ancestors[node].update(diff_expr_descendants)

    all_mega_ancestorz += optimal_greedy(G, diff_expressed)

freq = {k: len(v) for k, v in diff_ancestors.items()}
freq = sorted(freq.items(), key=lambda x: x[1], reverse=True)

hsa_names = FileDict(file='data/hsa_names.tsv',
                     pattern='hsa:(?P<key>\d+)\t(?P<value>.*)')

n_top_genes = 10
print('Top {} genes:'.format(n_top_genes))
hsa_names.preload(*[gene for gene, fr in freq[:n_top_genes]])
for gene, fr in freq[:n_top_genes]:
    print('{:>4} : {:>7} : {}'.format(fr, gene, hsa_names[gene]))
print()

print('Top {} genes in optimal-greedy algorithm:'.format(n_top_genes))
for gene, fr in all_mega_ancestorz.most_common(n_top_genes):
    print('{:>4} : {:>7} : {}'.format(fr, gene, hsa_names[gene]))
print()

print('Looking in the supergraph')
supergraph.graph['name'] = 'SuperGraph'
print(nx.info(supergraph))
diff_ancestors = defaultdict(set)
diff_expressed_nodes = supergraph.nodes & diff_expressed
# for each node (from diff expr) making a set of  other diff expt gene nodes reachable from it
for node in diff_expressed_nodes:
    diff_expr_descendants = nx.descendants(supergraph, node) & diff_expressed_nodes
    diff_ancestors[node].update(diff_expr_descendants)

freq_dict = {k: len(v) for k, v in diff_ancestors.items()}
freq = sorted(freq_dict.items(), key=lambda x: x[1], reverse=True)

print('Top {} genes:'.format(n_top_genes))
hsa_names.preload(*[gene for gene, fr in freq[:n_top_genes]])
for gene, fr in freq[:n_top_genes]:
    print('{:>4} : {:>7} : {}'.format(fr, gene, hsa_names[gene]))

all_mega_ancestorz_supergraph = optimal_greedy(supergraph, diff_expressed)

print('Top {} genes in optimal-greedy algorithm:'.format(n_top_genes))
for gene, fr in all_mega_ancestorz_supergraph.most_common(n_top_genes):
    print('{:>4} : {:>7} : {}'.format(fr, gene, hsa_names[gene]))
print()

unsuper = supergraph.to_undirected()
components_sizes = Counter()
isolated_nodes = set()
for component in nx.connected_components(unsuper):
    x = len(component)
    components_sizes[x] += 1
    if x == 1:
        isolated_nodes.update(component)
print('Connected components of the supergraph:')
print(' size : number')
for size, number in components_sizes.items():
    print('{:>5} : {}'.format(size, number))
print('\n{} isolated diff expr genes'.format(len(isolated_nodes & diff_expressed)))

giant_component_nodes = max(nx.connected_components(unsuper), key=len)
giant_component = nx.subgraph(supergraph, giant_component_nodes)
giant_component.graph['name'] = 'giant component'
print('\nLooking in the giant component')
# nodes with diff expr genes in this graph
diff_expressed_nodes = giant_component_nodes & diff_expressed
first_diff_ancestors = set()
for node in diff_expressed_nodes:
    diff_ancestors = nx.ancestors(giant_component, node) & diff_expressed_nodes
    if not diff_ancestors:
        first_diff_ancestors.add(node)
print(len(first_diff_ancestors), 'diff expr genes in the giant component without diff expr ancestors')
first_freq = sorted(first_diff_ancestors, key=lambda x: freq_dict[x], reverse=True)
print('Top {} genes:'.format(n_top_genes))
hsa_names.preload(*first_freq[:n_top_genes])
for gene in first_freq[:n_top_genes]:
    print('{:>4} : {:>7} : {}'.format(freq_dict[gene], gene, hsa_names[gene]))
