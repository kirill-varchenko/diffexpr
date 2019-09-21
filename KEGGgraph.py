from os.path import exists
import requests
import networkx as nx
from xml.etree import ElementTree
import itertools
import logging

path_to_kegg_data = '/home/kirill/sources/R/Dif_expression_profiles_project/kegg_data/'
link_to_kgml = 'https://www.kegg.jp/kegg-bin/download?entry={}&format=kgml'


def getgraph(pathway, save_local=True):
    """Return the pathway graph either from local data or from kegg site"""
    logger = logging.getLogger(__name__)
    xmlpath = path_to_kegg_data + pathway + '.xml'
    if exists(xmlpath):
        tree = ElementTree.parse(xmlpath)
        root = tree.getroot()
    else:
        try:
            r = requests.get(link_to_kgml.format(pathway), timeout=5)
            r.raise_for_status()
            root = ElementTree.fromstring(r.text)
        except requests.exceptions.HTTPError:
            logger.warning('Unable to download pathway xml: {}'.format(pathway))
            return None
        except requests.exceptions.ConnectTimeout:
            logger.warning('Unable to download pathway xml: {}'.format(pathway))
            return None
        except ElementTree.ParseError:
            logger.warning('Unable to parse pathway xml: {}'.format(pathway))
            return None
        except Exception:
            logger.warning('Unknown error getting pathway xml: {}'.format(pathway))
            return None

        if save_local:
            with open(xmlpath, 'w') as fo:
                fo.write(r.text)

    return parseKGML2(root)


def parseKGML(xmlpath):
    """Парсит только гены и связи между ними"""
    tree = ElementTree.parse(xmlpath)
    root = tree.getroot()

    G = nx.DiGraph(**root.attrib)

    entry_ids = {}

    for entry in root.findall('entry'):
        if entry.get('type') == 'gene':
            attr = entry.attrib
            local_id = attr.pop('id')
            entry_ids[local_id] = list(map(lambda x: x[4:], attr['name'].split(' ')))
            attr['graphicalname'] = entry[0].get('name')
            for name in entry_ids[local_id]:
                G.add_node(name, **attr)
        elif entry.get('type') == 'group':
            local_id = entry.get('id')
            entry_ids[local_id] = []
            for component in entry:
                if component.tag == 'component':
                    entry_ids[local_id].extend(entry_ids[component.get('id')])

    # connection patterns of rectangles (gene products)
    for relation in root.findall('relation'):
        subtypes = [subtype.get('name') for subtype in relation]
        for node1, node2 in itertools.product(entry_ids[relation.get('entry1')], entry_ids[relation.get('entry2')]):
            G.add_edge(node1, node2, subtypes=subtypes)

    return G


def parseKGML2(root):
    """Парсит всё"""
    G = nx.DiGraph(**root.attrib)

    entry_ids = {}

    for entry in root.findall('entry'):
        attr = entry.attrib
        local_id = attr.pop('id')
        entry_ids[local_id] = list(map(lambda x: x.replace('hsa:', ''), attr['name'].split(' ')))

        if entry.get('type') in ['gene', 'compound']:
            attr['graphicalname'] = entry[0].get('name')
            for name in entry_ids[local_id]:
                G.add_node(name, **attr)
        elif entry.get('type') == 'group':
            entry_ids[local_id] = []
            for component in entry:
                if component.tag == 'component':
                    entry_ids[local_id].extend(entry_ids[component.get('id')])

    # connection patterns of rectangles (gene products)
    for relation in root.findall('relation'):
        subtypes = [subtype.get('name') for subtype in relation]
        for node1, node2 in itertools.product(entry_ids[relation.get('entry1')], entry_ids[relation.get('entry2')]):
            G.add_edge(node1, node2, subtypes=subtypes, xmltype='relation')

    # connection patterns of circles (chemical compounds)
    for reaction in root.findall('reaction'):
        gene_id = reaction.get('id')
        reagents = {'substrate': [], 'product': []}
        for reagent in reaction:
            reagents[reagent.tag].extend(entry_ids[reagent.get('id')])
        for substr, gene in itertools.product(reagents['substrate'], entry_ids[gene_id]):
            G.add_edge(substr, gene, type=reaction.get('type'), xmltype='reaction')
        for gene, prod in itertools.product(entry_ids[gene_id], reagents['product']):
            G.add_edge(gene, prod, type=reaction.get('type'), xmltype='reaction')

    return G
