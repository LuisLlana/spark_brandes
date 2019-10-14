#!/usr/bin/python3

# -*- coding: utf-8 -*-

"""
Original Surce:
  https://networkx.github.io/documentation/development/_modules/networkx/algorithms/centrality/betweenness.html#edge_betweenness_centrality

Modified by
Cristobal Pareja Flores <cpareja@ucm.es>
Luis Llana Díaz <llana@ucm.es>
"""

#    Copyright (C) 2004-2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

"""
Betweenness centrality measures.

References
----------
.. [1]  A Faster Algorithm for Betweenness Centrality.
   Ulrik Brandes,
   Journal of Mathematical Sociology 25(2):163-177, 2001.
   http://www.inf.uni-konstanz.de/algo/publications/b-fabc-01.pdf
.. [2] Ulrik Brandes: On Variants of Shortest-Path Betweenness
   Centrality and their Generic Computation.
   Social Networks 30(2):136-145, 2008.
   http://www.inf.uni-konstanz.de/algo/publications/b-vspbc-08.pdf
.. [3] Ulrik Brandes and Christian Pich:
   Centrality Estimation in Large Networks.
   International Journal of Bifurcation and Chaos 17(7):2303-2318, 2007.
   http://www.inf.uni-konstanz.de/algo/publications/bp-celn-06.pdf
.. [4] Aric Hagberg, Dan Schult, Pieter Swart, 2016
   https://networkx.github.io/documentation/development/_modules/networkx/algorithms/centrality/betweenness.html
"""

"""
Cambios:

(1) The functions
    _accumulate_bacis
    _accumulate_endpoints
    _accumulate_edges
have been substituted by similar functions that compute the contributions
to centrality by without making the accumulation.

This change allows that the contributions of each node or edge could be
calculated in parallel.
"""

from heapq import heappush, heappop
from itertools import count
import networkx as nx
import random
import functools

# Obsolete name

def betweenness_centrality(G, k=None, normalized=True, weight=None,
                           endpoints=False,
                           seed=None):
    return node_betweenness_centrality(G, k, normalized, weight,
                           endpoints,seed)

def node_betweenness_centrality(G, k=None, normalized=True, weight=None,
                           endpoints=False,
                           seed=None):

    betweenness = dict.fromkeys(G, 0.0)  # b[v]=0 for v in G
    if k is None:
        nodes = G
    else:
        random.seed(seed)
        nodes = random.sample(G.nodes(), k)

    def node_betweenness_parameterized(G, weight, endpoints):
        def f(s):
            if weight is None:  # use BFS
                S, P, sigma = _single_source_shortest_path_basic(G, s)
            else:  # use Dijkstra's algorithm
                S, P, sigma = _single_source_dijkstra_path_basic(G, s, weight)
            # Calculate new contribution:
            if endpoints:
                res = _accumulate_endpoints__contribution_from_node(S, P, sigma, s)
            else:
                res = _accumulate_basic__contribution_from_node(S, P, sigma, s)
            return res
        return f

    betweenness_contr_due_to_a_node = node_betweenness_parameterized(G, weight, endpoints)

    """
    for s in nodes:
        # single source shortest paths
        new_contribution = process_a_node(s)
        update_centrality(betweenness, new_contribution)
    """
    all_contributions = map(betweenness_contr_due_to_a_node, nodes)
    contribution_together = functools.reduce(add_contributions, all_contributions)
    update_centrality(betweenness, contribution_together)
    # rescaling
    rescaled_betweenness = _rescale(betweenness, len(G),
                           normalized=normalized,
                           directed=G.is_directed(),
                           k=k)
    return rescaled_betweenness

"""
def _accumulate_basic(betweenness, S, P, sigma, s):
    new_contrib = _accumulate_basic__contribution_from_node(S, P, sigma, s)
    update_centrality(betweenness, new_contrib)
    return betweenness
"""

def edge_betweenness_centrality(G, normalized=True, weight=None):
    betweenness = dict.fromkeys(G, 0.0)               # b[v]=0 for v in G
    betweenness.update(dict.fromkeys(G.edges(), 0.0)) # b[e]=0 for e in G.edges()

    def edge_betweenness_parameterized(G, weight):
        def f(s):
            if weight is None:  # use BFS
                S, P, sigma = _single_source_shortest_path_basic(G, s)
            else:  # use Dijkstra's algorithm
                S, P, sigma = _single_source_dijkstra_path_basic(G, s, weight)
            # accumulation
            return _accumulate_edges__contribution_from_node(betweenness, S, P, sigma, s)
        return f

    betweenness_contr_due_to_an_edge = edge_betweenness_parameterized(G, weight)

    """
    for s in G:
        # single source shortest paths
        new_contribution = process_a_node(s)
        update_centrality(betweenness, new_contribution)
    """
    all_contributions = map(betweenness_contr_due_to_an_edge, G)
    contribution_together = functools.reduce(add_contributions, all_contributions)
    update_centrality(betweenness, contribution_together)
    # rescaling
    for n in G:  # remove nodes to only return edges
        del betweenness[n]
    betweenness = _rescale_e(betweenness, len(G),
                             normalized=normalized,
                             directed=G.is_directed())
    return betweenness

# obsolete name

def edge_betweenness(G, normalized=True, weight=None):
    return edge_betweenness_centrality(G, normalized, weight)


# helpers for betweenness centrality

def _single_source_shortest_path_basic(G, s):
    S = []
    P = {}
    for v in G:
        P[v] = []
    sigma = dict.fromkeys(G, 0.0)    # sigma[v]=0 for v in G
    D = {}
    sigma[s] = 1.0
    D[s] = 0
    Q = [s]
    while Q:   # use BFS to find shortest paths
        v = Q.pop(0)
        S.append(v)
        Dv = D[v]
        sigmav = sigma[v]
        for w in G[v]:
            if w not in D:
                Q.append(w)
                D[w] = Dv + 1
            if D[w] == Dv + 1:   # this is a shortest path, count paths
                sigma[w] += sigmav
                P[w].append(v)  # predecessors
    return S, P, sigma


def _single_source_dijkstra_path_basic(G, s, weight='weight'):
    # modified from Eppstein
    S = []
    P = {}
    for v in G:
        P[v] = []
    sigma = dict.fromkeys(G, 0.0)    # sigma[v]=0 for v in G
    D = {}
    sigma[s] = 1.0
    push = heappush
    pop = heappop
    seen = {s: 0}
    c = count()
    Q = []   # use Q as heap with (distance,node id) tuples
    push(Q, (0, next(c), s, s))
    while Q:
        (dist, _, pred, v) = pop(Q)
        if v in D:
            continue  # already searched this node.
        sigma[v] += sigma[pred]  # count paths
        S.append(v)
        D[v] = dist
        for w, edgedata in G[v].items():
            vw_dist = dist + edgedata.get(weight, 1)
            if w not in D and (w not in seen or vw_dist < seen[w]):
                seen[w] = vw_dist
                push(Q, (vw_dist, next(c), v, w))
                sigma[w] = 0.0
                P[w] = [v]
            elif vw_dist == seen[w]:  # handle equal paths
                sigma[w] += sigma[v]
                P[w].append(v)
    return S, P, sigma

# ----------------------------------------------------------------------------

"""
def _accumulate_basic(betweenness, S, P, sigma, s):
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            betweenness[w] += delta[w]
    return betweenness
"""

def _accumulate_basic__contribution_from_node(S, P, sigma, s):
    r"""Compute the contributions from each node to betweenness,
    considering the paths from a single source, s.
    This functions doesn't updates the parameter betweenness."""
    contributions = {}
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            add_to_dictionary(contributions, w, delta[w])
    return contributions

"""
def _accumulate_basic(betweenness, S, P, sigma, s):
    new_contrib = _accumulate_basic__contribution_from_node(S, P, sigma, s)
    update_centrality(betweenness, new_contrib)
    return betweenness
"""

# ----------------------------------------------------------------------------

"""
def _accumulate_endpoints(betweenness, S, P, sigma, s):
    betweenness[s] += len(S) - 1
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            betweenness[w] += delta[w] + 1
    return betweenness
"""

def _accumulate_endpoints__contribution_from_node(S, P, sigma, s):
    r"""Compute the contributions from each node to betweenness,
    considering the paths from a single source, s
    and taking in account endpoints.
    This functions doesn't updates the parameter betweenness."""
    contributions = {s: len(S) - 1}
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            add_to_dictionary(contributions, w, delta[w] + 1)
    return contributions

"""
def _accumulate_endpoints(betweenness, S, P, sigma, s):
    new_contrib = _accumulate_endpoints__contribution_from_node(S, P, sigma, s)
    update_centrality(betweenness, new_contrib)
    return betweenness
"""

# ----------------------------------------------------------------------------

"""
def _accumulate_edges(betweenness, S, P, sigma, s):
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            c = sigma[v] * coeff
            if (v, w) not in betweenness:
                betweenness[(w, v)] += c
            else:
                betweenness[(v, w)] += c
            delta[v] += c
        if w != s:
            betweenness[w] += delta[w]
    return betweenness
"""

def _accumulate_edges__contribution_from_node(betweenness, S, P, sigma, s):
    r"""Compute the contributions from each edge to betweenness,
    considering the paths from a single source, s.
    This functions doesn't updates the parameter betweenness."""
    contributions = {}
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            c = sigma[v] * coeff
            if (v, w) in betweenness:
                add_to_dictionary(contributions, (v, w), c)
            else:
                add_to_dictionary(contributions, (w, v), c)
            delta[v] += c
        if w != s:
            add_to_dictionary(contributions, w, delta[w])
    return contributions

"""
def _accumulate_edges(betweenness, S, P, sigma, s):
    new_contrib = _accumulate_edges__contribution_from_node(betweenness, S, P, sigma, s)
    update_centrality(betweenness, new_contrib)
    return betweenness
"""

# ----------------------------------------------------------------------------
# Helper functions

def add_to_dictionary(dictionary, key, value):
    if key in dictionary:
        dictionary[key] += value
    else:
        dictionary[key] = value

def add_contributions(dictionary_a, dictionary_b):
    new_dicc = dictionary_a.copy()
    for key in dictionary_b:
        add_to_dictionary(new_dicc, key, dictionary_b[key])
    return new_dicc

def update_centrality(dictionary, keys_values):
    for elem, value in keys_values.items():
        add_to_dictionary(dictionary, elem, value)

# ----------------------------------------------------------------------------

def _rescale(betweenness, n, normalized, directed=False, k=None):
    if normalized is True:
        if n <= 2:
            scale = None  # no normalization b=0 for all nodes
        else:
            scale = 1.0 / ((n - 1) * (n - 2))
    else:  # rescale by 2 for undirected graphs
        if not directed:
            scale = 1.0 / 2.0
        else:
            scale = None
    if scale is not None:
        if k is not None:
            scale = scale * n / k
        for v in betweenness:
            betweenness[v] *= scale
    return betweenness


def _rescale_e(betweenness, n, normalized, directed=False):
    if normalized is True:
        if n <= 1:
            scale = None  # no normalization b=0 for all nodes
        else:
            scale = 1.0 / (n * (n - 1))
    else:  # rescale by 2 for undirected graphs
        if not directed:
            scale = 1.0 / 2.0
        else:
            scale = None
    if scale is not None:
        for v in betweenness:
            betweenness[v] *= scale
    return betweenness

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

def graph_from_file(archivo, is_list_adyac=True, is_directed=False):
    """
    If adyac=True:
        Each line in the file contains a vertex (say an origin),
        and the set of adjacent vertices to this origin
    adyac=False
        Each line in the file contains a couple of
        adjacent vertices, origin and target
    """
    my_graph = nx.DiGraph() if is_directed else nx.Graph()
    archi = open(archivo,'r', encoding='latin1')
    if is_list_adyac:
        for line in archi:
            vertices_list = map(int, line.strip().split())
            my_graph.add_nodes_from(vertices_list)
            origen, targets = vertices_list[0], vertices_list[1:]
            my_graph.add_edges_from([(origen, target) for target in targets])
    else:
        for line in archi:
            origen, target = map(int, line.strip().split())
            my_graph.add_nodes_from([origen, target])
            my_graph.add_edge(origen, target)
    archi.close()
    return my_graph

def show_graph(my_graph):
    print("--------------------------------------------------------")
    print("  - Nodes: ", sorted(my_graph.nodes()))
    print("  - Edges: ", sorted(my_graph.edges()))
    for node, neighbors in my_graph.adjacency_iter():
        print(node, neighbors)
    print("--------------------------------------------------------")

def demo(archivo_grafo, is_list_adyac=True, is_directed=False):
    G = graph_from_file(archivo_grafo, is_list_adyac, is_directed)
    #show_graph(G)
    centr1 = betweenness_centrality(G)
    print(centr1)
    centr2 = edge_betweenness_centrality(G)
    print()
    print(centr2)
    #show_graph(G)

def estudiar(archivo_grafo, archivo_resultados, is_list_adyac=True, is_directed=False):
    G = graph_from_file(archivo_grafo, is_list_adyac, is_directed)
    #show_graph(G)
    f = open(archivo_resultados, "w")
    centr1 = betweenness_centrality(G)
    f.write(str(centr1))
    f.write("\n\n")
    centr2 = edge_betweenness_centrality(G)
    f.write(str(centr2))
    f.close()


import time

def demo_karate():
    estudiar("./datasets/karate/karate-edges.txt", "./datasets/karate/karate_spark.txt", is_list_adyac=False)

def demo_astro():
    estudiar("./datasets/ca-AstroPh/CA-AstroPh-edges.txt", "./datasets/ca-AstroPh/CA-AstroPh-edges_spark.txt", is_list_adyac=False, is_directed=True)

def demo_miserables():
    estudiar("./datasets/miserables/lesmis-edges.txt", "./datasets/miserables/lesmis-edges_spark.txt", is_list_adyac=False)

def demo_wiki_vote():
    estudiar("./datasets/wiki-vote/Wiki-Vote-edges.txt", "./datasets/wiki-vote/Wiki-Vote-edges_spark.txt", is_list_adyac=False, is_directed=True)

def lanzar_con_crono(esta_demo):
    inicio = time.time()
    print("Inicio: ", time.ctime())
    esta_demo()
    final = time.time()
    print("Final: ", time.ctime())
    print("Tiempo transcurrido: ", final-inicio, "segundos.")

"""
--------------------------------------------------------
  - Nodes:  ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
  - Edges:  [('0', '8'), ('1', '2'), ('1', '3'), ('1', '7'), ('1', '9'), ('2', '4'), ('3', '4'), ('3', '5'), ('5', '6'), ('7', '6'), ('8', '2'), ('8', '9')]
0 {'8': {}}
8 {'0': {}, '2': {}, '9': {}}
1 {'3': {}, '7': {}, '9': {}, '2': {}}
3 {'1': {}, '4': {}, '5': {}}
7 {'1': {}, '6': {}}
9 {'1': {}, '8': {}}
2 {'1': {}, '8': {}, '4': {}}
4 {'2': {}, '3': {}}
5 {'3': {}, '6': {}}
6 {'5': {}, '7': {}}
--------------------------------------------------------
{'0': 0.0, '8': 0.24537037037037038, '1': 0.46759259259259256, '3': 0.23148148148148143, '7': 0.1388888888888889, '9':
0.12037037037037038, '2': 0.2592592592592592, '4': 0.0648148148148148, '5': 0.05555555555555555, '6': 0.027777777777777776}
{('0', '8'): 0.2, ('8', '2'): 0.23333333333333334, ('8', '9'): 0.15925925925925924, ('1', '3'): 0.22222222222222218, ('1',
'7'): 0.26666666666666666, ('1', '9'): 0.23333333333333334, ('1', '2'): 0.22592592592592597, ('3', '4'): 0.14814814814814814,
('3', '5'): 0.2, ('7', '6'): 0.15555555555555556, ('2', '4'): 0.15555555555555556, ('5', '6'): 0.08888888888888889}
--------------------------------------------------------
  - Nodes:  ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
  - Edges:  [('0', '8'), ('1', '2'), ('1', '3'), ('1', '7'), ('1', '9'), ('2', '4'), ('3', '4'), ('3', '5'), ('5', '6'), ('7', '6'), ('8', '2'), ('8', '9')]
0 {'8': {}}
8 {'0': {}, '2': {}, '9': {}}
1 {'3': {}, '7': {}, '9': {}, '2': {}}
3 {'1': {}, '4': {}, '5': {}}
7 {'1': {}, '6': {}}
9 {'1': {}, '8': {}}
2 {'1': {}, '8': {}, '4': {}}
4 {'2': {}, '3': {}}
5 {'3': {}, '6': {}}
6 {'5': {}, '7': {}}
--------------------------------------------------------
"""

def demo_arboles():
    print("--------------------------------------------------------------")
    g1 = {"vertices": [0, 1], "edges": [(0, 1)]}
    print("--- n = ", len(g1["vertices"]))
    print(g1)
    G1 = graph_ad_hoc(g1)
    print(edge_betweenness_centrality(G1))
    print("--------------------------------------------------------------")
    g1 = {"vertices": [0, 1, 2], "edges": [(0, 1), (1, 2)]}
    print("--- n = ", len(g1["vertices"]))
    print(g1)
    G1 = graph_ad_hoc(g1)
    print(edge_betweenness_centrality(G1))
    print("--------------------------------------------------------------")
    g1 = {"vertices": [0, 1, 2, 3], "edges": [(0, 1), (1, 2), (2, 3)]}
    print("--- n = ", len(g1["vertices"]))
    print(g1)
    G1 = graph_ad_hoc(g1)
    print(edge_betweenness_centrality(G1))
    print("--------------------------------------------------------------")
    g1 = {"vertices": [0, 1, 2, 3], "edges": [(0, 1), (1, 2), (1, 3)]}
    print("--- n = ", len(g1["vertices"]))
    print(g1)
    G1 = graph_ad_hoc(g1)
    print(edge_betweenness_centrality(G1))
    print("--------------------------------------------------------------")
    g1 = {"vertices": [0, 1, 2, 3, 4], "edges": [(0, 1), (0, 2), (0, 3), (0, 4)]}
    print("--- n = ", len(g1["vertices"]))
    print(g1)
    G1 = graph_ad_hoc(g1)
    print(edge_betweenness_centrality(G1))
    print("--------------------------------------------------------------")
    g1 = {"vertices": [0, 1, 2, 3, 4], "edges": [(0, 1), (0, 2), (0, 3), (3, 4)]}
    print("--- n = ", len(g1["vertices"]))
    print(g1)
    G1 = graph_ad_hoc(g1)
    print(edge_betweenness_centrality(G1))
    print("--------------------------------------------------------------")
    g1 = {"vertices": [0, 1, 2, 3, 4], "edges": [(0, 1), (1, 2), (2, 3), (3, 4)]}
    print("--- n = ", len(g1["vertices"]))
    print(g1)
    G1 = graph_ad_hoc(g1)
    print(edge_betweenness_centrality(G1))
    print("--------------------------------------------------------------")


def estudiar_nodes(archivo_grafo, archivo_resultados, is_list_adyac=True, is_directed=False):
    G = graph_from_file(archivo_grafo, is_list_adyac, is_directed)
    f = open(archivo_resultados, "w")
    time_ini = time.time()
    centr = betweenness_centrality(G)
    time_end = time.time()
    f.write('tiempo {0}\n'.format(time_ini - time_end))
    for k in sorted(centr.keys()):
        f.write('{0}: {1}\n'.format(k, centr[k]))
    f.close()

import functools

def estudiar_edges(archivo_grafo, archivo_resultados, is_list_adyac=True, is_directed=False):
    G = graph_from_file(archivo_grafo, is_list_adyac, is_directed)
    f = open(archivo_resultados, "w")
    time_ini = time.time()
    centr = edge_betweenness_centrality(G)
    time_end = time.time()
    f.write('tiempo {0}\n'.format(time_ini - time_end))
    for k in sorted(centr.keys()):
        if type(k) is tuple:
            f.write('{0}: {1}\n'.format(k, centr[k]))
    f.write("\n\n")
    f.close()

def demo_karate_nodes():
    estudiar_nodes("/home/luis/inves/betweeness/datasets/karate/karate-edges.txt", "/home/luis/inves/betweeness/datasets/karate/karate_cris_nodes.txt", is_list_adyac=False)

def demo_karate_edges():
    estudiar_edges("/home/luis/inves/betweeness/datasets/karate/karate-edges.txt", "/home/luis/inves/betweeness/datasets/karate/karate_cris_edges.txt", is_list_adyac=False)


def demo_astro_nodes():
    estudiar_nodes("/home/luis/inves/betweeness/datasets/ca-AstroPh/CA-AstroPh-edges.txt", "/home/luis/inves/betweeness/datasets/ca-AstroPh/CA-AstroPh-cris-nodes.txt", is_list_adyac=False)

def demo_astro_edges():
    estudiar_edges("/home/luis/inves/betweeness/datasets/ca-AstroPh/CA-AstroPh-edges.txt", "/home/luis/inves/betweeness/datasets/ca-AstroPh/CA-AstroPh-cris-edges.txt", is_list_adyac=False)



demo_astro_nodes()
demo_astro_edges()
