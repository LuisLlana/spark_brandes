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


from heapq import heappush, heappop
from itertools import count
import networkx as nx
import random
import functools
import time
import operator




# Obsolete name

def betweenness_centrality(G, k=None, normalized=True, weight=None,
                           endpoints=False,
                           seed=None):
    return node_betweenness_centrality(G, k, normalized, weight,
                           endpoints,seed)

def node_betweenness_parameterized(G, weight, endpoints):
    def f(s):
        if weight is None:  # use BFS
            S, P, sigma = _single_source_shortest_path_basic(G.value, s)
        else:  # use Dijkstra's algorithm
            S, P, sigma = _single_source_dijkstra_path_basic(G.value, s, weight)
        # Calculate new contribution:
        if endpoints:
            return _accumulate_endpoints__contribution_from_node(S, P, sigma, s)
        else:
            return _accumulate_basic__contribution_from_node(S, P, sigma, s)
    return f

def node_betweenness_centrality(G, spark_context,
                                k=None, normalized=True, weight=None,
                                endpoints=False,
                                seed=None):
    brG = spark_context.broadcast(G)
    ncores = int(spark_context.getConf().get('total-executor-cores'))
    betweenness_contr_due_to_a_node = node_betweenness_parameterized(brG, weight, endpoints)

    paralNodes = spark_context.parallelize(G.nodes_iter(), ncores)
    betweenness = paralNodes.flatMap(betweenness_contr_due_to_a_node).\
                  reduceByKey(operator.add).collectAsMap()
    # Creo que el reduce no se se paraleliza, hay que buscar alternativas
    _rescale(betweenness, len(G),
             normalized=normalized,
             directed=G.is_directed(),
             k=k)
    return betweenness

"""
def _accumulate_basic(betweenness, S, P, sigma, s):
    new_contrib = _accumulate_basic__contribution_from_node(S, P, sigma, s)
    update_centrality(betweenness, new_contrib)
    return betweenness
"""

def edge_betweenness_parameterized(G, weight):
    def f(s):
        if weight is None:  # use BFS
            S, P, sigma = _single_source_shortest_path_basic(G.value, s)
        else:  # use Dijkstra's algorithm
            S, P, sigma = _single_source_dijkstra_path_basic(G.value, s, weight)
        # accumulation
        return _accumulate_edges__contribution_from_node(G.value, S, P, sigma, s)
    return f


def edge_betweenness_centrality(G, spark_context, normalized=True, weight=None):
    brG = spark_context.broadcast(G)
    ncores = int(spark_context.getConf().get('total-executor-cores'))
    betweenness_contr_due_to_an_edge = edge_betweenness_parameterized(brG, weight)
    paralNodes = spark_context.parallelize(G.nodes_iter(), numSlices=ncores)
    betweenness = paralNodes.flatMap(betweenness_contr_due_to_an_edge).reduceByKey(operator.add).collectAsMap()

    # Creo que el reduce no se se paraleliza, hay que buscar alternativas
    # rescaling
    _rescale_e(betweenness, len(G),
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
0        if w != s:
            betweenness[w] += delta[w]
    return betweenness
"""

def _accumulate_basic__contribution_from_node(S, P, sigma, s):
    r"""Compute the contributions from each node to betweenness,
    considering the paths from a single source, s.
    This functions doesn't updates the parameter betweenness."""
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            yield w, delta[w]

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
    yield s, len(S) - 1
    delta = dict.fromkeys(S, 0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            delta[v] += sigma[v] * coeff
        if w != s:
            yield w, delta[w] + 1

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

def _accumulate_edges__contribution_from_node(G, S, P, sigma, s):
    r"""Compute the contributions from each edge to betweenness,
    considering the paths from a single source, s.
    This functions doesn't updates the parameter betweenness."""
    edges = set(G.edges())
    delta = dict.fromkeys(S, 0.0)
    while S:
        w = S.pop()
        coeff = (1.0 + delta[w]) / sigma[w]
        for v in P[w]:
            c = sigma[v] * coeff
            if (v, w) in edges:
                yield (v, w), c
            else:
                yield (w, v), c
            delta[v] += c
        # if w != s:
        #     add_to_dictionary(contributions, w, delta[w])
    return

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
    new_dict = dictionary_a.copy()
    for key in dictionary_b:
        add_to_dictionary(new_dict, key, dictionary_b[key])
    return new_dict

# def add_contributions(dictionary_a, dictionary_b):
#     new_dicc = dictionary_a.copy()
#     for key in dictionary_b:
#         add_to_dictionary(new_dicc, key, dictionary_b[key])
#     return new_dicc

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
