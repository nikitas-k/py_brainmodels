#!/usr/bin/env python

import copy
import networkx as nx
import numpy as np

__author__ = "Nikitas Koussis"
__license__ = "GPL-3.0"
__email__ = "nikitas.koussis@gmail.com"

def linear_threshold(G, L, seeds, threshold=0.5):
    """Return the active nodes of each diffusion step by linear threshold model.
    Parameters
    ----------
    G : networkx graph
        Node weights
    L : networkx graph
        Length weights
    seeds: list of nodes
        The seed nodes of the graph
    Return
    ------
        i_edges : list of activated edges, as well as 'time' of activation
    Notes
    -----
    1. Each node is supposed to have an attribute "threshold".  If not, the
    default value is given (0.5).
    2. Each edge is supposed to have an attribute "influence".  If not, the
    default value is given (1/in_degree)
    3. This model does not currently diffuse over a certain number of steps -
    it simply diffuses until it cannot activate any more regions.
    4. 'Time' of diffusion is defined as the distance weight [L.edge(i,j)]
    divided by the fibre density weight [G.edge(i,j)]
    5. Primarily to be used for exploring undirected connectomes generated 
    from white matter tractography. Could also be used for other metrics, 
    such as multi-modal functional-structural networks.
    
    References
    ----------
    [1] Granovetter, Mark. Threshold models of collective behavior.
    The American journal of sociology, 1978.
    Examples
    --------
    >>> G = nx.Graph()
    >>> G.add_edges_from([(1,2), (1,3), (1,5), (2,1), (3,2), (4,2), (4,3), \
    >>>   (4,6), (5,3), (5,4), (5,6), (6,4), (6,5)])
    >>> L = nx.Graph()
    >>> L.add_edges_from([(1,2), (1,3), (1,5), (2,1), (3,2), (4,2), (4,3), \
    >>>   (4,6), (5,3), (5,4), (5,6), (6,4), (6,5)])
    >>> i_nodes = run_ltm.linear_threshold(G, L, [1], threshold=0.05)
"""

    if type(G) == nx.MultiGraph or type(G) == nx.MultiDiGraph:
        raise Exception( \
        "linear_threshold() is not defined for graphs with multiedges.")
    # make sure graphs are equal in size
    if len(list(G.nodes)) != len(list(L.nodes)):
        raise Exception("Graphs must have the same dimensions.")
        
    # make sure graphs are undirected
    if nx.is_undirected(G) is False or nx.is_undirected(L) is False:
        raise Exception("This function only works on undirected graphs.")
        
    # make sure the seeds are in the graph
    for s in seeds:
        if s not in G.nodes():
            raise Exception("seed", s, "is not in graph")

    # init thresholds
    for n in G.nodes():
        G.node[n]['threshold'] = threshold
        if G.node[n]['threshold'] > 1:
            raise Exception("node threshold:", G.node[n]['threshold'], \
            "cannot be larger than 1")

    # init influences
    deg = G.degree()
    for e in G.edges():
        G[e[0]][e[1]]['influence'] = 1.0 / deg[e[1]]
        if G[e[0]][e[1]]['influence'] > 1:
           raise Exception("edge influence:", G[e[0]][e[1]]['influence'], \
            "cannot be larger than 1")

    # perform diffusion
    A = copy.deepcopy(seeds)
    if steps <= 0:
        # perform diffusion until no more nodes can be activated
        return _diffuse_all(G, L, A, seeds)

def _diffuse_all(G, L, A, seeds):
    i_nodes = [ ]
    time = 0.0
    while True:
        A, i_edges, tau = _diffuse_one_round(G, L, A, time)
        if A == None:
            return i_nodes
        time += tau
        i_nodes.extend(i_edges)

def _diffuse_one_round(G, L, A, time):
    activated_nodes_of_this_round = set()
    len_old = len(A)
    last_active = A[-1]
    nbs = G.neighbors(last_active)
    i_edges = []
    taus = []
    for nb in nbs:
        if nb not in A:
            try:
                taus.append((nb, np.divide(L.edges[last_active, nb]['weight'], 
                G.edges[last_active, nb]['weight'])))   
                taus = sorted(taus, key = lambda x: x[1])
            except:
                continue
    for targ, tau in taus:
        actives = list(set(G.neighbors(targ)).intersection(set(A)))
        true_actives = []
        for active in actives:
            t = np.divide(L.edges[active, targ]['weight'], G.edges[active, targ]['weight'])
            if t <= tau:
                true_actives.append(active)
        if _influence_sum(G, true_actives, targ) >= G.node[targ]['threshold']:
            activated_nodes_of_this_round.add(targ)
            A.extend(list(activated_nodes_of_this_round))
            for a in true_actives:
                i_edges.append((a, targ, time+tau))
            return A, i_edges, tau

        else:
            continue

    if len(A) == len_old:
        return None, None, None




def _influence_sum(G, froms, to):
    influence_sum = 0.0
    for f in froms:
        influence_sum += G[f][to]['influence']
    return influence_sum


