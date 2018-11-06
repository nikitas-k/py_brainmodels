#!/usr/bin/env python

import copy
import networkx as nx
import numpy as np

__author__ = "Nikitas Koussis"
__license__ = "GPL-3.0"
__email__ = "nikitas.koussis@gmail.com"

def SIR_model(G, L, seeds, beta=0.05, gamma=0.005):
    """Return the active nodes of each diffusion step by SIR model.
    Parameters
    ----------
    G : networkx graph
        Node weights
    L : networkx graph
        Length weights
    seeds: list of nodes
        The seed nodes of the graph
	beta: float
	    The 'normal' probability of activation of that node given one neighbour
    gamma: float
	    Flat probability of a node being 'removed' i.e. deactivated.
    Return
    ------
        i_edges : list of activated edges, as well as 'time' of activation
    Notes
    -----
    1. Each node is supposed to have an attribute "beta".  If not, the
    default value is given (0.05).
    2. This model does not currently diffuse over a certain number of steps -
    it simply diffuses until it cannot activate any more regions.
    3. 'Time' of diffusion is defined as the distance weight [L.edge(i,j)]
    divided by the weighting of the graph G between each node.
    References
    ----------
    [1] Kermack and A. McKendrick, “A Contribution to the Mathematical Theory of Epidemics,”
	Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character, 
	vol. 115, no. 772, pp. 700–721, Aug. 1927
    Examples
    --------
    >>> G = nx.Graph()
    >>> G.add_edges_from([(1,2), (1,3), (1,5), (2,1), (3,2), (4,2), (4,3), \
    >>>   (4,6), (5,3), (5,4), (5,6), (6,4), (6,5)])
    >>> L = nx.Graph()
    >>> L.add_edges_from([(1,2), (1,3), (1,5), (2,1), (3,2), (4,2), (4,3), \
    >>>   (4,6), (5,3), (5,4), (5,6), (6,4), (6,5)])
    >>> layers = run_sir.SIR_model(G, L, [1], beta=0.05, gamma=0.005)
"""

    if type(G) == nx.MultiGraph or type(G) == nx.MultiDiGraph:
        raise Exception( \
        "linear_threshold() is not defined for graphs with multiedges.")
    if len(list(G.nodes)) != len(list(L.nodes)):
        raise Exception("Graphs must have the same dimensions.")

    # make sure the seeds are in the graph
    for s in seeds:
        if s not in G.nodes():
            raise Exception("seed", s, "is not in graph")

    # perform diffusion
    A = copy.deepcopy(seeds)
    if steps <= 0:
        # perform diffusion until no more nodes can be activated
        return _diffuse_all(G, L, A, seeds)

    for n in G.nodes():
	G.node[n]['beta'] = beta
	if G.node[n]['beta'] > 1:
	    raise Exception("beta:", beta, \
			    "cannot be larger than 1")

def _diffuse_all(G, L, A, seeds):
    i_nodes = [ ]
    r = []
    time = 0.0
    while True:
        A, i_edges, r_, tau = _diffuse_one_round(G, L, A, r, time, gamma)
        if A == None:
            return i_nodes, list(set(r))
        time += tau
        i_nodes.extend(i_edges)
	r.extend(r_)

def _diffuse_one_round(G, L, A, r, time, gamma):
    activated_nodes_of_this_round = set()
    len_old = len(A)
    try:
	last_active = A[-1]
    except:
	return None, None, None, None
    nbs = G.neighbors(last_active)
    i_edges = []
    taus = []
    r_ = []
    eventp = np.random.random_sample()
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
            try:
                t = np.divide(L.edges[active, targ]['weight'], G.edges[active, targ]['weight'])
            except:
                continue
            if t <= tau:
		true_actives.append(active)
        if eventp < G.node[targ]['beta'] * len(true_actives):
            activated_nodes_of_this_round.add(targ)
            A.extend(list(activated_nodes_of_this_round))
            for a in true_actives:
                i_edges.append((a, targ, time+tau))
        else:
            continue
	
    for nd in A:
	if eventp < gamma:
	    r_.append(nd)
	    A.remove(nd)

    if len(A) == len_old:
        return None, None, None
    
    return A, i_edges, r_, tau

