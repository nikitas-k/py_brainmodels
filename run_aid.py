from epydemic import *

import pandas as pd

import matplotlib
%matplotlib inline
%config InlineBackend.figure_format = 'svg'
import matplotlib.pyplot as plt
import seaborn
import pickle
import AID

import networkx
N = 500                 # order (number of nodes) of the network
kmean = 5                 # mean node degree
phi = (kmean + 0.0) / N   # probability of attachment between two nodes chosen at random

# create the network
g = networkx.erdos_renyi_graph(N, phi)

pInfected = 0.2
pInfect = 0.2
pRecover = 0.05
pRemoveEdge = 1.0
pAddEdge = 1.0
maxTime = 50000
model = AID()

params = dict()
params[AID.P_INFECTED] = pInfected
params[AID.P_INFECT] = pInfect
params[AID.P_RECOVER] = pRecover
params[AID.P_REMOVE_EDGE] = pRemoveEdge
params[AID.P_ADD_EDGE] = pAddEdge
params[AID.MAX_TIME] = maxTime

# run the simulation
f = SynchronousDynamics(model, g)
sync = f.set(params).run()

# save the results for later
with open('sync.pickle', 'wb') as handle:
    pickle.dump(sync, handle)


