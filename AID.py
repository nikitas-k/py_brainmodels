# Adapted from SIS Compartment Model from epydemic package
#
# epydemic Copyright (C) 2017 Simon Dobson 
#
# Copyright (C) 2020 Nikitas Koussis 
# 
# This file is modified from epydemic, epidemic network simulations in Python.
#
# This code is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with epydemic. If not, see <http://www.gnu.org/licenses/gpl.html>.



class AID(CompartmentedModel):

    '''The Susceptible-Infected-Susceptible :term:`compartmented model of disease`.
    Susceptible nodes are infected by infected neighbours, and recover back
    to the susceptible state (which allows future re-infection, unlike for
    :class:`SIR`). 
    Adaptive variant "Adaptive Information Diffusion Model"
    (c.f. Guo, D., Trajanovski, S., van de Bovenkamp, R., Wang, H., & Van Mieghem, P. (2013). 
    Epidemic threshold and topological structure of susceptible-infectious-susceptible 
    epidemics in adaptive networks. Physical Review E, 88(4), 042802.)
    
    Allows for adding neighbouring edges when a node has "information" i.e. infected
    and breaking links from susceptible nodes that aren't infected i.e. "no information".'''

    # the model parameters
    P_INFECTED = 'pInfected'           #: Parameter for probability of initially being infected.
    P_INFECT = 'pInfect'               #: Parameter for probability of infection on contact.
    P_RECOVER = 'pRecover'             #: Parameter for probability of recovery (returning to susceptible).
    P_REMOVE_EDGE = 'pRemoveEdge'      #: Parameter for probability of deleting a link given two susceptible nodes (default 1.0)
    P_ADD_EDGE = 'pAddEdge'            #: Parameter for probability of creating a link given one infected node (default 1.0)
    MAX_TIME = 'maxTime'               #: Parameter for the maximum sim time override setting (default 20000).
    
    # the possible dynamics states of a node for SIS dynamics
    SUSCEPTIBLE = 'S'         #: Compartment for nodes susceptible to infection.
    INFECTED = 'I'            #: Compartment for nodes infected.

    # the edges at which dynamics can occur
    SI = 'SI'                 #: Edge able to transmit infection.
    SS = 'SS'                 #: Edge able to be deleted if unused.

    def __init__( self ):
        super(AID, self).__init__()

    def build( self, params ):
        '''Build the SIS model.

        :param params: the model parameters'''
        pInfected = params[self.P_INFECTED]
        pInfect = params[self.P_INFECT]
        pRecover = params[self.P_RECOVER]
        pRemoveEdge = params[self.P_REMOVE_EDGE]
        maxTime = params[self.MAX_TIME]
        
        if maxTime != None:
            self.setMaximumTime(maxTime)
        
        self.a = AddDelete()

        self.addCompartment(self.SUSCEPTIBLE, 1 - pInfected)
        self.addCompartment(self.INFECTED, pInfected)

        self.trackEdgesBetweenCompartments(self.SUSCEPTIBLE, self.INFECTED, name = self.SI)
        self.trackEdgesBetweenCompartments(self.SUSCEPTIBLE, self.SUSCEPTIBLE, name = self.SS)
        self.trackNodesInCompartment(self.INFECTED)

        self.addEventPerElement(self.INFECTED, pRecover, self.recover)
        self.addEventPerElement(self.SI, pInfect, self.infect)
        
        self.addEventPerElement(self.SS, pRemoveEdge, self.delete_unused_edge)
        self.addEventPerElement(self.INFECTED, pAddEdge, self.add_neighboring_edges)

    def infect( self, t, g, e ):
        '''Perform an infection event. This changes the compartment of
        the susceptible-end node to :attr:`INFECTED`. It also marks the edge
        traversed as occupied.

        :param t: the simulation time (unused)
        :param g: the network
        :param e: the edge transmitting the infection, susceptible-infected'''
        (n, m) = e
        try:
            self.changeCompartment(g, n, self.INFECTED)
            self.markOccupied(g, e)
        except:
            return
        

    def recover( self, t, g, n ):
        '''Perform a recovery event. This changes the compartment of
        the node back to :attr:`SUSCEPTIBLE`, allowing re-infection.

        :param t: the simulation time (unused)
        :param g: the network
        :param n: the node'''
        try:
            self.changeCompartment(g, n, self.SUSCEPTIBLE)
        except:
            return
        
    def delete_unused_edge( self, t, g, e ):
        '''Delete edges that contain no information.'''
        (n, m) = e
        self.a.deleteEdge(g, n, m)
        
    def add_neighboring_edges( self, t, g, n ):
        '''Add edges from parents that contain information.'''
        self.a.addNeighborEdges(g, n)
