"""
BAmodel.py
Module to implement the algorithms to generate complex networks.

I Manco    26/03/17

Classes: - System

Derived Classes (of System): - Preferential
                             - Random
                             - Random_Walk
                             
Functions: - create_ensemble()
           - get_ensemble()
           
Create_ensemble() writes to a file the results for an ensemble of 
networks generated through preferential attachment, random attachment or random
walk. A destination directory needs to be specified.
"""

import networkx as nx
import numpy as np
import random
import os
import Networks_Analysis as net
reload(net)

class System: 
    
    def __init__(self, n, N, m):  
        """Initialise graph.
        n -- number of initial nodes
        N -- number of final nodes
        m -- number of edges to add to new vertice
        """
        n = int(n)
        self.N = int(N)
        self.m = int(m)
        self.n0 = int(n)
        self.n = range(int(n))
        self.e = []
        self.G = nx.Graph()
        for i in range(n):
            if i == range(n)[-1]:   
                self.e.append([i, range(n)[0]])
            else:
                self.e.append([i, i+1])
    
    def nodes(self):
        return self.n
    
    def edges(self, node = None):
        if node == None:
            return self.e
        #return list of edges of specific node
        else:
            edges = []
            for edge in self.e:
                if edge[0] == node:
                    edges.append(edge)
                if edge[1] == node:
                    edges.append(edge)
            return edges
    
    def check_edges(self, n1, n2):
        if [n1, n2] in self.edges(n1):
            return False
        if [n2, n1] in self.edges(n1):
            return False
            
    def import_graph(self):
        #self.G.add_nodes_from(self.nodes()) 
        self.G.add_edges_from(self.edges()) 
    
    def draw_graph(self):
        self.import_graph()
        nx.draw(self.G)

class Preferential(System):
    
    def iteration(self):
        for i in range(self.N - self.n0):
            new_node = self.nodes()[-1]+1
            #random_edges = np.random.choice(len(self.edges()), self.m)
            edges_new_node = []
            while len(edges_new_node) < self.m:
                random_edges = np.random.choice(len(self.edges()))
                random_node = np.random.choice(self.edges()[random_edges])
                if random_node not in edges_new_node:
                    self.e.append([new_node, random_node])
                    edges_new_node.append(random_node)
            self.n.append(new_node)     

class Random(System):     
             
    def iteration(self):
        for i in range(self.N - self.n0):
            new_node = self.nodes()[-1]+1
            random_node = random.sample(range(self.n[-1]), self.m)
            for node in random_node:
                self.e.append([new_node, node])
            self.n.append(new_node) 

class Random_Walk(System):
    
    def __init__(self, n, N, m, L):
        """Initialise graph.
        n -- number of nodes
        N -- number of final nodes
        m -- number of edges to add to new vertice
        L -- length of random walk
        """
        System.__init__(self, n, N, m)
        self.L = L
        
    def find_neighbour(self, node):
        neighbours = []
        edges = np.array(self.edges(node))
        edges = np.concatenate(edges)
        for n_node in edges:
            if n_node != node:
                neighbours.append(n_node)
        return neighbours
            
    def iteration(self):
        L = self.L
        for i in range(self.N - self.n0):
            new_node = self.nodes()[-1]+1
            random_node = random.sample(range(self.n[-1]), self.m)
            if L == 0:
                for node in random_node:
                    self.e.append([new_node, node])
            else:
                i = 0
                while i < L:
                    new_target = []
                    for node in random_node:
                        neighbours = self.find_neighbour(node)
                        random_neighbour = np.random.choice(neighbours)
                        new_target.append(random_neighbour)
                    random_node = new_target
                    i = i + 1
                for node_to_connect in random_node:
                    self.e.append([new_node, node_to_connect])
            self.n.append(new_node)

def create_ensemble(m, N, attachment, runs, L = 0):
    """Write results to file.
    N          -- number of final nodes
    m          -- number of edges to add to new vertice
    L          -- length of random walk
    attachment -- type of attachment (string)
    """
    #n = number of nodes in initial graph G_0
    n = m + 1
    k = []
    max_degree = []
    new_dir = '/Users/Ilaria/Documents/Imperial/Third Year/Networks/Data/'+str(attachment)
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    os.chdir(new_dir)
    for run in range(runs):
        if attachment == 'preferential':
            network = Preferential(n, N, m)    
        if attachment == 'random':
            network = Random(n, N, m)    
        if attachment == 'random walk':
            network = Random_Walk(n, N, m, L)  
        network.iteration()
        network.import_graph()
        for value in network.G.degree().values():
            k.append(value)
        max_degree.append(np.max(network.G.degree().values()))
    if attachment == 'random walk':
        f1 = open('Max_Degree_'+str(N)+'_'+str(m)+str(L), 'w')
        f = open('Degree_'+str(N)+'_'+str(m)+str(L), 'w')
    else:
        f1 = open('Max_Degree_'+str(N)+'_'+str(m), 'w')
        f = open('Degree_'+str(N)+'_'+str(m), 'w')
    for degree in k:
        f.write(str(degree) + '\n')
    for max_k in max_degree:
        f1.write(str(max_k) + '\n')
    os.chdir('/Users/Ilaria/Documents/Imperial/Third Year/Networks/Code')
 
def get_ensemble(m, N, attachment, analysis, runs, L = 0):
    """Read results from text file.
    N          -- number of final nodes
    m          -- number of edges to add to new vertice
    L          -- length of random walk
    attachment -- type of attachment (string)
    """
    if analysis == 'm':
        new_dir = '/Users/Ilaria/Documents/Imperial/Third Year/Networks/Data_to_use/'+str(attachment)+'_m'
    if analysis == 'finite size':
        new_dir = '/Users/Ilaria/Documents/Imperial/Third Year/Networks/Data_to_use/'+str(attachment)+'_N'
    if analysis == 'none':
        new_dir = '/Users/Ilaria/Documents/Imperial/Third Year/Networks/Data_to_use/'+str(attachment)
    os.chdir(new_dir)
    if attachment == 'random walk':
        f = open('Degree_'+str(N)+'_'+str(m)+str(L), 'r')
        f1 = open('Max_Degree_'+str(N)+'_'+str(m)+str(L), 'r')
    else:
        f = open('Degree_'+str(N)+'_'+str(m), 'r')
        f1 = open('Max_Degree_'+str(N)+'_'+str(m), 'r')
    k = []
    max_k = []
    for line in f:
        k.append(int(line))
    for line in f1:
        max_k.append(int(line))
    os.chdir('/Users/Ilaria/Documents/Imperial/Third Year/Networks/Code')
    if m == 1:
        a = 1.3
    if m == 2:
        a = 1.3
    if m == 4:
        a = 1.3
    if m == 8:
        a = 1.3
    if m == 16:
        a = 1.3
    if m == 32:
        a = 1.1
    else:
        a = 1.2
    results = net.Ensemble(a, k, max_k, N, m, attachment, runs, L)
    return results