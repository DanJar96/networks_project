import os, pickle, itertools
from time import sleep

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns


import networkx as nx
import pandas as pd
import numpy as np
import random
from random import sample
from typing import DefaultDict
from tqdm import tqdm
from pathlib import Path
import math

from itertools import count, permutations, combinations, chain
import copy
from IPython import display
import time
from varname.helpers import Wrapper


from PIL import Image
import imageio

## TO-DO
#- Add legends
#- Add plot on RHS that shows progression of number of infected individuals (lineplot)
#- Check why p-intra and p-inter prints wrong in the actual GIF.
#- Check the whole functioning to make sure disease accounting is correct.
# PSEUDO-CODE
# 1. Establish DJDM-class
class MJDM:
    
    def __init__(self,n : int):
        
        # Initializing number of nodes
        self.n = n
        
        # Initializing the nx graph
        self.graph = nx.Graph()
        
        # Adding nodes to graph
        self.graph.add_nodes_from(range(self.n+1))
        
        # Adding edges initializatio flag
        self.edge_initialization_flag = 0
        
        # Setting time period attribute
        self.timeperiod = 0
        
        # Creating the disease dynamic attributes
        susceptible = {node : {"susceptible" : 1} for node in self.graph.nodes}
        infected = {node : {"infected" : 0} for node in self.graph.nodes}
        recovered = {node : {"recovered" : 0} for node in self.graph.nodes}
        dead = {node : {"dead" : 0} for node in self.graph.nodes}
        time_since_infection = {node : {"time_since_infection" : 0} for node in self.graph.nodes}
        time_since_recovered = {node : {"time_since_recovered" : 0} for node in self.graph.nodes}
        
        # Creating community structure attributes
        community = {node : {"community" : 0} for node in self.graph.nodes}
        
        # Adding positioning for position tracking throughout graphs when called
        self.pos = 0
        self.posfix = 0
        
        # Adding them to the graph
        nx.set_node_attributes(self.graph, susceptible)
        nx.set_node_attributes(self.graph, infected)
        nx.set_node_attributes(self.graph, recovered)
        nx.set_node_attributes(self.graph, dead)
        nx.set_node_attributes(self.graph, time_since_infection)
        nx.set_node_attributes(self.graph, time_since_recovered)

        
    def infect(self,p_infection,no_infected,infection_period,recovery_period):
        """
        This function defines the disease in terms of how long a person
        is infected for and how long they stay in the recovered state for 
        before becoming susceptible again.
        """
        # Defining the disease in terms of what the probability of infection is,
        # given that two members of society meet each other with one person
        # being susceptible and the other being infected.
        self.p_infection = p_infection

        # Defining the disease in terms of how long a person is infected for
        # and how long they stay in recovered state before becoming 
        # susceptible again
        self.infection_period = infection_period
        self.recovery_period = recovery_period
        

        # Creating node status history accounting lists
        self.nodehistory = {}
        self.nodehistory['susceptible'] = []
        self.nodehistory['infected'] = []
        self.nodehistory['recovered'] = []

        
        # for node in self.graph.nodes:
        #     self.nodehistory[node] = {}
        #     self.nodehistory[node]['susceptible'] = []
        #     self.nodehistory[node]['infected'] = []
        #     self.nodehistory[node]['recovered'] = []


        # Randomly sampling two nodes that will become infected
        infected_nodes = random.sample(list(self.graph.nodes),no_infected)
        
        # Infecting the nodes
        for node in infected_nodes:
            # Current state
            self.graph.nodes[node]['susceptible'] = 0
            self.graph.nodes[node]['infected'] = 1
            self.graph.nodes[node]['recovered'] = 0
            
            # Start infection timer
            self.graph.nodes[node]['time_since_infection'] += 1

        # Recording current states...
        self.nodehistory['susceptible'].append(sum(nx.get_node_attributes(self.graph,'susceptible').values())/len(self.graph.nodes))
        self.nodehistory['infected'].append(sum(nx.get_node_attributes(self.graph,'infected').values())/len(self.graph.nodes))
        self.nodehistory['recovered'].append(sum(nx.get_node_attributes(self.graph,'recovered').values())/len(self.graph.nodes))

        # for node in self.graph.nodes:
        #     self.nodehistory[node]['susceptible'].append(self.graph.nodes[node]['susceptible'])
        #     self.nodehistory[node]['infected'].append(self.graph.nodes[node]['infected'])
        #     self.nodehistory[node]['recovered'].append(self.graph.nodes[node]['recovered'])

                
    def assign_community(self, no_communities, c_intra : float, c_inter : float):
        """
        This function assigns the probabilities of connecting with someone from
        inside your community and someone from the outside of your community.

        Additionally, it sets the community structure so that we can use it for
        disease propagation.
        """
        
        # Assigning to model, for later potential use..
        self.p_intra = c_intra / (len(self.graph.nodes) / no_communities) #approx
        self.c_intra = c_intra
        self.p_inter = float(c_inter / len(self.graph.nodes))
        
        # Assigning the number of communities
        self.no_communities = no_communities
        self.communities = list(range(1,no_communities+1))
        
        # Assign each node to a community
        for node in self.graph.nodes:
            # Sample a number from the list of communities
            self.graph.nodes[node]['community'] = int(random.sample(self.communities,1)[0])
    
        # If we have not yet initialized any community edges, do so now.
        if self.edge_initialization_flag == 0:

            # Accounting..
            self.edge_initialization_flag = 1

            # Creating INTRA-COMMUNITY, create community structure
            for community in self.communities:
                # Getting nodes that belong to that community
                nodes = [node for node in self.graph.nodes if self.graph.nodes[node]['community'] == community]
                
                # Creating potential set of edges
                potential_edges = list(combinations(nodes,r=2))
                
                # Establishing the right community probabilities
                intra_community_p = c_intra / (len(nodes) - 1)

                # Assigning with c_intra probability nodes within community
                # and saving these edges for later sampling for disease propagation

                self.edges = [edge for edge in potential_edges if np.random.uniform(0,1,1) < intra_community_p]


            if self.no_communities > 1:
                inter_community_connections = list(combinations(self.communities,r=2))
                
                # Creating INTER-COMMUNITY, connect nodes with probability c_inter
                for community_pair in inter_community_connections:
                    nodes_1 = [node for node in self.graph.nodes if self.graph.nodes[node]['community'] == community_pair[0]]
                    nodes_2 = [node for node in self.graph.nodes if self.graph.nodes[node]['community'] == community_pair[1]]
                    
                    # Getting potential edges
                    potential_inter_edges = list(itertools.product(nodes_1,nodes_2))
                    
                    # Assigning edges with probability c_inter.
                    self.edges += [edge for edge in potential_inter_edges if np.random.uniform(0,1,1) < self.p_inter]

        
    def advance_disease(self,  percolation : float):
        """
        This function advances the disease state by 1 time period. 
        In its current form, the function:
        
        1) Samples a subset of these edges, and activates them.
        2) Propagates disease along these edges if they are in contact, 
           w.p. self.p_infection
        """
        self.percolation = percolation
        # Start by clearing all the edges in the graph from the grapher function.
        # This is done because of the graphing structure.
        self.graph.remove_edges_from(list(self.graph.edges()))
        
            
        # Sampling from our edges that we have, with probability percolation
        self.graph.add_edges_from(random.sample(list(self.graph.edges),int(self.percolation * len(self.graph.edges))))
        
        # Identifying what nodes are already infected
        infected_nodes = [node for node in self.graph.nodes if self.graph.nodes[node]['infected'] == 1]
        
        # Advancing the old infected nodes by one time period
        for node in infected_nodes:
            self.graph.nodes[node]['time_since_infection'] += 1
            
            # Passing to recovery period if some have been infected
            if self.graph.nodes[node]['time_since_infection'] > self.infection_period:
                # Current state
                self.graph.nodes[node]['susceptible'] = 0
                self.graph.nodes[node]['infected'] = 0
                self.graph.nodes[node]['recovered'] = 1
                
                # Stopping infection and starting recovered timers
                self.graph.nodes[node]['time_since_infection'] = 0
                self.graph.nodes[node]['time_since_recovered'] += 1
                
         # Advancing the old recovered nodes by one time period
        recovered_nodes = [node for node in self.graph.nodes if self.graph.nodes[node]['recovered'] == 1]
        for node in recovered_nodes:
            self.graph.nodes[node]['time_since_recovered'] += 1
            
            # Passing to susceptibility those that have exceeded self.recovery_period
            if self.graph.nodes[node]['time_since_recovered'] > self.recovery_period:
               
                # Current state
                self.graph.nodes[node]['susceptible'] = 1
                self.graph.nodes[node]['infected'] = 0
                self.graph.nodes[node]['recovered'] = 0
                
                # Stopping recovered timer
                self.graph.nodes[node]['time_since_recovered'] = 0

        
        # Making sure that we discount edges that are already recovered.
        infected_nodes = [node for node in self.graph.nodes if self.graph.nodes[node]['infected'] == 1]
        
        # Filtering out the new infected nodes
        susceptible_nodes = [node for node in self.graph.nodes if (self.graph.nodes[node]['susceptible'] == 1)\
                                                              and (self.graph.nodes[node]['recovered'] == 0)]
        
        # Infecting edges, remember that we have the probability of infection to account for too!
        infecting_edges = [edge for edge in self.graph.edges \
                           if (edge[0] in infected_nodes and edge[1] in susceptible_nodes) \
                           or (edge[1] in infected_nodes and edge[0] in susceptible_nodes) \
                           and (np.random.uniform(0,1,1) < self.p_infection)
                                                                  ]

        # set infecting node color to different color, and pass this
        new_infected_nodes = []
        for edge in infecting_edges:
            for node in edge:
                if (node not in infected_nodes) & (self.graph.nodes[node]['susceptible'] == 1):
                    new_infected_nodes.append(node)
            
        # Infecting the new nodes     
        for node in new_infected_nodes:
            # Current state
            self.graph.nodes[node]['infected'] = 1
            self.graph.nodes[node]['susceptible'] = 0
            self.graph.nodes[node]['recovered'] = 0
            
            # Starting infection timer
            self.graph.nodes[node]['time_since_infection'] += 1

        # Updating node status history

        self.nodehistory['susceptible'].append(sum(nx.get_node_attributes(self.graph,'susceptible').values())/len(self.graph.nodes))
        self.nodehistory['infected'].append(sum(nx.get_node_attributes(self.graph,'infected').values())/len(self.graph.nodes))
        self.nodehistory['recovered'].append(sum(nx.get_node_attributes(self.graph,'recovered').values())/len(self.graph.nodes))        


        # for node in self.graph.nodes:
        #     self.nodehistory[node]['susceptible'].append(self.graph.nodes[node]['susceptible'])
        #     self.nodehistory[node]['infected'].append(self.graph.nodes[node]['infected'])
        #     self.nodehistory[node]['recovered'].append(self.graph.nodes[node]['recovered'])

        self.timeperiod += 1  
     
        return infecting_edges
    
    def get_fixed_positions(self):
        """
        This function takes the nodes of a graph, and creates a fixed layout so that
        it is easier to see how the disease is propagating in a collection of 
        communities in a village!
        """
        
        
        connections_list = []
        for community in self.communities:
            nodes_concerned = [node for node in self.graph.nodes if self.graph.nodes[node]['community'] == community]
            all_community_edges = list(combinations(nodes_concerned,r=2))
            connections = {edge : 10 for edge in all_community_edges}
            connections_list.append(connections)
        
        dict_of_connections = {}
        for con in connections_list:
            dict_of_connections.update(con)
        
        
        remaining_edges = [edge for edge in list(combinations(self.graph.nodes,r=2)) if edge not in dict_of_connections.keys()]
        non_connections = {edge : 0 for edge in remaining_edges}
        dict_of_connections.update(non_connections)
        
        self.graph.add_edges_from(dict_of_connections.keys())

        nx.set_edge_attributes(self.graph, dict_of_connections, "weight")
                
        return None

    
    def grapher(self,infecting_edges = 0, iteration = -1):
        
        """
        Takes a graph, and plots the current state of infection of the graph.
        """
                
        partition = 'infected'

        if self.posfix == 0:
            
            self.pos =  nx.spring_layout(self.graph, k = 1/np.sqrt(0.3))
            self.posfix = 1
            
        
        # Visualizing the infected nodes
        groups = set(nx.get_node_attributes(self.graph,partition).values())
        mapping = dict(zip(sorted(groups),count()))
        nodes = self.graph.nodes()
        colors = [mapping[self.graph.nodes[n][partition]] for n in nodes]
        
        color_map = {}
        susceptible_nodes = []
        infected_nodes = []
        recovered_nodes = []
        for node in self.graph:
            if  self.graph.nodes[node]['infected'] == 1:
                color_map[node] = 'red'
                infected_nodes.append(node)
            elif self.graph.nodes[node]['recovered'] == 1: 
                color_map[node] = 'green'
                recovered_nodes.append(node)
            elif self.graph.nodes[node]['susceptible'] == 1:
                color_map[node] = 'orange'
                susceptible_nodes.append(node)
        
#         # Initializing the figure
#         f = plt.figure()
#         f.set_figheight(15)
#         f.set_figwidth(15)
            
#         # Initializing the plots
#         ax1 = f.add_subplot(121)
#         ax2 = f.add_subplot(122)
        
        
        fig, ax = plt.subplots(1,1, figsize=(15,15))
        no_infected = sum(nx.get_node_attributes(self.graph,'infected').values())
        fig.suptitle(f"Village contagion, time-period {self.timeperiod}, {no_infected} infected villagers.", fontsize = 30)
        text = f"p-intra: {round(self.p_intra,6)}, p-inter: {round(self.p_inter,6)}, \n sociality: {self.sociality}, p_infection: {self.p_infection}"
        fig.text(.5, .05, text, ha='center', FontSize = 20)

        if infecting_edges == 0:
            nx.draw_networkx(self.graph,
                        node_size = 40,
                        with_labels = False,
                        width = 0.0,
                        node_color = color_map.values(),
                        pos = self.pos,
                        ax = ax)
                   
            
#             nx.draw_networkx_nodes(self.graph,
#                                    pos=self.pos,
#                                    nodelist=susceptible_nodes,
#                                    node_color='orange',
#                                    label='Susceptible',
#                                    node_size = 50,
#                                    ax = ax)

#             nx.draw_networkx_nodes(self.graph,
#                                    pos=self.pos,
#                                    nodelist=infected_nodes,
#                                    node_color='red',
#                                    label='Infected',
#                                    node_size = 50,
#                                    ax = ax)

#             nx.draw_networkx_nodes(self.graph,
#                                    pos=self.pos,
#                                    nodelist=recovered_nodes,
#                                    node_color='green',
#                                    label='Recovered',
#                                    node_size = 50,
#                                    ax = ax)

        else:
            nx.draw_networkx(self.graph,
                        node_size = 40,
                        with_labels = False,
                        width = 1,
                        node_color = color_map.values(),
                        pos = self.pos,
                        edgelist = infecting_edges,
                        edge_color = 'red',
                        ax = ax)
            
        if iteration != -1:
            plt.savefig(os.path.join(os.getcwd(),f"for_gif/{iteration}.png"))
            # plt.show()
        plt.close(fig)


def to_gif():

    images = []
    path = os.path.join(os.getcwd(),'for_gif/')
    filenames = os.listdir(path)
    print(filenames)
    filenames.sort(key = lambda x: int(x.split(".")[0]))

    for filename in filenames:
        newname = filename.split(".")[0]+".JPEG"
        im = Image.open(path+filename)
        rgb_im = im.convert('RGB')
        rgb_im.save(path+newname)

    filenames = os.listdir(path)
    filenames = [filename for filename in filenames if 'JPEG' in filename]
    filenames.sort(key = lambda x: int(x.split(".")[0]))

    for filename in filenames:
        images.append(imageio.imread(path+filename))
        imageio.mimsave(os.path.join(os.getcwd(), 'movie.gif'), images, duration = 3)



##############################################
# ALL FUNCTIONS BELOW ARE NOT CLASS-SPECIFIC!#
##############################################


def get_model1(n):
    """
    Returns a model with parameters:
    no_communities = 1
    c_intra = 0.05 * n (not used because no_communities = 1)
    c_inter = 0.05 * n

    """
    model  = MJDM(n)
    model.assign_community(no_communities = 1, c_intra = 0.05, c_inter = 0.05*n)
    return model

def get_model2(n):
    """
    Returns a model with parameters:
    no_communities = 1
    c_intra = 0.05 * n (not used because no_communities = 1)
    c_inter = 0.05 * n

    """
    model  = MJDM(n)
    model.assign_community(no_communities = 1, c_intra = 0.05, c_inter = 0.05*n)
    return model  

def get_model3(n):
    """
    Returns a model with parameters:
    no_communities = 1
    c_intra = 0.05 * n (not used because no_communities = 1)
    c_inter = 0.05 * n

    """
    model  = MJDM(n)
    model.assign_community(no_communities = 1, c_intra = 0.2, c_inter = 0.025*n)
    return model

def get_model4(n):
    """
    Returns a model with parameters:
    no_communities = 1
    c_intra = 0.05 * n (not used because no_communities = 1)
    c_inter = 0.05 * n

    """
    model  = MJDM(n)
    model.assign_community(no_communities = 1, c_intra = 0.2, c_inter = 0.1*n)
    return model


def simulation(parameter_combinations,no_simulations,time_periods):
    """
    This function takes in all combinations of parameters.

    It then runs simulations of all of these parameters.

    In the process, it records everything into a dict,
    which can be saved as a dataframe.

    Input (in parameter_combinations):
    n_nodes
    p_infections,
    n_infecteds,
    infection_periods,
    recovery_periods,
    percolation,

    """

    data_dict = {
        'model_number' : [],
        'parameter_combination' : [],
        'simulation_number' : [],
        'time_period' : [],
        'susceptible_pct' : [],
        'infected_pct' : [],
        'recovered_pct' : []
    }

    combinations = len(parameter_combinations)

    for idx,comb in enumerate(parameter_combinations):
        # Run a simulation for a certain parameter combination across all
        # four models

        for run in range(no_simulations):
            
            # Getting our models with their community structure
            model1 = get_model1(comb[0])
            model2 = get_model2(comb[0])
            model3 = get_model3(comb[0])
            model4 = get_model4(comb[0])
            
            # Assigning models to list 
            model_list = [model1,model2,model3,model4]

            # Infecting each model
            for model_id,model in enumerate(model_list):
                model.infect(p_infection = comb[0+1],\
                             no_infected = comb[1+1],\
                             infection_period = comb[2+1],\
                             recovery_period = comb[3+1])           
            
                for t in range(time_periods):
                    infecting_edges = model.advance_disease(percolation = comb[4+1])
                    if (len(infecting_edges) == 0) & (model.nodehistory['infected'][-1] == 0):
                        break
            
                tps = len(model.nodehistory['infected'])
                data_dict['model_number'].extend([model_id] * tps)
                data_dict['parameter_combination'].extend([comb] * tps)
                data_dict['simulation_number'].extend([run] *tps)
                data_dict['time_period'] += list(range(tps))
                data_dict['susceptible_pct'] += model.nodehistory['susceptible']
                data_dict['infected_pct'] += model.nodehistory['infected']
                data_dict['recovered_pct'] += model.nodehistory['recovered']     


        print(f"Finished round number {idx+1} out of {combinations}!")
    # Saving the data as a dataframe
    # df = pd.DataFrame(data_dict)

    return data_dict
