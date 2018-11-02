import networkx as nx
import numpy as np
import matplotlib.pyplot as plt 
from networkx.drawing.nx_agraph import graphviz_layout


def hierarchy_pos(G, root, width=10., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5, 
                  pos = None, parent = None):
    '''If there is a cycle that is reachable from root, then this will see infinite recursion.
       G: the graph
       root: the root node of current branch
       width: horizontal space allocated for this branch - avoids overlap with other branches
       vert_gap: gap between levels of hierarchy
       vert_loc: vertical location of root
       xcenter: horizontal location of root
       pos: a dict saying where all nodes go if they have been assigned
       parent: parent of this branch.'''
    if pos == None:
        pos = {root:(xcenter,vert_loc)}
    else:
        pos[root] = (xcenter, vert_loc)
    neighbors = G.neighbors(root)

    if len(neighbors)!=0:
        dx = width/len(neighbors) 
        nextx = xcenter - width/2 - dx/2
        for neighbor in neighbors:
            nextx += dx
            pos = hierarchy_pos(G,neighbor, width = dx, vert_gap = vert_gap, 
                                vert_loc = vert_loc-vert_gap, xcenter=nextx, pos=pos, 
                                parent = root)
    return pos
def grow_y_network_repeated(p_1arm, p_2arm, p_3arm, n_repeats = 1000):
	#repeatedly growing y networks for a given set of probabilities and aggregating statistics
	#on average number of nodes and likelihood of growth being terminated
	terminated_list = []
	n_nodes_list = []
	for i in range(n_repeats):
		G, terminated = grow_y_network(p_1arm, p_2arm, p_3arm)
		terminated_list.append(terminated)
		n_nodes_list.append(G.number_of_nodes())
	#print n_nodes_list
	average_nodes = sum(n_nodes_list) / float(len(n_nodes_list))
	termination_probability = float(terminated_list.count(True))/float(len(terminated_list))

	return average_nodes, termination_probability


def grow_y_network(p_1arm, p_2arm, p_3arm, n_iterations = 10):
	#writing a script to grow a nanotube network, the inputs will be the relative
	#abundance of 0,1,2, and 3 armed structures (nodes)
	#the graph will be allowed to grow until there are no suitable points to add additional edges
	#then the goal will be to analyze the size of the network (number of added nodes) as a function of the relative 
	#abundances of each species, this could definitely be extended to include L and T junctions if needed! 

	###############
	#hard coded relative abundances
	#p_3arm = .50
	#p_2arm = .30
	#p_1arm = .20

	connectivities = [1, 2, 3]
	probabilities = [p_1arm, p_2arm, p_3arm]

	#initializing a directed graph
	G = nx.DiGraph()

	#initializing the first node
	first_node_connectivity = np.random.choice(connectivities, p = probabilities)
	G.add_node(1, arms = first_node_connectivity)

	#adding the first shell of connecting nodes
	for i in range(first_node_connectivity):
		n_arms = np.random.choice(connectivities, p = probabilities)

		if n_arms == 1:
			G.add_node(G.number_of_nodes()+1, arms = 1, terminal = True)
			G.add_edge(1, G.number_of_nodes())
		else:
			G.add_node(G.number_of_nodes()+1, arms = n_arms)
			G.add_edge(1, G.number_of_nodes())

	#plt.subplot(121)
	#nx.draw(G, with_labels = True)
	#plt.show()

	#now looping through many iterations, each iteration:
	# 1. finds the nodes that were added during the last iterations
	# 2. based on the number of arms for each node that was added in the last iteration, add the appropriate
	#    number of nodes or terminate growth from a particular node 
	# 2.5. if no nodes allow for the addition of new nodes, terminate the growth process
	# 3. Add the appropriate edges to each new node
	terminated = False
	for i in range(n_iterations):
		#print 'growth iterataion: ', i 
		#need to loop over all non-terminal nodes added during the last iteration
		#in other words, loop over all nodes that do not have children and are not terminal

		#grabbing a list of nodes that do not have childern and are not terminals
		terminal_nodes_dict = nx.get_node_attributes(G, 'terminal')
		connectivity_dict = nx.get_node_attributes(G, 'arms')
		terminal_nodes = terminal_nodes_dict.keys()
		growth_nodes = [node for node, out_degree in G.out_degree_iter() if out_degree == 0 and node not in terminal_nodes]

		#if there are no growth nodes, the network is terminated
		if len(growth_nodes) == 0:
			#print 'network growth terminated!'
			terminated = True
			break
		for growth_node in growth_nodes:
			connectivity = connectivity_dict[growth_node]
			for arm in range(connectivity-1):
				n_arms = np.random.choice(connectivities, p = probabilities)

				if n_arms == 1:
					G.add_node(G.number_of_nodes()+1, arms = 1, terminal = True)
					G.add_edge(growth_node, G.number_of_nodes())
				else:
					G.add_node(G.number_of_nodes()+1, arms = n_arms)
					G.add_edge(growth_node, G.number_of_nodes())

	return G, terminated


#ok, idea is to use gillespie-type algorithm to explicitly track all species in our system
#instead of designating "species" this algorithm will explicitly account for all structures in the system
#all possible joining events will be considered and their rates calculated, then the gillespie algorithm can be applied

#system will be one giant graph and connected_component_subgraphs will be used to identify individual species

#all nodes will contain "arms" data on whether they are a 1, 2, or 3 arm Y junction. We can detect tube ends by finding nodes that have more arms
#than they do connections to other nodes

#network joining will occur by adding an edge between two tube ends (that is, an edge between two nodes that both have arms>connections)


def initialize_system():
	#using expt data we will initialize the system with the proper number of 1, 2, 3 armed structures
	system = nx.Graph()
	n_1arm_pos = 7
	n_2arm_pos = 13
	n_3arm_pos = 13

	n_1arm_neg = 13
	n_2arm_neg = 27
	n_3arm_neg = 27
	for i in range(n_1arm_pos):
		system.add_node(system.number_of_nodes()+1, arms = 1, adapter = 'pos')
	for i in range(n_2arm_pos):
		system.add_node(system.number_of_nodes()+1, arms = 2, adapter = 'pos')
	for i in range(n_3arm_pos):
		system.add_node(system.number_of_nodes()+1, arms = 3, adapter = 'pos')

	for i in range(n_1arm_neg):
		system.add_node(system.number_of_nodes()+1, arms = 1, adapter = 'neg')
	for i in range(n_2arm_neg):
		system.add_node(system.number_of_nodes()+1, arms = 2, adapter = 'neg')
	for i in range(n_3arm_neg):
		system.add_node(system.number_of_nodes()+1, arms = 3, adapter = 'neg')

	print system.number_of_nodes()

	return system

def compute_rtot(system):
	#given a system, compute the total reaction rate for the gillespie algorithm
	#this means consider all possible joining reactions, calculate the rate for each, and sum this 

	#first, we need a list of all network species: all connected subgraphs
	graphs = list(nx.connected_component_subgraphs(system))
	print len(graphs)

	rtot = 0
	#now we consider all pairwise joining reactions for these subgraphs, calculate the joining rate, and add it to the rtot
	for i in range(len(graphs)):
		for j in range (len(graphs)):
			if i == j:
				continue
			else:
				network_1 = graphs[i]
				network_2 = graphs[j]
				kjoin = joining_rate(network_1, network_2)
				rtot += kjoin
	print rtot 
	return rtot 





#we will need a method that considers two possible networks that can join and computes their joining rate
#this will either be a constant (for the constant model) or it will scale with the number of free ends exposed on the network

def joining_rate(network_1, network_2):
	#placeholder, simplest joining rate is a constant rate
	#note that even for the constant model this should return zero if the two networks cannot be joined
	#that is if either one does not have any nodes with arms>connections
	#

	#rather than compute the joining rate for all pairs every time we could try updating rtot and all the reaction probabilities 
	connectivity_dict1 = nx.get_node_attributes(network_1, 'arms')
	adapter_dict1 = nx.get_node_attributes(network_1, 'adapter')
	connectivity_dict2 = nx.get_node_attributes(network_2, 'arms')
	adapter_dict2 = nx.get_node_attributes(network_2, 'adapter')

	#count valid attachment points for network 1
	valid_positive_attachment_points_1 = 0
	valid_negative_attachment_points_1 = 0
	#has_valid_attachment_point_1 = False
	for node, narms in connectivity_dict1.iteritems():
		if network_1.degree(node) < narms and network_1.node[node]['adapter'] == 'pos':
			valid_positive_attachment_points_1 += narms - network_1.degree(node)
		if network_1.degree(node) < narms and network_1.node[node]['adapter'] == 'neg':
			valid_negative_attachment_points_1 += narms - network_1.degree(node)

	#count valid attachment points for network 2
	valid_positive_attachment_points_2 = 0
	valid_negative_attachment_points_2 = 0
	#has_valid_attachment_point_1 = False
	for node, narms in connectivity_dict2.iteritems():
		if network_2.degree(node) < narms and network_2.node[node]['adapter'] == 'pos':
			valid_positive_attachment_points_2 += narms - network_2.degree(node)
		if network_2.degree(node) < narms and network_2.node[node]['adapter'] == 'neg':
			valid_negative_attachment_points_2 += narms - network_2.degree(node)
			

	number_of_possible_connections = float(valid_positive_attachment_points_1*valid_negative_attachment_points_2 + valid_positive_attachment_points_1*valid_negative_attachment_points_2)
	kjoin = .005
	return kjoin*number_of_possible_connections

def compute_dt(rtot):
	#compute dt from the total reaction rate according to the exponential distribution
	dt = np.random.exponential(float(1.0/rtot))
	return dt

def choose_joining_reaction(rtot, system):
	#randomly choose a joining reaction with probabilities equal to the reaction rate over rtot
	#first, we need a list of all network species: all connected subgraphs
	graphs = list(nx.connected_component_subgraphs(system))
	print "old number of connected graphs", len(graphs)


	#now we consider all pairwise joining reactions for these subgraphs, calculate the joining rate, and divide by rtot to get the probability 
	joining_probabilities = []
	network_pairs = []
	for i in range(len(graphs)):
		for j in range (len(graphs)):
			if i == j:
				network_1 = graphs[i]
				network_2 = graphs[j]
				kjoin = 0 #a network will never join to itself in this model
			else:
				network_1 = graphs[i]
				network_2 = graphs[j]
				kjoin = joining_rate(network_1, network_2)
			joining_probability = float(kjoin)/rtot
			#print joining_probability
			joining_probabilities.append(joining_probability)
			network_pairs.append([network_1, network_2])
	network_pairs_index = range(len(network_pairs))

	chosen_pair_index = np.random.choice(network_pairs_index, p = joining_probabilities)
	print "chosen pair index: ", chosen_pair_index
	chosen_pair = network_pairs[chosen_pair_index]

	return chosen_pair[0], chosen_pair[1]

def perform_joining(network_1, network_2, system):
	#joing network_1 and network_2 together by adding an edge between them
	#this edge is added at a random available "tube endpoint"
	connectivity_dict1 = nx.get_node_attributes(network_1, 'arms')
	connectivity_dict2 = nx.get_node_attributes(network_2, 'arms')
	valid_connection = False
	while valid_connection == False: 
		network_1_node = np.random.choice(network_1.nodes())
		network_2_node = np.random.choice(network_2.nodes())
		network_1_node_narms = connectivity_dict1[network_1_node]
		network_2_node_narms = connectivity_dict2[network_2_node]
		if network_1.degree(network_1_node) < network_1_node_narms and network_2.degree(network_2_node) < network_2_node_narms and network_1.node[network_1_node]['adapter'] != network_2.node[network_2_node]['adapter'] :
			print "adding valid connection"
			valid_connection = True

	system.add_edge(network_1_node, network_2_node)

	return system



'''
p_1arms = np.arange(0.1, 0.9, 0.1)
p_23_splits = np.arange(0.0,1.1,0.1)

for p_1arm in p_1arms:
	average_nodes_list = []
	termination_probability_list = []
	p_3arm_list = []
	for p_23_split in p_23_splits:
		p_3arm = p_23_split * (1.0 - p_1arm)
		p_2arm = 1.0 - p_1arm - p_3arm

		average_nodes, termination_probability = grow_y_network_repeated(p_1arm, p_2arm, p_3arm)
		print p_3arm, average_nodes, termination_probability 

		average_nodes_list.append(average_nodes)
		termination_probability_list.append(termination_probability)
		p_3arm_list.append(p_3arm)

	plt.plot(p_3arm_list, average_nodes_list, label = 'p_1arm: '+"{:4.2f}".format(p_1arm), linewidth = 5)

plt.xlabel('probability of 3-arm structure')
plt.ylabel('total number of nodes after 10 iterations')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 50.0])
plt.legend()
plt.savefig('number_of_nodes_vs_p_3arm.pdf')

'''

'''plt.xlabel('probability of 3-arm structure')
plt.ylabel('probability of network termination before 10 iterations')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.legend()
plt.savefig('termination_probability_vs_p_3arm.pdf')'''

'''pos = hierarchy_pos(G,1)
#plt.subplot(121)
nx.draw(G, pos = pos, with_labels = True)
plt.savefig("network_"+str(p_1arm)+"_"+str(p_2arm)+"_"+str(p_3arm)+".png")
'''
def perform_gillespie_simulation(i):
	#initialize the system according the expt starting conditions
	system = initialize_system()
	number_of_joining_steps = system.number_of_nodes()-10 #we will stop just short of everything in the system being joined to itself
	average_nodes_per_network = [1]
	number_of_connected_graphs = [50]
	average_seeds_per_network = [1]
	average_network_size_of_a_seed_list = [1]
	time = [0]
	for step in range(number_of_joining_steps):
		rtot = compute_rtot(system)
		dt = compute_dt(rtot)
		print dt 
		previous_time = time[len(time)-1]
		time.append(previous_time + dt)
	
		#now we need to determine which joining reaction will occur
		#the probability for a given reaction is that reaction's rate over rtot
	
		joined_network_1, joined_network_2 = choose_joining_reaction(rtot, system)
	
	
		connectivity_dict1 = nx.get_node_attributes(joined_network_1, 'arms')
		connectivity_dict2 = nx.get_node_attributes(joined_network_2, 'arms')
		print connectivity_dict1, connectivity_dict2
	
		#now we actually perform the joining reaction by creating an edge between the two networks  
		system = perform_joining(joined_network_1, joined_network_2, system)
	
		graphs = list(nx.connected_component_subgraphs(system))
		print "new number of connected graphs", len(graphs)
		number_of_connected_graphs.append(len(graphs))
		average_seeds_per_network.append(50.0/float(len(graphs)))

		#now we will compute the average network size that a seed is in
		network_size_count = []
		for graph in graphs:
			for node in graph.nodes():
				network_size_count.append(graph.number_of_nodes())
		average_network_size_of_a_seed = float(sum(network_size_count))/float(len(network_size_count))

		average_network_size_of_a_seed_list.append(average_network_size_of_a_seed)

	print time
	print number_of_connected_graphs
	print average_seeds_per_network
	print average_network_size_of_a_seed_list



	plt.plot(time, average_network_size_of_a_seed_list, label = 'Gillespie algorithm run '+str(i+1), linewidth = 5, alpha = .5)
	

for i in range(1):
	perform_gillespie_simulation(i)

expt_times = [0.5, 2.5, 4.5, 6.5, 8.5]
expt_seeds_per_network = [2.45, 4.68, 6.84, 19.33, 17.32]
expt_errors = [.045, .204, .280, 1.12, 1.14]

#expt_average_network_size_of_seed = [1.51, 2.42, 3.61, 4.34]

plt.errorbar(expt_times, expt_seeds_per_network, expt_errors, label = 'experiment', linewidth = 5, color = 'black' )
plt.xlim([0.0, 10])
plt.ylim([0.0, 30])
plt.legend(loc = 4, fontsize = 8)
plt.xlabel('time (hours)')
plt.ylabel('average network size that a seed is in')
plt.savefig("gillespie_results.pdf")











	

