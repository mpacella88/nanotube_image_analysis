import networkx as nx
import numpy as np
import matplotlib.pyplot as plt 
#from networkx.drawing.nx_agraph import graphviz_layout
#from datetime import datetime
import pdb
#startTime = datetime.now()

#optimizing performance 


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
	n_1arm_pos = 150 #7
	n_2arm_pos = 400 #13
	n_3arm_pos = 550 #13

	n_1arm_neg = 150 #13
	n_2arm_neg = 400 #27
	n_3arm_neg = 550#27
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

	print "initial number of nodes: ", system.number_of_nodes()

	return system

def compute_rtot_first(system):
	#given a system, compute the total reaction rate for the gillespie algorithm
	#this means consider all possible joining reactions, calculate the rate for each, and sum this 

	#first, we need a list of all network species: all connected subgraphs
	graphs = list(nx.connected_component_subgraphs(system))
	#print len(graphs)

	rtot = 0
	#now we consider all pairwise joining reactions for these subgraphs, calculate the joining rate, and add it to the rtot
	for i in range(len(graphs)):
		for j in range (i+1, len(graphs)):
			if i == j:
				continue
			else:
				network_1 = graphs[i]
				network_2 = graphs[j]
				kjoin = joining_rate(network_1, network_2)
				print kjoin 
				rtot += kjoin
	#print rtot 
	return rtot 


def update_rtot(old_rtot, network_1, network_2, new_network, system_old):
	#given a joining reaction, update rtot by subtracting all terms involving either of the joining partners 
	#and add terms involving joining to the new network

	
	#we need a list of all subgraphs not including network_1, network_2, or new_network
	graphs = list(nx.connected_component_subgraphs(system_old))
	print len(graphs)

	rtot = old_rtot
	
	for i in range(len(graphs)):
		#subtract all joining terms for network_1
		network = graphs[i]
		kjoin1 = joining_rate(network_1, network)
		rtot -= kjoin1
		print "subtracting: ", kjoin1
		
		#subtract all joining terms for network_2
		kjoin2 = joining_rate(network_2, network)
		rtot -= kjoin2
		print "subtracting: ", kjoin2

		#add joining terms for new_network
		kjoin_new = joining_rate(new_network, network)
		rtot+=kjoin_new 
		print "adding: ", kjoin_new 
		
	network_1_correction = joining_rate(network_1, network_1) #this is to avoid double counting the interaction of the two networks with themselves 
	network_2_correction = joining_rate(network_2, network_2) #this is to avoid double counting the interaction of the two networks with themselves
	rtot += network_1_correction
	rtot += network_2_correction
	#print "adding network 1 2 correction: ", network_1_2_correction
	new_network_correction_1 = joining_rate(new_network, network_1) #this is to avoid counting the interaction of the new network with the two old networks 
	new_network_correction_2 = joining_rate(new_network, network_2)
	rtot -= new_network_correction_1
	rtot -= new_network_correction_2


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
	#when a joining event occurs rtot is reduced because there is one less possible reaction
	#the joining events corresponding to joining with either of the partners are eliminated
	#and a new set of joining events for joining to the new network are added  
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
			
	#print "valid pos 1: ", valid_positive_attachment_points_1
	#print "valid neg 1: ", valid_negative_attachment_points_1
	#print "valid pos 2: ", valid_positive_attachment_points_2
	#print "valid neg 2: ", valid_negative_attachment_points_2
	number_of_possible_connections = float(valid_positive_attachment_points_1*valid_negative_attachment_points_2 + valid_negative_attachment_points_1*valid_positive_attachment_points_2)
	#print "number of possible connections: ", number_of_possible_connections
	kjoin = 1.0/(float(network_1.number_of_nodes()**3)*float(network_2.number_of_nodes()**3))
	return kjoin*number_of_possible_connections

def compute_dt(rtot):
	#compute dt from the total reaction rate according to the exponential distribution
	dt = np.random.exponential(float(1.0/rtot))
	return dt

def choose_joining_reaction(rtot, system, joining_matrix, joining_matrix_index_dictionary):
	#we also need to consider intra-network joining reactions here.......

	#randomly choose a joining reaction with probabilities equal to the reaction rate over rtot
	#first, we need a list of all network species: all connected subgraphs
	graphs = list(nx.connected_component_subgraphs(system))
	print "old number of connected graphs", len(graphs)

	#now we will use the joining matrix to determine the probabilities for all pairwise joining reactions
	pairwise_probabilities = joining_matrix.ravel()/np.sum(joining_matrix)
	ravel_indices = []
	for i in range(len(pairwise_probabilities)):
		ravel_indices.append(i)

	chosen_ravel_index = np.random.choice(ravel_indices, p = pairwise_probabilities)
	chosen_matrix_indices = np.unravel_index(chosen_ravel_index, joining_matrix.shape)

	#now that we have the joining matrix indices we will use the dictionary to convert those to actual networks
	print chosen_matrix_indices[0]
	index_of_network_1 = chosen_matrix_indices[0]
	index_of_network_2 = chosen_matrix_indices[1]
	lowest_node_in_network_1 = joining_matrix_index_dictionary[index_of_network_1]
	lowest_node_in_network_2 = joining_matrix_index_dictionary[index_of_network_2]

	for network in graphs:
		if network.has_node(lowest_node_in_network_1):
			network_1 = network
			break
	for network in graphs:
		if network.has_node(lowest_node_in_network_2):
			network_2 = network
			break


	print "lowest node in network 1:", lowest_node_in_network_1
	print "lowest node in network 2:", lowest_node_in_network_2
	
	#pdb.set_trace()
	return network_1, network_2, index_of_network_1, index_of_network_2

def perform_joining(network_1, network_2, system):
	#joing network_1 and network_2 together by adding an edge between them
	#this edge is added at a random available "tube endpoint"
	print network_1
	print network_2

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
		#pdb.set_trace()
	system.add_edge(network_1_node, network_2_node)
	new_network = nx.union(network_1, network_2)
	new_network.add_edge(network_1_node, network_2_node)
	

	return system, new_network



def initialize_joining_matrix(system):
	#initializing the 2D joining matrix for the first time
	#subsequently, this joining matrix can simply be updated by using a separate routine that does not involve large loops
	#joining matrix only handles inter-network joining, a separate routine will be needed to calculate intra-network joining

	#initialize matrix, it will be square with dimensions equal to the number of nodes in the initial system
	joining_matrix = np.zeros((system.number_of_nodes(),system.number_of_nodes()))
	graphs = list(nx.connected_component_subgraphs(system))
	#now we will populate each matrix element by calculating the pairwise joining rates for all nodes in the system
	#this will be a symmetric matrix because joining is symmetric, the diagonal elements will be zero because a node cannot join to itself
	#because it is symmetric we will only calculate/store the elements in the upper diagonal

	#we will also create a dictionary that associates an index in the joining_matrix with a network in the system graph, this will need to be updated
	#the dictionary key will be an integer: the lowest valued node number contained in the network
	#for this initialization step the key is simpley the node number
	joining_matrix_index_dictionary = {}


	for i in range(len(graphs)):
		lowest_node_number = min(graphs[i].nodes())
		joining_matrix_index_dictionary[i] = lowest_node_number 
		for j in range (i+1, len(graphs)):
			network_1 = graphs[i]
			network_2 = graphs[j]
			kjoin = joining_rate(network_1, network_2)
			joining_matrix[i,j] = kjoin 
			#rtot += kjoin
	return joining_matrix, joining_matrix_index_dictionary

def compute_rtot_from_joining_matrix(joining_matrix):
	#computing rtot from the initial joining matrix
	#sum all elements in the matrix
	matrix_sum = np.sum(joining_matrix)
	return matrix_sum

def update_joining_matrix(joining_matrix, joining_matrix_index_dictionary, joined_network_1, joined_network_2, new_network, index_1, index_2, system):
	#this will update the joining matrix after two networks have been joined and will also update the joining matrix dictionary
	#first things first, find the lowest node number for each of the joined networks, these are used to find the correct elements of
	#the joining matrix to edit
	lowest_node_1 = min(joined_network_1.nodes())
	lowest_node_2 = min(joined_network_2.nodes())

	#now we want to remove the row and column associated with network 1
	joining_matrix = np.delete(joining_matrix, index_1, 0) #deletes the row for network 1
	joining_matrix = np.delete(joining_matrix, index_1, 1) #deletes the column for network 1

	#now, the indices of all networks with an index larger than index 1 will have decreased by one
	#we will update the joining matrix index dictionary accordingly 
	updated_dict_1 = {}
	for key, value in joining_matrix_index_dictionary.iteritems():
		if key< index_1:
			#this means this entry is unaffected by the deletion
			updated_dict_1[key] = value 
		if key> index_1:
			#this means that the key should be reduced by 1 because of the deletion
			new_key = key - 1
			updated_dict_1[new_key] = value
		#this dictionary is now one key smaller than the origianl b/c key==index 1 was not covered

	#do the same thing for index 2 but using the updated joining matrix and dict

	#now we want to remove the row and column associated with network 1
	print index_2
	#IF INDEX 2 WAS GREATER THAN INDEX 1 IT NEEDS TO BE SHIFTED BEFORE THE CUT IS MADE!!!
	if index_2>index_1:
		index_2-=1
	joining_matrix = np.delete(joining_matrix, index_2, 0) #deletes the row for network 2
	joining_matrix = np.delete(joining_matrix, index_2, 1) #deletes the column for network 2

	#now, the indices of all networks with an index larger than index 1 will have decreased by one
	#we will update the joining matrix index dictionary accordingly 
	updated_dict_2 = {}
	for key, value in updated_dict_1.iteritems():
		if key< index_2:
			#this means this entry is unaffected by the deletion
			updated_dict_2[key] = value 
		if key> index_2:
			#this means that the key should be reduced by 1 because of the deletion
			new_key = key - 1
			updated_dict_2[new_key] = value
		#this dictionary is now one key smaller than the origianl b/c key==index 2 was not covered

	#now we will add a row and column for the new network that was created by the joining event!
	#this will be the rightmost column and bottommost row
	#the elements in this new row and column are determined by calculating all of the pairwise
	#joining rates for the new network with each of the old networks

	#first insert a column of zeroes on the right
	new_column = np.zeros(np.size(joining_matrix,1))
	joining_matrix = np.insert(joining_matrix,np.size(joining_matrix,1),new_column, axis=1)
	#print "joining matrix shape: ", np.size(joining_matrix,0), np.size(joining_matrix,1)
	#next insert a row of zeroes on the bottom
	new_row = np.zeros(np.size(joining_matrix,1))
	#print "new row: ", new_row
	joining_matrix = np.insert(joining_matrix,np.size(joining_matrix,0),new_row, axis=0)
	#now we are ready to insert new values for the joining rates into the matrix

	#we will calculate the values for the new row/column
	#we will loop through all of the networks in the system, excluding the reaction of the new network with itself, and then we will calculate the joining rates
	graphs = list(nx.connected_component_subgraphs(system))

	column = np.size(joining_matrix,1)-1
	new_network_lowest_node_number = min(new_network.nodes())
	for i in range(len(graphs)):
		network = graphs[i]
		#before we can add this to the joining_matrix we need to find the index for the network being considered
		lowest_node_number = min(network.nodes())
		if lowest_node_number == new_network_lowest_node_number:
			#print "considering joining of new network with itself"
			kjoin = 0
			#break
		else:
			kjoin = joining_rate(new_network, network)
			#print "joinning of the new network to the other networks is being accounted: ", kjoin 
			for key, value in updated_dict_2.iteritems(): 
				if value == lowest_node_number:
					index = key
					break
					#pdb.set_trace()
			joining_matrix[index,column] = kjoin 

	
	#now the joining matrix is completely up to date!
	#finally we will update the dictionary
	updated_dict_2_original = updated_dict_2.copy()
	updated_dict_2[column]=new_network_lowest_node_number
	joining_matrix_index_dictionary = updated_dict_2
	#pdb.set_trace()
	return joining_matrix, joining_matrix_index_dictionary



def perform_gillespie_simulation(i):
	#np.random.seed(1337)
	#initialize the system according the expt starting conditions
	system = initialize_system()

	#this can be adjusted to reflect the time measured in our experiments
	number_of_joining_steps = system.number_of_nodes()-100 #we will stop just short of everything in the system being joined to itself
	
	#these are stats that we will track over time
	average_nodes_per_network = [1]
	number_of_connected_graphs = [system.number_of_nodes()]
	average_seeds_per_network = [1]
	average_network_size_of_a_seed_list = [1]

	time = [0]

	#going to use a matrix to store joining probabilities and to compute and update rtot
	#just using an array/loops is not working out
	#this should be more efficient as well
	#this matrix will only apply to inter-network joining events
	#there will be a separate calculation for intra-network joining events and both will be accounted for in rtot
	#joining_matrix = np.zeros((system.number_of_nodes(),system.number_of_nodes()))
	joining_matrix, joining_matrix_index_dictionary = initialize_joining_matrix(system)
	print joining_matrix
	print joining_matrix_index_dictionary

	#need to be careful that system and joining_matrix are kept in-sync and that the indices in system connected graphs match those in joining matrix 
	#we can use is_isomorphic along with a node match function to determine if two networks are equivalent, this should circumvent the need to keep track of
	#different numbering in system vs the joining_matrix
	#need a dictionary to convert joining matrix numbering to system networks

	rtot = compute_rtot_from_joining_matrix(joining_matrix)

	print "rtot before loop starts: ", rtot 
	#joining_reactions_list, probabilities_list = generate_joining_reactions_list(system)
	for step in range(number_of_joining_steps):
		#print joining_matrix
		print "updated rtot: ", rtot
		#brute_force_rtot = compute_rtot_first(system)
		#print "brute force rtot: ", brute_force_rtot
		if rtot == 0.0:
			#no more reactions can occur
			break
		dt = compute_dt(rtot)
		#print dt 
		previous_time = time[len(time)-1]
		time.append(previous_time + dt)



		#now we need to determine which joining reaction will occur
		#the probability for a given reaction is that reaction's rate over rtot
		joined_network_1, joined_network_2, index_1, index_2 = choose_joining_reaction(rtot, system, joining_matrix, joining_matrix_index_dictionary)

		#now we actually perform the joining reaction by creating an edge between the two networks  
		system, new_network = perform_joining(joined_network_1, joined_network_2, system)

		new_joining_matrix, new_joining_matrix_index_dictionary = update_joining_matrix(joining_matrix, joining_matrix_index_dictionary, joined_network_1, joined_network_2, new_network, index_1, index_2, system)
		joining_matrix = new_joining_matrix
		joining_matrix_index_dictionary = new_joining_matrix_index_dictionary

		rtot_new = compute_rtot_from_joining_matrix(joining_matrix)
		rtot = rtot_new #updating rtot for the next iteration

		graphs = list(nx.connected_component_subgraphs(system))
		print "new number of connected graphs", len(graphs)
		number_of_connected_graphs.append(len(graphs))
		#average_seeds_per_network.append(50.0/float(len(graphs)))

		#now we will compute the average network size that a seed is in
		network_size_count = []
		for graph in graphs:
			for node in graph.nodes():
				network_size_count.append(graph.number_of_nodes())
		average_network_size_of_a_seed = float(sum(network_size_count))/float(len(network_size_count))

		average_network_size_of_a_seed_list.append(average_network_size_of_a_seed)

		#printing out data, do this every 200 time steps 
		if step%20 == 0:

				print "printing current simulation time"
				f1=open('simulation_time.dat', 'a')
				print >>f1, previous_time 
				f1.close()

				print "average network size of a seed"
				f1=open('average_network_size_of_a_seed.dat','a')
				print >>f1, average_network_size_of_a_seed
				f1.close()

				print "printing network sizes squared"
				f1=open(str(previous_time)+"_network_size_count_squared.dat", 'w+')
				for size in network_size_count:
					print >>f1, size
				f1.close()

				print "printing network sizes"
				f1=open(str(previous_time)+"_network_size_count.dat", 'w+')
				for graph in graphs:
					print >>f1, graph.number_of_nodes()
				f1.close()



	#print time
	#print number_of_connected_graphs
	#print average_seeds_per_network
	#print average_network_size_of_a_seed_list

	'''print "average_network_size_of_a_seed"
	f1=open('average_network_seeds.dat','w+')
	for network_seed in average_network_size_of_a_seed_list:
		print >>f1, network_seed
	f1.close()'''



'''	plt.plot(time, average_network_size_of_a_seed_list, label = 'Gillespie algorithm run '+str(i+1), linewidth = 5, alpha = .5)
	

for i in range(1):
	perform_gillespie_simulation(i)

expt_times = [0.5, 2.5, 4.5, 6.5, 8.5]
expt_seeds_per_network = [2.45, 4.68, 6.84, 19.33, 17.32]
expt_errors = [.045, .204, .280, 1.12, 1.14]

#expt_average_network_size_of_seed = [1.51, 2.42, 3.61, 4.34]

plt.errorbar(expt_times, expt_seeds_per_network, expt_errors, label = 'experiment', linewidth = 5, color = 'black' )
#plt.xlim([0.0, 10])
#plt.ylim([0.0, 30])
plt.legend(loc = 4, fontsize = 8)
plt.xlabel('time (hours)')
plt.ylabel('average network size that a seed is in')
plt.savefig("gillespie_results.pdf")
print datetime.now() - startTime 

#5*5*9 + 10*6

#6*6*9'''

if __name__ == "__main__":
	perform_gillespie_simulation(1)






	

