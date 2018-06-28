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
















	

