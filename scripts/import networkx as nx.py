import networkx as nx
import numpy as np
import matplotlib.pyplot as plt 


#writing a script to grow a nanotube network, the inputs will be the relative
#abundance of 0,1,2, and 3 armed structures (nodes)
#the graph will be allowed to grow until there are no suitable points to add additional edges
#then the goal will be to analyze the size of the network (number of added nodes) as a function of the relative 
#abundances of each species, this could definitely be extended to include L and T junctions if needed! 

###############
#hard coded relative abundances
p_3arm = .33
p_2arm = .33
p_1arm = .34

connectivities = [1, 2, 3]
probabilities = [p_1arm, p_2arm, p_3arm]

#initializing a directed graph
G = nx.DiGraph()

#initializing the first node
first_node_connectivity = np.random.choice(connectivities, probabilities)
G.add_node(1, arms = first_node_connectivity, terminal = False)

#adding the first shell of connecting nodes
for i in range(first_node_connectivity):
	n_arms = np.random.choice(connectivities, probabilities)

	if n_arms == 1:
		G.add_node(G.numer_of_nodes()+1, arms = 1, terminal = True)
		G.add_edge(1, G.numer_of_nodes()+1)
	else:
		G.add_node(G.numer_of_nodes()+1, arms = n_arms, terminal = False)
		G.add_edge(1, G.numer_of_nodes()+1)

plt.subplot(121)
nx.draw(G, with_labels = True)
plt.show()

#now looping through many iterations, each iteration:
# 1. finds the nodes that were added during the last iterations
# 2. based on the number of arms for each node that was added in the last iteration, add the appropriate
#    number of nodes or terminate growth from a particular node 
# 2.5. if no nodes allow for the addition of new nodes, terminate the growth process
# 3. Add the appropriate edges to each new node
#for i in range(10000000000):
	#need to loop over all non-terminal nodes added during the last iteration
	#in other words, loop over all nodes that do not have children and are not terminal

	#grabbing a list of nodes that do not have childern


