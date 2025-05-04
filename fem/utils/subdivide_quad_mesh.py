"""
Created by Maurizio Chiaramonte 
	  
Copyright (c) 2017 Maurizio Chiaramonte

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
"""

import numpy as np
import matplotlib.pyplot as plt 

def subdivide_mesh(coordinates, connectivity):
	
	print(len(connectivity))
	# A dictionary indexed by the sorted indeces of
	# the edge nodes that stores the edge mid node number 
	edge_list = {}

	# The edges connectivity array
	edge_connectivity = [] 

	# The new connectivity array
	new_connectivity = [] 
	
	# The total number of nodes
	num_nodes = len( coordinates )

	# Add a center node for each element 
	for e,conn_e in enumerate(connectivity):
		coordinates = np.r_[coordinates, [1./4*np.sum( coordinates[ conn_e ], axis=0 ) ] ] 
	
	# Add the edge nodes 
	for e,conn_e in enumerate(connectivity): 

		# The number of nodes per element
		num_nodes_element = len(conn_e )

		# Create a connectivity array for edges
		edge_connectivity.append([])

		# Loop over all edges
		for i in range( num_nodes_element ):

			# Get the sorted edge nodes
			edge = tuple(np.sort(( conn_e[ i ] , conn_e[ (i+1)%num_nodes_element ] )))
			
			# If the edge is already in the dictionary
			if edge in edge_list: 
				
				# Add the edge label to the edge connectivity
				edge_connectivity[-1].append( edge_list[edge]  )
				continue 

			# Add to the edge dictionary the edge label and the coordinate of the midpoint
			edge_list[edge] = len( coordinates )

			# Add mid node to coordinates list 
			coordinates = np.r_[ coordinates, [1./2*( coordinates[ conn_e[ i%num_nodes_element ] ]  \
				+  coordinates[ conn_e[ (i+1)%num_nodes_element ] ] ) ]  ]

			# Add to the edge connectivity of the element 
			edge_connectivity[-1].append( edge_list[edge] )

	# For each edge add a mid node 
	for e,conn_e in enumerate(connectivity): 
		edge_conn_e = edge_connectivity[e] 
		center_node = num_nodes + e 
		new_connectivity.append([ center_node, edge_conn_e[3], conn_e[0] , edge_conn_e[0] ]) 
		new_connectivity.append([ center_node, edge_conn_e[0], conn_e[1] , edge_conn_e[1] ]) 
		new_connectivity.append([ center_node, edge_conn_e[1], conn_e[2] , edge_conn_e[2] ]) 
		new_connectivity.append([ center_node, edge_conn_e[2], conn_e[3] , edge_conn_e[3] ]) 
	
	# Assign the new connectivity to the connectivity array
	connectivity = np.array(new_connectivity)

	return coordinates, connectivity


