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

# @param[in] nodes the coordinates of the nodes
# @param[in] index the index of the basis function
# @param[in] order the order of the polynomial basis
# @param[in] x the coordinate where to evaluate the basis
# @return the value of the index-th basis function at point x
def lagrange_basis(nodes, index, order, x): 

	# Check that we have the right number of nodes
	assert len(nodes) == (order + 1)

	# Create the initial value of the function
	ell = 1.

	# Loop over all nodes
	for i in range(0,order+1):
		index = np.floor(index)
		index = int(index)
		# If this node is the same as the 
		# support node of the basis function, skip it
		if i == index:
			continue
		# Otherwise perform the multiplication
		#print('i is:', i, 'Index is: ', index)
		ell *= ( x - nodes[i])/(nodes[index] - nodes[i])
	
	return ell	

# @param[in] nodes the coordinates of the nodes
# @param[in] index the index of the basis function
# @param[in] order the order of the polynomial basis
# @param[in] x the coordinate where to evaluate the basis
# @return the value of the index-th basis function at point x
def d_lagrange_basis(nodes, index, order, x): 
	index = int(index)
	# Check that we have the right number of nodes
	assert len(nodes) == (order + 1)

	# Create the initial value of the function
	d_ell = 0

	# Loop over all nodes
	for i in range(0,order+1):
		
		# If this node is the same as the 
		# support node of the basis function, skip it
		if i == index:
			continue

		# Otherwise perform the multiplication
		d_ell += 1./(x - nodes[i]) 
	
	# Multiply by the corresponding Lagrange basis
	d_ell *= lagrange_basis(nodes, index, order, x)
		
	return d_ell	

