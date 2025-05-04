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

from utils import lagrange_basis, d_lagrange_basis

def hexahedral_base( poly_order, base_index, x ):
	'''
	Parameters
	----------
	poly_order: int 
		the order of the polynomial interpolant
	base_index: int
		the index of the basis function

	Returns
	-------
	The value of the base function corresponding to 
	base_index defined sover parameteric space 
	$[-1,1]\times[-1,1]$ evaluated at point x

	The bases are numbered consequentially. Eg. for
	quadratic elements 
	6--7--8
	|     |
	3  4  5
	|     |
	0--1--2
	'''

	# The nodes of the one dimensional polynomial 
	# over the parametric interval [-1,1]
	param_one_dim = np.linspace(-1., 1., poly_order + 1 )

	# Get the index in the xi_1 direction 
	base_index_xi_1 = base_index%(poly_order + 1 )

	# Get the index in the xi_2 direction 
	base_index_xi_2 = base_index/(poly_order + 1)

	return lagrange_basis(param_one_dim,base_index_xi_1,poly_order,x[0])*\
						lagrange_basis(param_one_dim,base_index_xi_2,poly_order,x[1])

def grad_hexahedral_base( poly_order, base_index,x ):
	'''
	Parameters
	----------
	poly_order: int 
		the order of the polynomial interpolant
	base_index: int
		the index of the basis function
	
	Returns
	-------
	The value of the gradient of the base function corresponding to 
	base_index defined sover parameteric space 
	$[-1,1]\times[-1,1]$ evaluated at point x
	
	The bases are numbered consequentially. Eg. for
	quadratic elements 
	6--7--8
	|     |
	3  4  5
	|     |
	0--1--2
	'''

	# The nodes of the one dimensional polynomial 
	# over the parametric interval [-1,1]
	param_one_dim = np.linspace(-1., 1., poly_order + 1 )

	# Get the index in the xi_1 direction 
	base_index_xi_1 = base_index%(poly_order + 1 )

	# Get the index in the xi_2 direction 
	base_index_xi_2 = base_index/(poly_order + 1)

	# Return the gradient
	return np.array([ d_lagrange_basis(param_one_dim,base_index_xi_1,poly_order,x[0])*\
						lagrange_basis(param_one_dim,base_index_xi_2,poly_order,x[1]),\
						lagrange_basis(param_one_dim,base_index_xi_1,poly_order,x[0])*\
						d_lagrange_basis(param_one_dim,base_index_xi_2,poly_order,x[1]) ])
