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

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from numpy import linalg
from basis_functions import *
from utils import *

class element:

	# Define some convenient constant attributes
	num_vertex = 4 # the four vertices of a quadrilateral element

	# The parametric and physical space dim
	space_dim_param = 2
	space_dim = 2
 
	## @brief the constructor for a quadrilateral element 
	#  @param[in] element_index the index of the element
	#  @param[in] coordinates the coordinates of the mesh nodes (passed by reference)
	#  @param[in] connectivity the connectivity of the mesh (passed by reference)
	def __init__(self, element_index, coordinates, connectivity, poly_order):

		# The element index
		self.element_index = element_index 

		# Get pointers to coordinates and connectivity arrays
		# note here that these are passed by reference
		self.coordinates = coordinates 
		self.connectivity = connectivity 

		# The polynomial order 
		self.poly_order = poly_order 

		# The number of degrees of freedom
		self.num_dofs =  pow(self.poly_order + 1, self.space_dim_param) 

		# Get the coordinates of the nodes in the parametric domain
		# this are the nodes of the exterior vertices, edge nodes,
		# as well as interior nodes
		self.nodes_crds_param = np.linspace(-1., 1.,  poly_order +1 )
		X,Y = np.meshgrid( self.nodes_crds_param, self.nodes_crds_param )
		self.nodes_crds_param =  np.c_[ X.flatten() , Y.flatten() ]

		# Construct the array containing the physical coordinates 
		# corresponding to the coordinates of each degree of freedom
		self.interpolate_interior_nodes()
	
	## @brief only relevant for higher order elements this method interpolates
	#  linearly the mid and interior nodes on the physical domain. If we want to 
	#  specify the value we can set the nodes using the 
	def interpolate_interior_nodes(self):

		# Here we are simply interpolating linearly the physical cordinates of the interior nodes
		# this could be more general but for the purpose of the homework it will suffice
		crds_vertex = self.coordinates[ self.connectivity[ self.element_index ] ]
		self.nodes_crds_phys = np.zeros( self.nodes_crds_param.shape )

		# For each of all nodes (vertex edge or internal) 
		for a in range( len( self.nodes_crds_phys ) ): 

			# Interpolate the physical coordinates of edge and interior nodes linearly
			for i in range( self.num_vertex ):
				self.nodes_crds_phys[ a ] += hexahedral_base(1,i,self.nodes_crds_param[a])*\
						crds_vertex[ abs( 3*(i//2) - i%2)  ]

		return

	## @brief it overloads the assumption that the edge nodes lie along a line 
	#  allowing for curved boundaries
	#  @param[in] nodes_crds_phys the array of nodal coordinates of dimension 
	#  pow(poly_order+1,space_dim) x space_dim
	def set_interior_nodes_coordinates(self, nodes_crds_phys ):

		assert nodes_crds_phys.shape == ( pow(self.poly_order+1, self.space_dim) , self.space_dim ),\
		'There must be a coordinate of space dim for each node, including vertex nodes'
		
		self.nodes_crds_phys = nodes_crds_phys
		return 

	## @brief The map from the parametric domain to the physical domain
	#  @param[in] y the vector values position in the parametric domain
	def get_map(self, y ):

		# Compute the map  
		x = 0

		# Here we are using an isoparametric map
		for i in range(self.num_dofs):
			# Get the Gauss point position in physical space 
			x += hexahedral_base(self.poly_order,i,y)*self.nodes_crds_phys[i]
		
		return x

	## @brief The derivative of the map from the parametric domain to the physical domain
	#  @param[in] y the vector values position in the parametric domain
	#  @param[in] jacobian whether to compute the jacobian or not.
	def get_dmap(self, y, jacobian=False ):

		d_x = np.zeros( (2,2) )
        
		# Here we are using an isoparametric map
		for i in range(self.num_dofs):
			#print('shape of is:', self.nodes_crds_phys.shape)
			# Get the jacobian of the mapping at the gauss point
			d_x += np.outer(self.nodes_crds_phys[i],\
					grad_hexahedral_base(self.poly_order,i, y ) )            
				

		if jacobian: 

			# Get the jacobian
			jacobian = np.linalg.det(d_x)
			return jacobian


		return d_x

	## @brief returns the value of a base function at a secific point
	#	The bases are numbered consequentially. Eg. forquadratic elements 
	# 6--7--8
	# |     |
	# 3  4  5
	# |     |
	# 0--1--2
	# @param[in] base_index the index of the base 
	# @param[in] y the parametric coordinated 
	# @returns the value of the shape function at y
	def get_base_function_val( self, base_index, y ):
		return hexahedral_base(self.poly_order,base_index,y)

	## @brief returns the value of the gradient with respect to the parametric coordinates
	#  of a base function at a secific point
	#  @param[in] base_index the index of the base 
	#  @param[in] y the parametric coordinated 
	#  @returns the value of the gradient of the shape function at y
	def get_base_function_grad_parametric( self, base_index, y ):
		return grad_hexahedral_base(self.poly_order,base_index,y)

	## @brief returns the value of the gradient with respect to phisical coordinates
	#  of a base function at a secific point
	#  @param[in] base_index the index of the base 
	#  @param[in] y the parametric coordinated 
	#  @returns the value of the gradient of the shape function at y
	def get_base_function_grad( self, base_index, y ):

		# Get the derivative of the mapping from parametric to physical
		d_x = np.zeros( (2,2) )
		d_x = self.get_dmap( y )
		# Get the inverse of the tangent map
		d_x_inv = np.linalg.inv( d_x )
		
		# Get the gradient of the basis functions with respect to parametric coordinates
		grad_phi_param = self.get_base_function_grad_parametric( base_index, y )

		# Compute gradient with respect to physical coordinates
		grad_phi = np.dot( d_x_inv.T, grad_phi_param )
		
		return grad_phi

	## @returns returns the element total number of degrees of freedom
	def get_num_dofs(self):
	
		return self.num_dofs

	## @brief gets the quadrature for a quadrilateral element as a result of tensor products 
	#  of Legendre quadrature in one dimension
	#  @param[in] quadrature_order the order of the quadrature
	#  @retur a tuple containing the array of gauss point coordinates and weights
	def get_quadrature(self, quadrature_order=0  ):

		# If the quadrature order is not specified use one that is sufficient in 
		# correlation to the degree of the polynomial interpolant
		if not quadrature_order:
			quadrature_order = int(np.ceil( (self.poly_order+1)//2 )+1)

		# Get the quadrature rule for the interval [-1,1]
		gauss_points, gauss_weights = np.polynomial.legendre.leggauss( quadrature_order  )

		# Get the quadrature rule for the interval [-1,1]x[-1,1]
		X,Y = np.meshgrid( gauss_points, gauss_points )
		gauss_points = np.c_[ X.flatten() , Y.flatten() ]
		gauss_weights = np.outer( gauss_weights, gauss_weights ).flatten()

        
		return gauss_points, gauss_weights


