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
sys.path.append(os.getcwd()+'/../../')

import matplotlib.pyplot as plt
from quadrilateral_element import *


if __name__ == "__main__":
	''' 
	tests the above implementation of the local to global map
	'''

	# Define the order of the polynomial interpolant
	poly_order = 3

	# Define the half width of the domain
	w = 2.

	# Create the mesh of a square
	coordinates = np.array([ (-w,-w), (w,-w), (w,w), (-w,w) ])  	
	
	# Connectivity 
	connectivity = np.array([[0,1,2,3]])

	# Subdivide mesh
	#coordinates, connectivity = subdivide_mesh(coordinates, connectivity)

	l2g = local_to_global_map( poly_order, connectivity )

	# Print out 
	print 'Total dofs : %i Should be: %i ' %( l2g.get_total_dof() , pow(poly_order+1,2) )

	xs = np.linspace(-1,1,poly_order + 1 )
	X,Y = np.meshgrid( xs,xs )

	crds_param = np.c_[ X.flatten(), Y.flatten() ]
	
	for e,conn_e in enumerate(connectivity):

		crds = coordinates[ conn_e[0] ] + (crds_param + 1 )/2.*( coordinates[ conn_e[2] ] - coordinates[ conn_e[0] ] )
		plt.plot( crds[:,0],crds[:,1] , 'o' )

		for i in range(len(crds)): 
			plt.text( crds[i,0] , crds[i,1], '%s'%l2g.get_global_dof( e, i ) )
			#plt.text( crds[i,0] - 0.1 , crds[i,1], '%i'%i )

	plt.show()
