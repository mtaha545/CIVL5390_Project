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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri

import numpy as np

def plot_contour( coordinates, values, space_dim= 2 ):

	assert space_dim == 2, 'only works for two space dim'

	x = coordinates[:,0] 
	y = coordinates[:,1] 
	
	# Create the Triangulation just for plotting
	triang = tri.Triangulation(x, y)

	#	Plot contour
	plt.figure()
	plt.gca().set_aspect('equal')
	plt.tricontourf(triang, values)
	
	return
    
        
def plot_quad_mesh( coordinates, connectivity ):
	
	for conn_e in connectivity: 
		
		crds = coordinates[ conn_e ] 

		plt.plot( crds[:,0], crds[:,1] , '-k',linewidth=1)

	return

