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
import time

# @brief assembles the local arrays into the global arrays
# @param[in] local_to_global_map the local to global map
# @param[in] element the element object
# @param[in] ke the element stiffness matrix
# @param[in] fe the element stiffness matrix
# @param[out] K the reference to the global stiffness matrix
def assemble_local(local_to_global_map, element, ke, me, fe, K, M, F ):

	# Loop over all degrees of freedom 
	for i in range(element.num_dofs):

		# Get the global index of the local i dof
		i_global = local_to_global_map.get_global_dof( element.element_index , i )
		#i2_global = local_to_global_map.get_global_dof( element.element_index , i )
        
		# Add the contribution to the source vector
		F[i_global] += fe[i]

		for j in range(element.num_dofs):

			# Get the global index of the local j dof
			j_global = local_to_global_map.get_global_dof( element.element_index , j )

			# Add the contribution of the local stiffness to the global stiffness matrix		
			K[i_global, j_global] += ke[ i, j]
			M[i_global, j_global] += me[ i, j]

	return


# @brief assembles the local arrays into the global arrays
# @param[in] local_to_global_map the local to global map
# @param[in] element_operation the element operation 
# @param[in] elements the aray containing the lement objects 
# @param[in] ke the element stiffness matrix
# @param[in] fe the element stiffness matrix
# @param[out] K the reference to the global stiffness matrix
# @param[out] F the reference to the global source vector
def assemble_global(local_to_global_map, element_operation, elements, K, M,F ):

	time_elt_stiff = 0
	time_assemble = 0

	# Loop over all elements to assemble the global stiffness matrix
	for element in elements:

		# The element stiffness matrix and source vector
		t0 = time.time()
		ke, me, fe = element_operation.get_element_arrays(element)
		time_elt_stiff += time.time() - t0

		# Assemble them 
		t0 = time.time()
		assemble_local(local_to_global_map, element, ke, me, fe, K, M, F )
		time_assemble += time.time() - t0
	print (time_elt_stiff , time_assemble)

	return

