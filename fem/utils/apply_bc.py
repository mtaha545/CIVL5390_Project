# @brief applies boundary conditions where appropriate
# @param[in] dofs the degree of freedom indeces of the boundary conditions
# @param[in] values thevalues of the boundary conditions
# @param[out] K the global stiffness matrix 
# @param[out] F the global forcing vector
def apply_bc(dofs, values, K, M, F ):

	# Loop over all degrees of freedom
	for i in range(len(dofs)):
		
		# Zero out the corresponding row 
		K[dofs[i],:] = 0 

		# Place one on the diagonal 
		K[dofs[i],dofs[i]] = 1
        
		# Zero out the corresponding row 
		M[dofs[i],:] = 0 

		# Place one on the diagonal 
		M[dofs[i],dofs[i]] = 1
		# Place the value of the boundary in the forcing vector
		F[dofs[i]] = values[i]
	return
