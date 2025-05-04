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
sys.path.append(os.getcwd()+'/../')

# @brief the parent class for element operations
class element_operation: 
    
	# @brief the method should compute the element arrays 
	# @param[in] the element object
	def get_element_arrays(self, element ):

		raise NotImplementedError()
