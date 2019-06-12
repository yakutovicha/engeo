"""## PUT DESCRIPTION HERE##"""

from __future__ import print_function
import numpy as np

class TuTraSt:
	KCAL_MOL_TO_HARTREE = 0.00159360109742136
	KJ_MOL_TO_HARTREE = 0.0003808798033989866
	EV_TO_HARTREE = 0.03674930495120813
	def __init__(self, data, step=0.1, units='kcal/mol'):

		data = np.copy(data)
		# dealing with units
		if units == 'a.u.':
			conversion = 1.0
		elif units == 'kcal/mol':
			conversion = self.KCAL_MOL_TO_HARTREE
		elif units == 'kJ/mol':
			conversion = self.KJ_MOL_TO_HARTREE
		elif units == 'eV':
			conversion = self.EV_TO_HARTREE
		self.step = step * conversion
		data *= conversion

		# create flat index array
		self.x_size, self.y_size, self.z_size = data.shape
		xindeces = np.arange(self.x_size)
		yindeces = np.arange(self.y_size)
		zindeces = np.arange(self.z_size)
		coordinates = np.stack(np.meshgrid(xindeces, yindeces, zindeces)).transpose((2,1,3,0))
		coordinates = coordinates.reshape(self.x_size*self.y_size*self.z_size, 3)

		# create flat Energy array
		flat_e = data.flatten()

		# TODO: move this to the tests?
		# check if old and new energy array match
		for i, indx in enumerate(coordinates):
			if flat_e[i] != data[tuple(indx)]:
			    print(flat_e[i], data[tuple(indx)], indx, i)

		# sort energy array and index array according to the energy values
		self.forward_indx_permutation = np.argsort(flat_e) # order indeces by energy
		self.reverse_indx_permutation = np.argsort(self.forward_indx_permutation) # create back permutation array
		self.flat_e = flat_e[self.forward_indx_permutation] # apply new ordering to the energy array
		self.sorted_coordinates = coordinates[self.forward_indx_permutation] # apply the same ordering to the
		# coordinates array

	def _chunk_energy_array(self):
		energy_array = [self.flat_e[0]]
		for energy in self.flat_e[1:]:
			if energy - min(energy_array) < self.step:
				energy_array.append(energy)
			else:
				yield np.array(energy_array)
				energy_array = [energy]
		

	def find_bassins(self):
		"""WARNING: The function is still under development"""
		bassins = [] # [[bassin 0 point flat indices], [bassin 1 point flat indices], ...]

		# step 1: find all points in the range [E_current; E_current+step]
		for flat_indeces in self._chunk_energy_array():
			# look in the `self.sorted_coordinates` and try to find neighboring points of known bassins

			for bssn in bassins:
				# growing one layer
				neighbor_coordinates = get_neighbors(bssn)
				neighbor_indeces = coordinates_to_indices(neighbor_coordinates)

				# check for transition state
				# if TS  found - store the index
 

	def get_neighbors(self, bassin):
		"""Gets coordinates of the bassin neighboring points.

		:param bassin: grid point indeces in the sorted arrays
		:type bassin: 1D array
		"""
		all_neighbors = np.array([])
		bassin_points_coordinates = self.sorted_coordinates[bassin]

		# neighbor operations
		# first element: 0 means 'x', 1, means 'y', 2 means 'z'
		# second element: add or subtract 1.
		neighbor_operations = [(0,  1), (0, -1), (1,  1), (1, -1), (2,  1), (2, -1)]

		for opn in neighbor_operations:
			new_nbrs = np.copy(bassin_points_coordinates)
			new_nbrs[:,opn[0]] += opn[1]
			all_neighbors = np.append(all_neighbors, new_nbrs, axis=0) if all_neighbors.size > 0 else new_nbrs

		# apply periodic boundary conditions (PBC)
		all_neighbors[:,0] %= self.x_size
		all_neighbors[:,1] %= self.y_size
		all_neighbors[:,2] %= self.z_size

		# make sure the all neighbors coordinates are unique
		all_neighbors = np.unique(all_neighbors, axis=0)

		# collect neighbors and bassins coordinates in one array
		# append self.sorted_coordinates[bassin] twice to make them non-unique
		all_neighbors = np.append(all_neighbors, bassin_points_coordinates, axis=0)
		all_neighbors = np.append(all_neighbors, bassin_points_coordinates, axis=0)

		# extract only the unique ones
		unq, count = np.unique(all_neighbors, return_counts=True, axis=0)
		return unq[count == 1] # return elements that appear one time only

	def coordinates_to_indices(self, coords):
		"""Converts list of coordinates [[x1, y1, z1], ...] to the list of indeces[i1, i2, ...].

		:param coords: a list of point coordinates
		:type: 2D array, [[x1, y1, z1], ...]

		:returns: a list of indeces of points in sorted array
		:rtype: 1D array

		In essensce this is what happens here:
		[] - an array
		() - a permutation
             0  1  2  3             0  1  2  3
		1. [ 7, 5, 8, 1] <--->  3. (2, 1, 3, 0)
                 |                      | 
                 |                      |
		         |                      |
            0  1  2  3              0  1  2  3
        2. (3, 1, 0, 2)  <--->  4. [1, 5, 7, 8]
        
        steps:
        * select element from 1., let it be 8 (comes as `coords` array)
        * take its index: 2
        * go to permutation 3. and find element at position 2: 3
        * go to position 3 in 4.: 8
		"""
		# TODO: add a check that it is 2D array
		indxs = np.array([c[0]*self.x_size + c[1]*self.y_size + c[2] for c in coords]) # indeces in the initial array
		return self.reverse_indx_permutation[indxs] # indeces in the sorted array
