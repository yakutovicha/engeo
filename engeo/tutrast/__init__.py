"""## PUT DESCRIPTION HERE##"""

from __future__ import print_function
import sys
import numpy as np

class Basin:
    """Class that works with basins. It only deals with the basin's geometry. Does not know
    anything about the data stored in it."""
    neighbor_operations = [(0,  1), (0, -1), (1,  1), (1, -1), (2,  1), (2, -1)]
    def __init__(self, indeces, shape, sorted_coordinates, coordinates_to_indices):
        self.indeces = set(indeces)
        self.do_not_have_all_neighbors = set(indeces)
        self.x_size = shape[0]
        self.y_size = shape[1]
        self.z_size = shape[2]
        self.coordinates_to_indices = coordinates_to_indices
        self.sorted_coordinates = sorted_coordinates

    def append(self, indeces):
        self.indeces.update(indeces)

    def choose_neighbors_from_list(self, points):
        chosen = []
        for indx, pnt in enumerate(points):
            if self.neighbors(pnt) & self.indeces: # if neighbors of a point and elements of the basin have points in
            # common, then the point is at the border of the basin
                chosen.append(indx)
        return points[chosen], np.delete(points, chosen)

    def neighbors(self, point):
        point_coords = self.sorted_coordinates[point]
        neighbors_coords = np.array([point_coords] * 6)
        for i, opn in enumerate(self.neighbor_operations):
            neighbors_coords[i][opn[0]] += opn[1]

        # apply periodic boundary conditions (PBC)
        neighbors_coords[:,0] %= self.x_size
        neighbors_coords[:,1] %= self.y_size
        neighbors_coords[:,2] %= self.z_size
        return set(self.coordinates_to_indices(neighbors_coords))


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

        # sort energy array and index array according to the energy values
        self.forward_indx_permutation = np.argsort(flat_e) # order indeces by energy
        self.reverse_indx_permutation = np.argsort(self.forward_indx_permutation) # create back permutation array

        # make the array of energies and array of coordinates sorted by energy.
        self.flat_e = flat_e[self.forward_indx_permutation] # apply new ordering to the energy array
        self.sorted_coordinates = coordinates[self.forward_indx_permutation] # apply the same ordering to the

        # coordinates array
        self.basins = [] # [[bassin 0 point flat indices], [bassin 1 point flat indices], ...]

    def _chunk_energy_and_indeces(self):
        index_min=0
        for i, energy in enumerate(self.flat_e):
            if energy - self.flat_e[index_min] > self.step:
                yield self.flat_e[index_min], np.array(range(index_min, i))
                index_min = i


    def find_basins(self):
        """TODO: The function is still under development"""
        # flat_indeces reference to all the energies in the range [e_current; e_current+step]
        total_len = 0
        for e_current, flat_indeces in self._chunk_energy_and_indeces():
            n_basins = len(self.basins)
            n_points = len(flat_indeces)
            total_len += n_points
            print("n_basins", len(self.basins), "n_points", n_points, "total_points", total_len)  
            #if (total_len > 40000):
            #    sys.exit(0)

            # add points to the existing basins until it's possible
            while True: # infinite loop
                added_elements = False
                for bsn in self.basins:
                    neigh_indxs, flat_indeces = bsn.choose_neighbors_from_list(flat_indeces)
                    if neigh_indxs.size > 0:
                        bsn.append(neigh_indxs)
                        added_elements = True
                    if flat_indeces.size == 0:
                        break
                if not added_elements or flat_indeces.size == 0:
                    break

            # form clusters from the leftover elements in the flat_indeces array
            for cluster in self.form_clusters(flat_indeces):
                self.basins.append(Basin(
                                    indeces = cluster,
                                    shape = (self.x_size, self.y_size, self.z_size),
                                    sorted_coordinates = self.sorted_coordinates,
                                    coordinates_to_indices = self.coordinates_to_indices,
                                )
                            )

            # TODO: merge basins if necessary
        tot = 0
        for i, b in enumerate(self.basins):
            print("basin_{} {} elements".format(i, len(b.indeces)))
            tot += len(b.indeces)
        print(tot)

    def form_clusters(self, indeces):
        """The function gets a list of point indeces and clusters them if they are neighbors.

        :param indeces: list of indeces
        :type indeces: list
        """
        
        # neighbor operations
        # first element: 0 means 'x', 1, means 'y', 2 means 'z'
        # second element: add or subtract 1.
        neighbor_operations = [(0,  1), (0, -1), (1,  1), (1, -1), (2,  1), (2, -1)]

        # list of clusters
        clusters = []
        for point, point_index in zip(self.sorted_coordinates[indeces], indeces):
            point = tuple(point) # make it tuple to be hashable
            goes_in_clusters = [] # to which clusters the point should go
            for nclstr, cluster in enumerate(clusters):
                for opn in neighbor_operations:
                    candidate = list(point)
                    candidate[opn[0]]+= opn[1]

                    # apply periodic boundary conditions (PBC)
                    candidate[0] %= self.x_size
                    candidate[1] %= self.y_size
                    candidate[2] %= self.z_size

                    # if a neighbor of the current poin is in a cluster, then the point can be added to the cluster
                    if tuple(candidate) in cluster:
                        goes_in_clusters.append(nclstr)

            goes_in_clusters = list(set(goes_in_clusters)) # make goes_in_clusters contain only unique elements
            # if the point is connected to one cluster only, then attach it
            if len(goes_in_clusters) == 1:
                clusters[goes_in_clusters[0]][point] = point_index
            # if the point is connected to two and more clusters, then attach it and merge all the clusters
            elif len(goes_in_clusters) > 1:
                clusters[goes_in_clusters[0]][point] = point_index
                for i in goes_in_clusters[1:]: # merge the clusters
                    clusters[goes_in_clusters[0]] = {**clusters[goes_in_clusters[0]], **clusters[i]}
                for i in sorted(goes_in_clusters[1:], reverse=True): # delete old clusters
                    clusters.pop(i)
            else:
                clusters.append({point:point_index})

        for cluster in clusters:
            yield [indx for indx in cluster.values()]

    def coordinates_to_indices(self, coords):
        """Converts list of coordinates [[x1, y1, z1], ...] to the list of indeces[i1, i2, ...] in the sorted array

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
        # this variant turned out to be the fastest:
        #indxs = np.array([c[0] * self.y_size * self.z_size + c[1] * self.z_size  + c[2] for c in coords])
        #indxs = list(map(lambda c: c[0] * self.y_size * self.z_size + c[1] * self.z_size  + c[2], coords))
        #print(type(coords))
        #indxs = np.array(map(lambda c: c[0] * self.y_size * self.z_size + c[1] * self.z_size  + c[2], coords))
        indxs = coords[:, 0] * self.y_size * self.z_size + coords[:, 1] * self.z_size + coords[:, 2]
        return self.reverse_indx_permutation[indxs] # indeces in the sorted array
