"""## PUT DESCRIPTION HERE##"""

from __future__ import print_function
import sys
import numpy as np

class Basin:
    """Class that works with basins. It only deals with the basin's geometry. Does not know
    anything about the data stored in it."""
    neighbor_operations = [(0,  1), (0, -1), (1,  1), (1, -1), (2,  1), (2, -1)]
    def __init__(self, indices, ttst, tolerance=0.1):
        self.indices = set(indices)
        self.do_not_have_all_neighbors = set(indices)
        self.x_size = ttst.x_size
        self.y_size = ttst.y_size
        self.z_size = ttst.z_size
        self.coordinates_to_indices = ttst.coordinates_to_indices
        self.sorted_values = ttst.flat_values
        self.sorted_coordinates = ttst.sorted_coordinates
        self.tolerance = tolerance
        self.ts = [] # list of transition states bind to this basin

    def append_ts(self, ts):
        self.ts.append(ts)

    def append(self, indices):
        self.indices.update(indices)

    def choose_neighbors_from_list(self, points):
        chosen = []
        for indx, pnt in enumerate(points):
            if self.neighbors(pnt) & self.indices: # if neighbors of a point and elements of the basin have points in
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

    def is_neighbor_with(self, bsn):
        bsn1, bsn2 = (self, bsn) if len(bsn.indices) > len(self.indices) else (bsn, self) # take the smaller basin for the sake of efficeincy
        return True if bsn1.get_neighbors_indices() & bsn2.indices else False

    def get_neighbors_indices(self):
        """Gets coordinates of the basin neighboring points."""
        potential_neighbors = set()
        for indx in self.indices:
            potential_neighbors.update(self.neighbors(indx))
        return potential_neighbors - self.indices

    def merge_basin(self, basin):
        self.append(basin.indices)

    def close_by_energy(self, basin):
        return True if abs(self.max_value - basin.min_value) < self.tolerance or abs(self.min_value
            - basin.max_value) < self.tolerance else False
    @property
    def is_definitely_basin(self):
        return True if self.max_value - self.min_value > self.tolerance else False

    @property
    def min_value(self):
        return self.sorted_values[self.min_index]
    
    @property
    def max_value(self):
        return self.sorted_values[self.max_index]
    
    @property
    def min_index(self):
        return min(self.indices)

    @property
    def max_index(self):
        return max(self.indices)

class TransitionState:
    # TODO: the TransitionState class is still under development, do not test
    def __init__(self, basin1, basin2, ttst):
        self.border = []
        self.basin1 = basin1
        self.basin2 = basin2
        self.ttst = ttst
        basin1.append_ts(self)
        basin2.append_ts(self)
        ts_indeces_bsn1 = basin2.get_neighbors_indices() & basin1.indices

        # TODO, add a test that TS class makes the same border even if one exchanges
        # basin1 and basin2
        for border_point_1 in ts_indeces_bsn1:
            nbrs = basin1.neighbors(border_point_1)
            for border_point_2 in nbrs & basin1.indices:
                self.border.append((border_point_1, border_point_2))
    @property
    def area(self):
        """Computes the area of the Transition State"""
        surface = 0

        # extracting basins' indices
        basin1_points = np.array(self.border)[:,0]
        basin2_points = np.array(self.border)[:,1]

        # getting their coordinates
        basin1_coords = self.ttst.sorted_coordinates[basin1_points] # p1_coords should be a numpy array
        basin2_coords = self.ttst.sorted_coordinates[basin2_points] # p2_coords should be a numpy array

        # Getting the difference betweeen coordinates. Should be 1 or -1 for one coordinate only.
        # The other coordinates should become 0
        res = np.abs(basin1_coords - basin2_coords)
        
        # computing the surface
        surface = res[:,0] * self.ttst.y_length * self.ttst.z_length + \
                  self.ttst.x_length * res[:,1] * self.ttst.z_length + \
                  self.ttst.x_length * self.ttst.y_length * res[:,1]
        return np.sum(surface)

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

        # TODO: replace with real values
        self.x_length = 0.1
        self.y_length = 0.2
        self.z_length = 0.3

        # tolerance to merge basins
        self.basins_tolerance = self.step

        # create flat index array
        self.x_size, self.y_size, self.z_size = data.shape
        xindices = np.arange(self.x_size)
        yindices = np.arange(self.y_size)
        zindices = np.arange(self.z_size)
        coordinates = np.stack(np.meshgrid(xindices, yindices, zindices)).transpose((2,1,3,0))
        coordinates = coordinates.reshape(self.x_size*self.y_size*self.z_size, 3)

        # create flat Energy array
        flat_values = data.flatten()

        # sort energy array and index array according to the energy values
        self.forward_indx_permutation = np.argsort(flat_values) # order indices by energy
        self.reverse_indx_permutation = np.argsort(self.forward_indx_permutation) # create back permutation array

        # make the array of energies and array of coordinates sorted by energy.
        self.flat_values = flat_values[self.forward_indx_permutation] # apply new ordering to the energy array
        self.sorted_coordinates = coordinates[self.forward_indx_permutation] # apply the same ordering to the

        # coordinates array
        self.basins = [] # [[basin 0 point flat indices], [basin 1 point flat indices], ...]
        self.ts = [] # Transition States

    def _chunk_energy_and_indices(self):
        index_min=0
        for i, energy in enumerate(self.flat_values):
            if energy - self.flat_values[index_min] > self.step:
                yield self.flat_values[index_min], np.array(range(index_min, i))
                index_min = i


    def find_basins(self):
        """TODO: The function is still under development"""
        # flat_indices reference to all the energies in the range [e_current; e_current+step]
        total_len = 0
        for e_current, flat_indices in self._chunk_energy_and_indices():
            n_basins = len(self.basins)
            n_points = len(flat_indices)
            total_len += n_points
            print("n_basins", len(self.basins), "n_points", n_points, "total_points", total_len)  
            #if (total_len > 40000):
            #    sys.exit(0)

            # (I) add points to the existing basins until it's possible
            while True: # infinite loop
                added_elements = False
                for bsn in self.basins:
                    neigh_indxs, flat_indices = bsn.choose_neighbors_from_list(flat_indices)
                    if neigh_indxs.size > 0:
                        bsn.append(neigh_indxs)
                        added_elements = True
                    if flat_indices.size == 0:
                        break
                if not added_elements or flat_indices.size == 0:
                    break

            # (II) form clusters from the leftover elements in the flat_indices array
            # and add them as new basins
            for cluster in self.form_clusters(flat_indices):
                self.basins.append(Basin(indices=cluster, ttst=self,tolerance=self.basins_tolerance))

            # (III) merge too small basins into the big ones
            to_delete = []
            for n_bsn1, bsn1 in enumerate(self.basins):
                if not bsn1.is_definitely_basin: # if the difference between smallest and largest values is too small
                    for n_bsn2 in range(n_bsn1+1, len(self.basins)): # loop over the remaining basins
                        bsn2 = self.basins[n_bsn2]
                        if bsn1.is_neighbor_with(bsn2) and bsn1.close_by_energy(bsn2): # merge only if two basins are
                        # neighbors and are close energetically
                            bsn2.merge_basin(bsn1)
                            to_delete.append(bsn1)
                            break

            # (IV) delete basins that were merged in the previous step
            for bsn in to_delete:
                print("removing basin: ", bsn)
                self.basins.remove(bsn)

        # TODO: delete, when not needed for testing
        tot = 0
        for i, b in enumerate(self.basins):
            print("basin_{} {} elements".format(i, len(b.indices)))
            tot += len(b.indices)
        print(tot)

    def find_transition_states(self):
        for n_bsn1, bsn1 in enumerate(self.basins):
            for n_bsn2 in range(n_bsn1+1, len(self.basins)): # loop over the other basins
                bsn2 = self.basins[n_bsn2]
                if bsn1.is_neighbor_with(bsn2):
                    self.ts.append(TransitionState(bsn1, bsn2, self))



    def form_clusters(self, indices):
        """The function gets a list of point indices and clusters them if they are neighbors.

        :param indices: list of indices
        :type indices: list
        """
        
        # neighbor operations
        # first element: 0 means 'x', 1, means 'y', 2 means 'z'
        # second element: add or subtract 1.
        neighbor_operations = [(0,  1), (0, -1), (1,  1), (1, -1), (2,  1), (2, -1)]

        # list of clusters
        clusters = []
        for point, point_index in zip(self.sorted_coordinates[indices], indices):
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
        """Converts list of coordinates [[x1, y1, z1], ...] to the list of indices[i1, i2, ...] in the sorted array

        :param coords: a list of point coordinates
        :type: 2D array, [[x1, y1, z1], ...]

        :returns: a list of indices of points in sorted array
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
        indxs = coords[:, 0] * self.y_size * self.z_size + coords[:, 1] * self.z_size + coords[:, 2]
        return self.reverse_indx_permutation[indxs] # indices in the sorted array
