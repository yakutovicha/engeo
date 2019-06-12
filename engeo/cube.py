from __future__ import print_function

import numpy as np
def _getline(cube):
    """
    Read a line from cube file where first field is an int 
    and the remaining fields are floats.
    
    params:
        cube: file object of the cube file
    
    returns: (int, list<float>)
    """
    l = cube.readline().strip().split()
    return int(l[0]), map(float, l[1:])

def read_cube(fname):
    """
    Read cube file into numpy array

    params:
        fname: filename of cube file

    returns: (data: np.array, metadata: dict)
    """
    meta = {}
    with open(fname, 'r') as cube:
        cube.readline(); cube.readline() # ignore comments 
        natm, meta['org'] = _getline(cube)
        nx, meta['xvec'] = _getline(cube)
        ny, meta['yvec'] = _getline(cube)
        nz, meta['zvec'] = _getline(cube)
        meta['atoms'] = [_getline(cube) for i in range(natm)]
        data = np.zeros((nx*ny*nz))
        idx = 0
        for line in cube:
            for val in line.strip().split():
                data[idx] = float(val)
                idx += 1
    data = np.reshape(data, (nx, ny, nz))
    return data, meta

