import numpy as np
def test_tutrast_creation():
    from engeo.cube import read_cube
    from engeo.tutrast import TuTraSt

    data, meta = read_cube("test/grid.cube")
    ttst = TuTraSt(data, units="a.u.")

    # check if old and new energy array match
    for i, coord in enumerate(ttst.sorted_coordinates):
        assert data[tuple(coord)] == ttst.flat_e[i]

    assert abs(ttst.flat_e[0] - 1.0) < 1e-10
    # ttst = TuTraSt(data, units="kcal/mol")
    # assert data[tuple(0,1,1)] == 2
    ttst = TuTraSt(data)
    assert abs(ttst.flat_e[0] - 0.00159360109742136) < 1e-10

    # check if old and new energy array do not match (because of conversion)
    for i, coord in enumerate(ttst.sorted_coordinates):
        assert data[tuple(coord)] != ttst.flat_e[i]

    ttst = TuTraSt(data, units="kJ/mol")
    assert abs(ttst.flat_e[0] - 0.0003808798033989866) < 1e-10

    ttst = TuTraSt(data, units="eV")
    assert abs(ttst.flat_e[0] - 0.03674930495120813) < 1e-10

def test_coordinates_to_indices_and_back():
    from engeo.cube import read_cube
    from engeo.tutrast import TuTraSt
    data, meta = read_cube("test/grid.cube")
    ttst = TuTraSt(data)
    coords = np.array([
        # x=0
        [0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 4],
        [0, 1, 0], [0, 1, 1], [0, 1, 2], [0, 1, 3], [0, 1, 4],
        [0, 2, 0], [0, 2, 1], [0, 2, 2], [0, 2, 3], [0, 2, 4],
        [0, 3, 0], [0, 3, 1], [0, 3, 2], [0, 3, 3], [0, 3, 4],
        # x=1
        [1, 0, 0], [1, 0, 1], [1, 0, 2], [1, 0, 3], [1, 0, 4],
        [1, 1, 0], [1, 1, 1], [1, 1, 2], [1, 1, 3], [1, 1, 4],
        [1, 2, 0], [1, 2, 1], [1, 2, 2], [1, 2, 3], [1, 2, 4],
        [1, 3, 0], [1, 3, 1], [1, 3, 2], [1, 3, 3], [1, 3, 4],
        # x=2
        [2, 0, 0], [2, 0, 1], [2, 0, 2], [2, 0, 3], [2, 0, 4],
        [2, 1, 0], [2, 1, 1], [2, 1, 2], [2, 1, 3], [2, 1, 4],
        [2, 2, 0], [2, 2, 1], [2, 2, 2], [2, 2, 3], [2, 2, 4],
        [2, 3, 0], [2, 3, 1], [2, 3, 2], [2, 3, 3], [2, 3, 4],
    ])
    indxs = ttst.coordinates_to_indices(coords)
    coords2 = ttst.sorted_coordinates[indxs]
    assert (coords == coords2).all()

def test_get_neighbors():
    from engeo.cube import read_cube
    from engeo.tutrast import TuTraSt, Basin

    data, meta = read_cube("test/grid.cube")
    ttst = TuTraSt(data)
    nbrs = np.array([[0, 0, 1],
       [0, 0, 4],
       [0, 1, 0],
       [0, 3, 0],
       [1, 0, 0],
       [2, 0, 0]])
    # checking that the neighbours of the particular cube file
    # of the zeroth entry are correctly found

    b = Basin(indeces=[0], shape=data.shape, sorted_coordinates=ttst.sorted_coordinates, coordinates_to_indices=ttst.coordinates_to_indices)
    print(ttst.sorted_coordinates[0])
    assert len(b.get_neighbors_indeces()) == 6

    for elem in ttst.sorted_coordinates[list(b.get_neighbors_indeces())]:
        assert elem in nbrs


    #the following test fails, no matter whether no, one or two lines of [1, 1, 0] are included
    #data, meta = read_cube("test/Cube_z_1.cube")
    #ttst = TuTraSt(data)
    #nbrs = np.array([[1, 1, 0],
    #   [1, 0, 1],
    #   [1, 2, 1],
    #   [0, 1, 0],
    #   [2, 1, 0]])

    # checking that the neighbours of the particular 2D cube file 
    # of the (1,1,0) entry are correctly found 
    #assert (ttst.get_neighbors(bassin=[0]) == nbrs).all()

    # data, meta = read_cube("test/Cube_z_2.cube")
    # ttst = TuTraSt(data)
    # nbrs = np.array([[0, 0, 1],
    #   [0, 1, 0],
    #   [0, 3, 0],
    #   [1, 0, 0],
    #   [2, 0, 0]])

    # checking that the neighbours of the particular 2D cube file 
    # of the zeroth entry are correctly found although the z-neighbour is found twice
    # assert (ttst.get_neighbors(bassin=[0]) == nbrs).all()


    # FAILS so far for an unknown reason
    #data, meta = read_cube("test/Cube_z2_11.cube")
    #ttst = TuTraSt(data)
    #nbrs = np.array([[1, 1, 1],
    #  [1, 0, 0],
    #  [1, 2, 0],
    #  [0, 1, 0],
    #  [2, 1, 0]])

    # checking that the neighbours of the particular 2D cube file 
    # of the zeroth entry are correctly found 
    #assert (ttst.get_neighbors(bassin=[0]) == nbrs).all()
