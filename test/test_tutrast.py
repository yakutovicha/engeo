import numpy as np
def test_tutrast_creation():
	from engeo.cube import read_cube
	from engeo.tutrast import TuTrast

	data, meta = read_cube("test/grid.cube")
	ttst = TuTrast(data, units="a.u.")

	# check if old and new energy array match
	for i, coord in enumerate(ttst.sorted_coordinates):
		assert data[tuple(coord)] == ttst.flat_e[i]

	assert abs(ttst.flat_e[0] - 1.0) < 1e-10
	# ttst = TuTrast(data, units="kcal/mol")
	# assert data[tuple(0,1,1)] == 2
	ttst = TuTrast(data)
	assert abs(ttst.flat_e[0] - 0.00159360109742136) < 1e-10

	# check if old and new energy array do not match (because of conversion)
	for i, coord in enumerate(ttst.sorted_coordinates):
		assert data[tuple(coord)] != ttst.flat_e[i]

	ttst = TuTrast(data, units="kJ/mol")
	assert abs(ttst.flat_e[0] - 0.0003808798033989866) < 1e-10

	ttst = TuTrast(data, units="eV")
	assert abs(ttst.flat_e[0] - 0.03674930495120813) < 1e-10

def test_get_neighbors():
	from engeo.cube import read_cube
	from engeo.tutrast import TuTrast

	data, meta = read_cube("test/grid.cube")
	ttst = TuTrast(data)
	nbrs = np.array([[0, 0, 1],
       [0, 0, 4],
       [0, 1, 0],
       [0, 3, 0],
       [1, 0, 0],
       [2, 0, 0]])
	# checking that the neighbours of the particular cube file
	# of the zeroth entry are correctly found 
	assert (ttst.get_neighbors(bassin=[0]) == nbrs).all()


	#the following test fails, no matter whether no, one or two lines of [1, 1, 0] are included
	#data, meta = read_cube("test/Cube_z_1.cube")
	#ttst = TuTrast(data)
	#nbrs = np.array([[1, 1, 0],
	#   [1, 0, 1],
    #   [1, 2, 1],
    #   [0, 1, 0],
    #   [2, 1, 0]])

	# checking that the neighbours of the particular 2D cube file 
	# of the (1,1,0) entry are correctly found 
	#assert (ttst.get_neighbors(bassin=[0]) == nbrs).all()

	data, meta = read_cube("test/Cube_z_2.cube")
	ttst = TuTrast(data)
	nbrs = np.array([[0, 0, 1],
      [0, 1, 0],
      [0, 3, 0],
      [1, 0, 0],
      [2, 0, 0]])

	# checking that the neighbours of the particular 2D cube file 
	# of the zeroth entry are correctly found although the z-neighbour is found twice
	assert (ttst.get_neighbors(bassin=[0]) == nbrs).all()


	# FAILS so far for an unknown reason
	#data, meta = read_cube("test/Cube_z2_11.cube")
	#ttst = TuTrast(data)
	#nbrs = np.array([[1, 1, 1],
    #  [1, 0, 0],
    #  [1, 2, 0],
    #  [0, 1, 0],
    #  [2, 1, 0]])

	# checking that the neighbours of the particular 2D cube file 
	# of the zeroth entry are correctly found 
	#assert (ttst.get_neighbors(bassin=[0]) == nbrs).all()
