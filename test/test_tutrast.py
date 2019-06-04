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