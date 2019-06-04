def test_tutrast_creation():
	from engeo.cube import read_cube
	from engeo.tutrast import TuTrast

	data, meta = read_cube("test/grid.cube")
	ttst = TuTrast(data)

	# check if old and new energy array match
	for i, coord in enumerate(ttst.sorted_coordinates):
		assert data[tuple(coord)] == ttst.flat_e[i]
