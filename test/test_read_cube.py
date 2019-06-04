def test_read_cube_file():
	from engeo.cube import read_cube
	data, meta = read_cube("test/grid.cube")

	assert data.shape == (3, 4, 5)
	assert abs(data[0, 0, 0] - 1.0) < 1e-10
