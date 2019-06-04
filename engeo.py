from __future__ import print_function
import sys
from engeo.cube import read_cube
from engeo.tutrast import TuTrast

if len(sys.argv) != 2:
    print("Usage: engeo.py <cube_file>")
    sys.exit(1)

# reading the cube file
data, meta = read_cube(sys.argv[1])

# initializing the TuTrast object
ttst = TuTrast(data)

# find basins
ttst.find_bassins()
