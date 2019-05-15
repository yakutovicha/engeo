from __future__ import print_function
import sys
from engeo.cube import read_cube

if len(sys.argv) != 2:
    print("Usage: engeo.py <cube_file>")
    sys.exit(1)

data, meta = read_cube(sys.argv[1])
