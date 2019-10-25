"""
Reads ascii vertex and element files, writes a pydec mesh and displays it
"""

import scipy
from decagt import SimplicialComplex
from matplotlib.pyplot import triplot, show
import numpy as np

vertices = scipy.loadtxt("v.txt")
elements = scipy.loadtxt("s.txt",dtype='int32') - 1
mymesh = SimplicialComplex(vertices,elements)
mymesh.build_complex()

v = np.array(mymesh.vertices)

triplot(v[:,0], v[:,1], triangles=mymesh.simplex)
show()

