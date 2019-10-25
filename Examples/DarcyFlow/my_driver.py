# Darcy flow in a triangle mesh of a planar domain, with constant
# velocity prescribed on boundary. 

# Reference:

# Numerical method for Darcy flow derived using Discrete Exterior Calculus
# A. N. Hirani, K. B. Nakshatrala, J. H. Chaudhry
# See arXiv:0810.3434v3 [math.NA] on http://arxiv.org/abs/0810.3434

from numpy import mat, zeros, sort, asarray, loadtxt, array, dot, \
     concatenate, sign, vstack, argmax, nonzero
import numpy as np
from numpy.linalg import norm, det
from scipy.sparse import bmat
from scipy.sparse.linalg import spsolve
from matplotlib.pylab import figure, gca, triplot, show
from decagt import SimplicialComplex, DiscreteExteriorCalculus #, simplex_quivers, signed_volume
import numpy as np
import scipy


# Input parameters
velocity = np.array([1, 0])

# Read the mesh
vertices = np.loadtxt('vertices.txt')
triangles = np.loadtxt('triangles.txt', dtype='int') - 1

# Make a simplicial complex from it
sc = SimplicialComplex(vertices,triangles)
sc.build_complex()
# Nk is number of k-simplices
N1 = sc.num_simplices[1]
N2 = sc.num_simplices[2]
# Permeability is k > 0 and viscosity is mu > 0
k = 1; mu = 1
# The matrix for the full linear system for Darcy in 2D, not taking into
# account the boundary conditions, is :
# [-(mu/k)star1 d1^T ]
# [    d1          Z ] 
# where Z is a zero matrix of size N2 by N2. 
# The block sizes are 
#   N1 X N1    N1 X N2
#   N2 X N1    N2 X N2

dec = DiscreteExteriorCalculus(sc)
dec.compute_hodge_stars()

d1 = dec.boundary_matrices[1].T; star1 = dec.hodge_stars[1]

# print(np.array(dec.simplices[1]))

A = bmat([[(-mu/k)*star1, d1.T],
          [d1, None]], format='csr')
b = zeros(N1 + N2) # RHS vector

all_fluxes = zeros(N1)
# Find boundary and internal edges
boundary_edges = np.array(sc.simplices[1])[
    (d1.T * np.ones(sc.num_simplices[2])).nonzero()[0]]
boundary_edges = np.array(boundary_edges)
b1 = sc.boundary_matrices[1]

boundary_indices = np.dot(b1.toarray(), np.ones(N2)).nonzero()[0]

num_boundary_edges = len(boundary_indices)
# internal_edges = set(sc[1].simplex_to_index.keys()) - set(boundary_edges)
internal_indices = list(set(list(range(N1))) - set(boundary_indices))
num_internal_edges = sc.num_simplices[1] - num_boundary_edges
# Assume triangles oriented the same way so can look at any triangle
s = sign(det(vertices[triangles[0,1:]] - vertices[triangles[0,0]]))
for i, e in enumerate(boundary_edges):
    evector = (vertices[e[1]] - vertices[e[0]])
    normal = array([-evector[1], evector[0]])
    all_fluxes[boundary_indices[i]] = -s * (1/norm(evector)**2) * \
                              dot(velocity, normal) * \
                              abs(det(vstack((normal, evector))))
pressures = zeros(N2)
# Adjust RHS for known fluxes and pressures
b = b - A * concatenate((all_fluxes,pressures))

# #Remove entries of b corresponding to boundary fluxes and known pressure
# #Pressure at the right most triangle circumcenter is assumed known
# pressure_indices = list(range(N1, N1+N2))
# pressure_indices.remove(N1 + argmax(sc.circumcenter[:,0]))
# entries_to_keep = concatenate((internal_indices, pressure_indices))
# b = b[entries_to_keep]
# # Remove corresponding rows and columns of A
# A = A[entries_to_keep][:,entries_to_keep]
# u = spsolve(A,b)

# print("u is")
# print(u)
