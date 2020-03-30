"""
Solve the resonant cavity problem with Whitney forms.

References:
    Douglas N. Arnold and Richard S. Falk and Ragnar Winther
    "Finite element exterior calculus: from Hodge theory to numerical
    stability"
    Bull. Amer. Math. Soc. (N.S.), vol. 47, No. 2, pp. 281--354
    DOI : 10.1090/S0273-0979-10-01278-4

"""
from decagt import SimplicialComplex, FiniteElementExteriorCalculus
from numpy import loadtxt
from scipy import real, zeros
from scipy.linalg import eig
from matplotlib.pylab import quiver, figure, triplot, show
import numpy as np

# Read in mesh data from files and construct complex
vertices = loadtxt('vertices.txt', dtype=float)
triangles = loadtxt('triangles.txt', dtype=int)
sc = SimplicialComplex(vertices,triangles)
fem = FiniteElementExteriorCalculus(sc)
fem.compute_hodge_stars()

# Construct stiffness and mass matrices 
K = fem.boundary_matrices[1] * fem.hodge_stars[2] * fem.boundary_matrices[1].T
M = fem.hodge_stars[1];


# Eliminate Boundaries from matrices
# boundary_edges = sc.boundary()
# non_boundary_edges = set(sc[1].simplex_to_index.keys()) - set(boundary_edges)
b1 = fem.boundary_matrices[1]
N1 = fem.num_simplices[1];
N2 = fem.num_simplices[2];
non_boundary_indices = np.dot(b1.toarray(), np.ones(N2)).nonzero()[0]

non_boundary_indices = np.setdiff1d(np.array(range(N1)), non_boundary_indices)


# Eliminate boundary conditions
K = K[non_boundary_indices,:][:,non_boundary_indices]
M = M[non_boundary_indices,:][:,non_boundary_indices]

# Compute eigenvalues and eigenvectors
# (could use sparse eigenvalue solver instead)
eigenvalues, eigenvectors = eig(K.todense(), M.todense())

# Plot eigenvalues
# NUM_EIGS = 50 # Number of eigenvalues to plot
# values = sorted([x for x in real(eigenvalues) if x > 1e-10])[0:NUM_EIGS]
# ax = figure().gca()
# ax.set_title('First ' + str(len(values)) + ' Eigenvalues\n\n')
# ax.hold(True)
# ax.plot(values,'ko')


# # Plot the eigenvector 1-cochain as a vector field
# N = 2 # Which non-zero eigenvector to plot?
# non_zero_values = real(eigenvectors[:,list(eigenvalues).index(values[N])])
# all_values = zeros((sc[1].num_simplices,))
# all_values[non_boundary_indices] = non_zero_values
# bases, arrows = simplex_quivers(sc,all_values)
# ax = figure().gca()
# ax.set_title('Mode #' + str(N+1))
# ax.quiver(bases[:,0],bases[:,1],arrows[:,0],arrows[:,1])
# ax.triplot(sc.vertices[:,0], sc.vertices[:,1], sc.simplices)
# ax.axis('equal')

# show()

