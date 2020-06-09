import decagt
import numpy as np
import time
import matplotlib.pyplot as plt

# Function to plot and find slope
def plot(x, y, ylabel="", xlabel="", title="", figname="fig.png"):
	x = np.array(x)
	y = np.array(y)

	fig = plt.figure()
	plt.loglog(x,y)
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.title(title)
	fig.savefig(figname)

	slope_ary = []
	for i in range(1, y.shape[0]):
		slope, intercept = np.polyfit(np.log(x[i-1:i+1]), np.log(y[i-1:i+1]), 1)
		slope_ary.append(slope)
	print("slope: ", slope_ary)

	plt.show()


def test1():
	# List to store error values
	errors = []

	#max edge length for 3d cube meshes (hardcoded for faster computation)
	edge_len = [1.0, 0.5, 0.25, 0.125]

	# vertices = np.loadtxt('tests/meshes/3d/cube_hierarchy/vertices0.txt')
	# tets = np.loadtxt('tests/meshes/3d/cube_hierarchy/tets0.txt', dtype='int')

	verts = np.loadtxt('tests/meshes/testMesh2/vertices.txt')
	tets = np.loadtxt('tests/meshes/testMesh2/tets.txt', dtype='int')
	tets = [tets]

	# Choosing only the 25th tetrahedron
	# verts = []
	# element = tets[25]
	# tets = [[0, 1, 2, 3]]
	# for i in element:
	# 	verts.append(vertices[i])


	# Creating DECAGT objects
	sc = decagt.SimplicialComplex(verts, tets)
	fem = decagt.FiniteElementExteriorCalculus(sc)

	# Looping over Bernstein_Bezier Polynomial degrees (n = degree)
	for n in range(1, 21):
		
		error = 0.0
		error = fem.bb_error_1(n, fem.simplices, fem.vertices, fem.num_simplices, 4)
		# errors.append(error)

		# print("For Polynomial of Degree " + str(n))
		print("error: " + str(error))

	return np.arange(1, 21, 1), errors


def test2(degree_analytical):

	#max edge length for 3d cube meshes (hardcoded for faster computation)
	max_edge_len = [1.0, 0.5, 0.25]
	
	arr = list()
	if degree_analytical == 0:
		arr = [1, 2, 3, 4, 5, 10]
		# arr = [4, 5, 10]
	else:
		arr = np.arange(1, degree_analytical + 1, 1)

	for n in arr:
		errors = []

		for j in range(3):
			vertices = np.loadtxt('tests/meshes/3d/cube_hierarchy/vertices' + str(j)+'.txt')
			triangles = np.loadtxt('tests/meshes/3d/cube_hierarchy/tets' + str(j) + '.txt', dtype='int')

			sc = decagt.SimplicialComplex(vertices, triangles)
			sc.build_complex()
			fem = decagt.FiniteElementExteriorCalculus(sc)
			
			error = 0.0
			error = fem.bb_error(n, sc.simplices, sc.vertices, sc.num_simplices, 8)
			errors.append(error)

			print("For mesh " + str(j+1))
			print("error: " + str(error))
			print("max edge len: " + str(max_edge_len[j]))
			print("")

		plot(max_edge_len, errors, "Error", "Max Edge Length", "Test 2", str(n)+"_"+str(j))

	return max_edge_len, errors


def test3():


	vertices = np.loadtxt('tests/meshes/3d/cube_hierarchy/vertices0.txt')
	tets = np.loadtxt('tests/meshes/3d/cube_hierarchy/tets0.txt', dtype='int')

	# Choosing only the 25th tetrahedron
	verts = []
	element = tets[25]
	tets = [[0, 1, 2, 3]]
	for i in element:
		verts.append(vertices[i])


	# Creating DECAGT objects
	sc = decagt.SimplicialComplex(verts, tets)
	fem = decagt.FiniteElementExteriorCalculus(sc)

	# Looping over Bernstein_Bezier Polynomial degrees (n = degree)
	for n in [1, 2, 5, 10]:
		# List to store error values
		errors = []

		for q in range(1, 11):
		
			# error = 0.0
			error = fem.bb_error(n, fem.simplices, fem.vertices, fem.num_simplices, q)
			errors.append(error)

			print("Polynomial Degree:", n)
			print("For Quadrature Order:", q)
			print("error: " + str(error))

		plot(np.arange(1, 11, 1), errors, "Error", "Quadrature Order", "Test 3")

	return np.arange(1, 11, 1), errors


x, y = test1()
# plot(x, y, "error", "Degree of Polynomial", "Test 1")

# x, y = test2(0)

# x, y = test3()

