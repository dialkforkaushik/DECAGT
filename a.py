import decagt
import numpy as np
import time
import matplotlib.pyplot as plt

y = []
x = [1.0, 0.5, 0.25, 0.125] #max edge length for 3d cube meshes (hardcoded)
q_order = 2
for j in range(4):
	U = []
	vertices = np.loadtxt('tests/meshes/3d/cube_hierarchy/vertices' + str(j)+'.txt')
	triangles = np.loadtxt('tests/meshes/3d/cube_hierarchy/tets' + str(j) + '.txt', dtype='int')

	sc = decagt.SimplicialComplex(vertices, triangles)
	sc.build_complex()

	####################### Code to find max edge length #################################
	# pts = []
	# x = []
	# max_vol = -1
	# for i in sc.simplices[1]:
	# 	a = np.array(sc.vertices[i[0]]) - np.array(sc.vertices[i[1]])
	# 	vol = np.linalg.norm(a, 2)
	# 	if vol > max_vol:
	# 		max_vol = vol
	# x.append(max_vol)
	#######################################################################################

	for i in range(len(sc.vertices)):
		U.append(decagt.get_analytical_soln(sc.vertices[i]))
	
	error = 0.0
	error = decagt.error_0(U,q_order,sc.simplices,sc.vertices,sc.num_simplices)
	y.append(error)

	print("For mesh " + str(j+1))
	print("error: " + str(error))
	print("max edge len: " + str(x[j]))
	print("")

plt.loglog(x,y)
plt.ylabel("error")
plt.xlabel("Max Edge length")
plt.title("q_order = " + str(q_order))
plt.show()
