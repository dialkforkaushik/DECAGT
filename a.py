import decagt
import numpy as np
import time
import matplotlib.pyplot as plt

y = []
# x = []
x = [1.0, 0.5, 0.25, 0.125] #max edge length for 3d cube meshes (hardcoded for faster computation)
q_order = 3
for j in range(4):
	U = []
	vertices = np.loadtxt('tests/meshes/3d/cube_hierarchy/vertices' + str(j)+'.txt')
	triangles = np.loadtxt('tests/meshes/3d/cube_hierarchy/tets' + str(j) + '.txt', dtype='int')

	sc = decagt.SimplicialComplex(vertices, triangles)
	sc.build_complex()

	####################### Code to find max edge length #################################
	# pts = []
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
	error = decagt.quadratic_error_0_bb_mass(U,q_order,sc.simplices,sc.vertices,sc.num_simplices)
	y.append(error)

	print("For mesh " + str(j+1))
	print("error: " + str(error))
	print("max edge len: " + str(x[j]))
	print("")

x = np.array(x)
y = np.array(y)

fig = plt.figure()
plt.loglog(x,y)
plt.ylabel("error")
plt.xlabel("Max Edge length")
plt.title("q_order = " + str(q_order))
fig.savefig('fig.png')

slope_ary = []
for i in range(1, y.shape[0]):
	slope, intercept = np.polyfit(np.log(np.sort(x[i-1:i+1])), np.log(np.sort(y[i-1:i+1])), 1)
	slope_ary.append(slope)
print("slope: ", slope_ary)
