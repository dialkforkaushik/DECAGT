#include "definitions.h"
#include "simplicial_complex.h"
#include "core_utils.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <Eigen/Eigen>
#include <algorithm>
#include <set>
#include <cmath>
#include <limits>
#include <dlfcn.h>
#include <stdlib.h>

#ifdef PYTHON
	#include <pybind11/pybind11.h>
	#include <pybind11/stl.h>
	#include <pybind11/numpy.h>
#endif

int factorial(int &n) {
	if (n == 0) {
		n = 1;
	}
	else if (n > 0) {
		for (int i = n - 1; i > 0; --i) {
	        n = n * i;   
	    }
	}
	else {
		return FAILURE;
	}
    
    return SUCCESS;
}

int l2_norm(double norm,
			VectorD v) {

	norm = 0.0;
	size_t size = v.size();
	for(size_t i = 0; i < size; ++i) {
		norm += v[i]*v[i];
	}
	norm = sqrt(norm);

	return SUCCESS;
}


double get_simplex_volume(Vector2D &vertices) {

	size_t complex_dimension = vertices.size() - 1;
	size_t embed_dim = vertices[0].size();

	DenMatD B = DenMatD::Ones(complex_dimension + 2, complex_dimension + 2);
	for (size_t i = 0; i < complex_dimension + 1; ++i) {
		for (size_t j = 0; j < complex_dimension + 1; ++j) {
			EigVectorD v(embed_dim);
			for (size_t k = 0; k < embed_dim; ++k) {
				v.coeffRef(k) = vertices[i][k] - vertices[j][k];
			}
			B.coeffRef(i+1, j+1) = pow(v.norm(), 2);
		}	
	}
	B.coeffRef(0, 0) = 0.0;

	int fac = complex_dimension;
	factorial(fac);

	double num = pow(-1, complex_dimension + 1) * B.determinant();
	double den = pow(2, complex_dimension) * pow(fac, 2);
	double vol = sqrt(num/den);

	return vol;
}


int unsigned_volume(Vector2D &pts,
					double &vol) {

    int row = pts.size();
    int col = pts[0].size();

    DenMatD A (row - 1, col);
    MapEigVectorD B (pts[0].data(), col);

    for (int i = 1; i < row; ++i) {
    	MapEigVectorD C(pts[i].data(), col);
        for (int j = 0; j < col; ++j) {
            A.row(i-1) = C - B;
        }
    }

    int fact = row - 1;
    factorial(fact);
    vol = (A * A.transpose()).determinant();
    vol = (std::sqrt(std::abs(vol)))/fact;

    return SUCCESS;
}


int signed_volume(Vector2D &pts,
				  double &vol) {

    int row = pts.size();
    int col = pts[0].size();

    DenMatD A (row - 1, col);
    MapEigVectorD B (pts[0].data(), col);

    for (int i = 1; i < row; ++i) {
    	MapEigVectorD C(pts[i].data(), col);
        for (int j = 0; j < col; ++j) {
            A.row(i-1) = C - B;
        }
    }

    int fact = row - 1;
    factorial(fact);
    vol = A.determinant()/fact;

    return SUCCESS;
}

inline DenMatD circumcenter_barycentric(DenMatD &pts_matrix, DenMatD &bary_coords) {
	int row = pts_matrix.rows();
	int col = pts_matrix.cols();

    DenMatD mat1 = 2 * (pts_matrix * pts_matrix.transpose());
    DenMatD mat2 = DenMatD::Ones(row, 1);
    DenMatD mat3 = DenMatD::Ones(1, row);
    DenMatD mat4 = DenMatD::Zero(1, 1);

    int A_rows = mat1.rows() + mat3.rows();
    int A_cols = mat1.cols() + mat2.cols();
    DenMatD A (A_rows, A_cols);
    A.topLeftCorner(mat1.rows(), mat1.cols()) = mat1;
	A.topRightCorner(mat2.rows(), mat2.cols()) = mat2;
	A.bottomLeftCorner(mat3.rows(), mat3.cols()) = mat3;
	A.bottomRightCorner(mat4.rows(), mat4.cols()) = mat4;

	mat1 = pts_matrix;
	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {
			mat1.coeffRef(i, j) = mat1.coeffRef(i, j) * mat1.coeffRef(i, j);
		}
	}
	mat2.resize(1, mat1.rows());
	mat3 = DenMatD::Ones(1, 1);
	for (int i = 0; i < mat2.cols(); ++i) {
		mat2.coeffRef(0, i) = mat1.row(i).sum();
	}
	
	DenMatD b (mat2.cols() + 1, 1);
	b.topRows(mat2.cols()) = mat2.transpose();
	b.bottomRows(1) = mat3.transpose();
	DenMatD x = A.colPivHouseholderQr().solve(b);
	bary_coords = x.topRows(x.rows() - 1);

	return bary_coords;
}


inline int get_circumcenter(VectorD &center,
				     		double &radius,
				     		Vector2D &pts) {

	int row = pts.size();
    int col = pts[0].size();

    DenMatD pts_matrix (row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            pts_matrix.coeffRef(i, j) = pts[i][j];
        }
    }

	DenMatD bary_coords;
	circumcenter_barycentric(pts_matrix, bary_coords);
	DenMatD temp_matrix (0, 0);
	temp_matrix = bary_coords.transpose() * pts_matrix;
	for (int i = 0; i < temp_matrix.cols(); ++i) {
		center.push_back(temp_matrix.coeffRef(0, i));
	}
	
	temp_matrix.resize(1, col);
	for (int i = 0; i < col; ++i) {
		temp_matrix.coeffRef(0, i) = pts_matrix.coeffRef(0, i) - center[i];
	}
	radius = temp_matrix.norm();

    return  SUCCESS;
}


int calculate_dual_volume(Vector2D &dual_volume,
						Vector2D &temp_centers,
						Vector2D &vertices,
						Vector3I &simplices,
						VectorMap3I &adjacency1d,
						VectorI &pts,
						VectorD &highest_dim_circumcenter,
						int dim,
						int index,
						size_t complex_dimension) {

	VectorD center;
	double radius;
	double vol;

	Vector2D vertex_pts;

	
	if (dim == complex_dimension) {
		temp_centers.push_back(highest_dim_circumcenter);
	}
	else {
		for (int i = 0; i < dim + 1; ++i) {
			vertex_pts.push_back(vertices[pts[i]]);
		}
		get_circumcenter(center,
					radius,
					vertex_pts);

		temp_centers.push_back(center);
	}	

	unsigned_volume(temp_centers, vol);

	#ifdef MULTICORE
		#pragma omp critical
	#endif

	dual_volume[dim][index] += vol;
	
	if (dim == 0) {
		return SUCCESS;
	}

	for (auto it = adjacency1d[dim][index][0].begin(); it != adjacency1d[dim][index][0].end(); ++it) {
		calculate_dual_volume(dual_volume,
							  temp_centers,
							  vertices,
							  simplices,
							  adjacency1d,
							  simplices[dim - 1][it->first],
							  highest_dim_circumcenter,
							  dim - 1,
							  it->first,
							  complex_dimension);

		temp_centers.pop_back();
	}

	return SUCCESS;
} 


int count_columns(std::string &line,
				  int &columns) {

    std::stringstream s;
    s << line;
    
    columns = 0;    
    double value;
    
    while(s >> value)
        columns++;
    
    return SUCCESS;
}


int compute_combinations(Vector2I &vector,
						 int &N,
						 int &K,
                		 Vector3I &combination) {

	std::string bitmask(K, 1);
    bitmask.resize(N, 0);
 
    do {
        Vector2I vec;
        Vector2I sorted_vec;
        for (int i = 0; i < N; ++i)
        {
            if (bitmask[i]) {
            	vec.push_back(vector[i]);
            }
        }
        
        sorted_vec = vec;
        std::sort(sorted_vec.begin(), sorted_vec.end());
        combination.push_back(sorted_vec);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    std::sort(combination.begin(), combination.end());
	combination.erase(std::unique(combination.begin(), combination.end()), combination.end());


	return SUCCESS;
}


int compute_combinations(VectorI &vector,
		         		 int &N,
		         		 int &K,
		         		 Vector2I &combination) {

    std::string bitmask(K, 1);
    bitmask.resize(N, 0);
 
    do {
        VectorI vec;
        VectorI sorted_vec;
        for (int i = 0; i < N; ++i)
        {
            if (bitmask[i]) {
            	vec.push_back(vector[i]);
            }
        }

        sorted_vec = vec;
	    std::sort(sorted_vec.begin(), sorted_vec.end());
        combination.push_back(sorted_vec);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    std::sort(combination.begin(), combination.end());
	combination.erase(std::unique(combination.begin(), combination.end()), combination.end());

	return SUCCESS;
}


int compute_combinations_mesh(Vector2I &simplex,
			         		  int &N,
			         		  int &K,
			         		  Vector2I &combination) {
    
	VectorI vector;
	int n = simplex.size();
    
    for (int i = 0; i < n; ++i) {
    	std::string bitmask(K, 1);
	    bitmask.resize(N, 0);
	    vector.clear();
	    vector = simplex[i];
	    VectorI vec;
	    VectorI sorted_vec;
	 
	    do {
	        vec.clear();
	        for (int i = 0; i < N; ++i)
	        {
	            if (bitmask[i]) {
	            	vec.push_back(vector[i]);
	            }
	        }

	        sorted_vec = vec;
	        std::sort(sorted_vec.begin(), sorted_vec.end());
	        combination.push_back(sorted_vec);
	    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }

    std::sort(combination.begin(), combination.end());
	combination.erase(std::unique(combination.begin(), combination.end()), combination.end());

	return SUCCESS;
}


int get_combinations_mesh(Vector2I &simplex,
						  Vector2I &combination,
						  int k) {

	int n = simplex[0].size();
	if (k == -1) {
		k = n - 1;
	}

	compute_combinations_mesh(simplex, n, k, combination);

	
	return SUCCESS;
}

int get_combinations_simplex(VectorI &simplex,
							 Vector2I &combination,
							 int k) {

	if (k == 0) {
		VectorI vec;
		combination.push_back(vec);

		return SUCCESS;
	}

	int n = simplex.size();
	if (k == -1) {
		k = n - 1;
	}

	compute_combinations(simplex, n, k, combination);
	
	return SUCCESS;
}


int get_combinations_simplex(Vector2I &simplex,
							 Vector3I &combination,
							 int k) {

	int n = simplex.size();

	if (k > n) {
		size_t simplex_size = simplex.size();
		for (int i = 0; i < simplex_size; ++i) {
			combination.push_back(Vector2I());
			combination[i].push_back(simplex[i]);
			}

		return SUCCESS;
		
	}

	if (k == -1) {
		k = n - 1;
	}

	compute_combinations(simplex, n, k, combination);
	
	return SUCCESS;
}


int is_rotated(VectorI v1,
               VectorI v2) { 

    int n = v1.size();

    std::string str_cat = ""; 
    for (int i = 0 ; i < n ; ++i) 
        str_cat = str_cat + "-" + std::to_string(v1[i]); 

    str_cat = str_cat + str_cat;
  
    std::string curr_str = ""; 
    for (int j = 0 ; j < n - 1 ; ++j) 
        curr_str = curr_str + "-" + std::to_string(v2[j]);

    if (str_cat.find(curr_str) == std::string::npos) {
        return FAILURE; 
    }
  
    return SUCCESS; 
}


int k_cochain_norm(EigVectorD U,
				   std::string weight,
				   VectorSpmatD hodge_stars,
				   int k) {
	double norm;
	norm = U.transpose() * hodge_stars[k] * U;
	norm = sqrt(norm);

	return SUCCESS;
}


int get_barycentric(VectorD &point,
					VectorI &simplex,
					Vector2D &vertices,
					EigVectorD &bary_coords) {

    size_t col = simplex.size();
    size_t row = vertices[0].size() + 1;

    DenMatD A = DenMatD::Ones(row, col);
    EigVectorD b = DenMatD::Ones(row, 1);

    for (size_t i = 0; i < col; ++i) {
    	VectorD pts = vertices[simplex[i]];
    	for (size_t j = 0; j < row - 1; ++j) {
    		A.coeffRef(j+1,i) = pts[j];
    	}
    }

    for (size_t i = 0; i < row - 1; ++i) {
    	b.coeffRef(i+1) = point[i];
    }

    bary_coords = A.colPivHouseholderQr().solve(b);

    return SUCCESS;
}


int get_simplices_from_sub_simplex(VectorI &simplex_simplices,
								   VectorMap3I &adjacency1d,
								   int sub_simplex_index,
								   size_t dim,
								   size_t &complex_dimension) {

	if (dim == complex_dimension) {
		simplex_simplices.push_back(sub_simplex_index);
		return SUCCESS;
	}

	for (auto it = adjacency1d[dim][sub_simplex_index][1].begin(); it != adjacency1d[dim][sub_simplex_index][1].end(); ++it) {
		get_simplices_from_sub_simplex(simplex_simplices,
									   adjacency1d,
									   it -> first,
									   dim + 1,
									   complex_dimension);
	}

	simplex_simplices.erase(std::unique(simplex_simplices.begin(), simplex_simplices.end()), simplex_simplices.end());
	return SUCCESS;
}


int get_closest_simplex_to_point(VectorD &point,
								 VectorI &simplex,
								 Vector2D &vertices,
								 Vector3I &simplices,
								 VectorMap3I &adjacency1d,
								 size_t embedding_dimension,
								 size_t complex_dimension) {
	VectorD distances;
	double distance;
	size_t vertices_size = vertices.size();
	for (size_t i = 0; i < vertices_size; ++i) {
		distance = 0;
		for (size_t j = 0; j < embedding_dimension; ++j) {
			distance += pow(vertices[i][j] - point[j], 2);
		}
		distances.push_back(distance);
	}

	auto minimum = std::min_element(distances.begin(), distances.end());
	int vertex = std::distance(distances.begin(), minimum);

	VectorI vertex_simplices;
	get_simplices_from_sub_simplex(vertex_simplices,
								   adjacency1d,
								   vertex,
								   0,
								   complex_dimension);

	size_t vertex_simplices_size = vertex_simplices.size();
	EigVectorD bary_coords;

	#ifdef MULTICORE
		#pragma omp parallel for private(vertex_simplices_size, bary_coords)
	#endif
	for(size_t i = 0; i < vertex_simplices_size; ++i) {
		get_barycentric(point,
						simplices[complex_dimension][vertex_simplices[i]],
						vertices,
						bary_coords);

		size_t bary_coords_size = bary_coords.size();
		int flag = 0;
		for (size_t j = 0; j < bary_coords_size; ++j) {
			if (bary_coords[j] <= 0.0) {
				flag = 1;
			}
		}
		if (flag == 0) {
			#ifdef MULTICORE
				#pragma omp critical
			#endif
			simplex = simplices[complex_dimension][vertex_simplices[i]];
			
			// return SUCCESS;
		}
	}

	return SUCCESS;
}


double get_analytical_soln(VectorD &vec) {

	void* handle = dlopen("./function.so", RTLD_LAZY);
	if (!handle) {
		system("g++ ./function.cc -o ./function.so -shared -fPIC");
		handle = dlopen("./function.so", RTLD_LAZY);
		if (!handle) {
	        std::cerr << "Cannot open file: " << dlerror() << '\n';
	        return FAILURE;
	    }
    }

    dlerror();
    function_t function = (function_t) dlsym(handle, "function");
    const char *dlsym_error = dlerror();
    if (dlsym_error) {
        std::cerr << "Cannot load symbol 'function': " << dlsym_error << '\n';
        dlclose(handle);
        return FAILURE;
    }

    double output = function(vec);
    return output;
}


int read_quadratures(Vector2D &nodes,
					 VectorD &weights,
			  		 std::string data) {

	std::ifstream f1(data);
	std::string l;
	getline(f1, l);
	std::istringstream iss1(l);
	int num;
	iss1 >> num;

	std::fstream file;
	int row = 0;
	int cols;

	file.open(data, std::ios::in);
	if (file.is_open()) {
		getline(file, l);
		while (file.good()) {
			std::string line;
			getline(file, line);
			count_columns(line, cols);
			VectorD row_vector1(cols);
			
			if (row < num) {
				if (not line.empty()) {
					std::istringstream iss(line);
					nodes.push_back(row_vector1);
					
					for (int col = 0; col < cols; ++col) {
						double n;
						iss >> n;
						nodes[row][col] = n;
					}

					row += 1;
				}
			}
			else {
				if (not line.empty()) {
					std::istringstream iss(line);
					double n;
					iss >> n;
					weights.push_back(n);
				}
			}
		}
	}    
	else
		throw "[INPUT ERROR] Unable to open file.";
    
	file.close();

	return SUCCESS;

}


double error_0(VectorD &U,
				int q_order,
				Vector3I &simplices,
				Vector2D &vertices,
				VectorI &num_simplices) {

	size_t N = num_simplices.size();
	double E = 0.0;
	size_t embed_dim = vertices[0].size();

	Vector2D nodes;
	VectorD weights;
	std::string data = "./data/quadrature/d" + std::to_string(N-1) + "o" + std::to_string(q_order) + ".txt";
	read_quadratures(nodes,
					 weights,
					 data);
	size_t nodes_size = nodes.size();

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t i = 0; i < num_simplices[N-1]; ++i) {
		double e = 0;
		double interpolated_U;

		for(size_t j = 0; j < nodes_size; ++j) {
			VectorD vec(embed_dim, 0.0);
			interpolated_U = 0.0;
			for(size_t k = 0; k < N; ++k) {
				interpolated_U += U[simplices[N-1][i][k]] * nodes[j][k];
				for(size_t l = 0; l < embed_dim; ++l) {
					vec[l] += vertices[simplices[N-1][i][k]][l] * nodes[j][k];
				}
			}
			e += weights[j] * pow(interpolated_U - get_analytical_soln(vec), 2);
		}

		
		Vector2D pts;
		for(size_t j = 0; j < N; ++j) {
			pts.push_back(vertices[simplices[N-1][i][j]]);
		}
		double vol = get_simplex_volume(pts);
		#ifdef MULTICORE
			#pragma omp critical
		#endif
		E += sqrt(vol*e);
	}

	return E;
}


int print_vector(Vector2D &vec) {
	size_t vec_size = vec.size();
    for (size_t i = 0; i < vec_size; ++i) {
    	size_t vec_i_size = vec[i].size();
        for (size_t j = 0; j < vec_i_size; ++j) 
            std::cout << vec[i][j] << " ";

        std::cout << std::endl << std::flush;
    }

    return SUCCESS;
}


int print_vector(Vector2I &vec) {
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[i].size(); ++j) 
            std::cout << vec[i][j] << " ";

        std::cout << std::endl << std::flush;
    }

    return SUCCESS;
}


int print_vector(Vector3I &vec) {
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[i].size(); ++j) {
        	for (int k = 0; k < vec[i][j].size(); ++k) {
        		std::cout << vec[i][j][k] << " ";
        	}
        	std::cout << std::endl << std::flush;
        }

        std::cout << std::endl << std::flush;
    }

    return SUCCESS;
}


int print_vector(VectorD &vec) {
    for (int i = 0; i < vec.size(); ++i)
        std::cout << vec[i] << " ";

    std::cout << std::endl << std::flush;

    return SUCCESS;
}

int print_vector(VectorI &vec) {
    for (int i = 0; i < vec.size(); ++i)
        std::cout << vec[i] << " ";

    std::cout << std::endl << std::flush;

    return SUCCESS;
}

int print_vector(VectorMap3I &vec) {
	for (int i = 0; i < vec.size(); ++i)
	{
		for (int j = 0; j < vec[i].size(); ++j)
		{
			for (int k = 0; k < 2; ++k)
			{
				for (auto it = vec[i][j][k].begin(); it != vec[i][j][k].end(); ++it) {
					std::cout<<it->first;	
				}
				std::cout<<"\n";
			}
			std::cout<<"\n";
		}
		std::cout<<"\n";
	}

	return SUCCESS;
}


