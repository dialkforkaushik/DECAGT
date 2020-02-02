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

DenMatD circumcenter_barycentric(DenMatD &pts_matrix, DenMatD &bary_coords) {
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


int get_circumcenter(VectorD &center,
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
						int dim,
						int index) {

	VectorD center;
	double radius;

	Vector2D vertex_pts;

	for (int i = 0; i < dim + 1; ++i) {
		vertex_pts.push_back(vertices[pts[i]]);
	}
	get_circumcenter(center,
				radius,
				vertex_pts);

	temp_centers.push_back(center);

	double vol;
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
							  dim - 1,
							  it->first);

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
		for (int i = 0; i < simplex.size(); ++i) {
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


int print_vector(Vector2D &vec) {
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[i].size(); ++j) 
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


