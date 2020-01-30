#include "definitions.h"
#include "simplicial_complex.h"
#include "core_utils.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <Eigen/Eigen>
#include <set>


int SimplicialComplex::build_complex() {
	
	compute_simplices();
    compute_elements();
    compute_simplex_simplices();
	compute_boundary_matrices();
	compute_adjacency1d();

	return SUCCESS;
}


int SimplicialComplex::compute_simplices() {
	int n = simplex[0].size();
	int N = simplex.size();
	simplices.reserve(n);

    simplex_sorted = simplex;
	for (int i = 0; i < N; ++i) {
		std::sort(simplex_sorted[i].begin(), simplex_sorted[i].end());
	}

    Vector2I complex;
    simplices.push_back(simplex);

	num_simplices.push_back(simplex.size());

	for (int i = 1; i < n; ++i) {
		get_combinations_mesh(simplices[0], complex);
		simplices.insert(simplices.begin(), complex);
		num_simplices.insert(num_simplices.begin(), complex.size());
		complex.clear();
	}

    complex_dimension = simplices.size() - 1;
    embedding_dimension = vertices[0].size();
        
	return SUCCESS;
}



int SimplicialComplex::compute_boundary_matrices() {
	int n = complex_dimension + 1;
	boundary_matrices.reserve(n - 1);

	SpMatI boundary_matrix;
	SpMatIC boundary_matrix_col_major;
	VectorTripletI triplet;

	for (int i = 1; i < n; ++i) {
		Vector2I complex_highdim = simplices[i];
		Vector2I complex_lowdim = simplices[i - 1];

		int num_complex_highdim = complex_highdim.size();
		int num_complex_lowdim = complex_lowdim.size();

		boundary_matrix.resize(num_complex_lowdim,
                               num_complex_highdim);
		boundary_matrix_col_major.resize(num_complex_lowdim,
                               num_complex_highdim);
		#ifdef MULTICORE
			#if MULTICORE
				#pragma omp parallel for shared(triplet)
			#endif
		#endif

		for (int j = 0; j < num_complex_highdim; ++j) {
			Vector2I combinations;
            get_combinations_simplex(complex_highdim[j], combinations);
            int N = combinations.size();

            int x;

			for (int k = 0; k < N; ++k) {
				int index = elements[combinations[k]];

				if (complex_highdim[0].size() > 2) {
					if (is_rotated(complex_highdim[j],
                                   combinations[k]) == 1)
						x = -1;
					else
						x = 1;
				}
				else {
					if(combinations[k][0] == complex_highdim[j][0] && combinations[k][1] == complex_highdim[j][1])
						x = 1;
					else
						x = -1;
				}
				#ifdef MULTICORE
					#if MULTICORE 
						#pragma omp critical
					#endif
				#endif
				triplet.push_back(TripletI(index, j, x));
			}
		}
		boundary_matrix.setFromTriplets(triplet.begin(), triplet.end());
		boundary_matrix_col_major.setFromTriplets(triplet.begin(), triplet.end());

		boundary_matrix.makeCompressed();
		boundary_matrices.push_back(boundary_matrix);
		boundary_matrix_col_major.makeCompressed();
		boundary_matrices_col_major.push_back(boundary_matrix_col_major);

		boundary_matrix.resize(0, 0);
		boundary_matrix_col_major.resize(0, 0);
		triplet.clear();
	}

	return SUCCESS;
}


int SimplicialComplex::compute_adjacency1d() {

	MapI simplex_map;
	for (int i = 0; i < complex_dimension + 1; ++i) {
		adjacency1d.push_back(VectorMap2I());

		for (int j = 0; j < num_simplices[i]; ++j) {
			adjacency1d[i].push_back(VectorMapI());
			if (i == 0) {
				simplex_map.clear();
				adjacency1d[i][j].push_back(simplex_map);

				for (SpMatI::InnerIterator it(boundary_matrices[i],j); it; ++it) {
					if (it.row() == j) {
						simplex_map[it.col()] = j;
					}
				}
				
				adjacency1d[i][j].push_back(simplex_map);
				simplex_map.clear();
			}
			else if (i == complex_dimension) {
				simplex_map.clear();

				int N = simplex_simplices[j][complex_dimension - 1].size();
				for (int k = 0; k < N; ++k) {
					simplex_map[simplex_simplices[j][complex_dimension - 1][k]] = j;
				}

				adjacency1d[i][j].push_back(simplex_map);
				simplex_map.clear();
				adjacency1d[i][j].push_back(simplex_map);
			}
			else {
				simplex_map.clear();
				for (SpMatIC::InnerIterator it(boundary_matrices_col_major[i-1],j); it; ++it) {
					if (it.col() == j) {
						simplex_map[it.row()] = j;
					}
				}

				adjacency1d[i][j].push_back(simplex_map);
				simplex_map.clear();

				for (SpMatI::InnerIterator it(boundary_matrices[i],j); it; ++it) {
					if (it.row() == j) {
						simplex_map[it.col()] = j;
					}
				}
				adjacency1d[i][j].push_back(simplex_map);
			}
		}
	}

	return SUCCESS;
}


int SimplicialComplex::compute_adjacency2d() {

	for (int i = 0; i < complex_dimension + 1; ++i) {
		adjacency2d.push_back(VectorMap2I());

		for (int j = 0; j < simplices[i].size(); ++j) {
			adjacency2d[i].push_back(VectorMapI());
			MapI simplex_map;
			if (i == 0) {
				auto vec = boundary_matrices[i].row(j);
				adjacency2d[i][j].push_back(simplex_map);
				for (int k = 0; k < vec.size(); ++k)
				{
			        if (vec.coeff(k) != 0) {
			        	for (auto it = adjacency1d[i+1][k][0].begin(); it != adjacency1d[i+1][k][0].end(); ++it) {
			        		if (it->first != j) {
								simplex_map[it->first] = j;
			        		}
						}
					}
				}
				adjacency2d[i][j].push_back(simplex_map);
			}
			else if (i == complex_dimension) {
				auto vec = boundary_matrices[i - 1].col(j);
				for (int k = 0; k < vec.size(); ++k)
				{
			        if (vec.coeff(k) != 0) {
			        	for (auto it = adjacency1d[i-1][k][1].begin(); it != adjacency1d[i-1][k][1].end(); ++it) {
							if (it->first != j) {
								simplex_map[it->first] = j;
			        		}
						}
					}
				}
				adjacency2d[i][j].push_back(simplex_map);
			}
			else {
				auto vec = boundary_matrices[i-1].col(j);
				for (int k = 0; k < vec.size(); ++k)
				{
			        if (vec.coeff(k) != 0) {
			        	for (auto it = adjacency1d[i-1][k][1].begin(); it != adjacency1d[i-1][k][1].end(); ++it) {
							if (it->first != j) {
								simplex_map[it->first] = j;
			        		}
						}
					}
				}
				adjacency2d[i][j].push_back(simplex_map);
				simplex_map.clear();

				auto vec1 = boundary_matrices[i].row(j);
				for (int k = 0; k < vec1.size(); ++k)
				{
			        if (vec1.coeff(k) != 0) {
			        	for (auto it = adjacency1d[i+1][k][0].begin(); it != adjacency1d[i+1][k][0].end(); ++it) {
							if (it->first != j) {
								simplex_map[it->first] = j;
			        		}
						}
					}
				}
				adjacency2d[i][j].push_back(simplex_map);
			}
		}
	}

	return SUCCESS;
}


inline int SimplicialComplex::compute_elements() {
	for (int i = 0; i < complex_dimension + 1; ++i) {
		int N = simplices[i].size();
		for (int j = 0; j < N; ++j) {
			elements[simplices[i][j]] = j;
		}
	}

	return SUCCESS;
}


inline int SimplicialComplex::compute_simplex_simplices() {
	int n = simplex.size();
	int N = simplex[0].size();

	simplex_simplices.reserve(n);

    Vector2I complex;

	for (int i = 0; i < n; ++i) {
		simplex_simplices.push_back(Vector2I());
		VectorI vec;
		vec.reserve(1);
		vec.push_back(i);
		simplex_simplices[i].push_back(vec);
		vec.clear();

		for (int j = N-1; j > 0; --j) {
			get_combinations_simplex(simplices[complex_dimension][i], complex, j);
			VectorI vec;
			for (int k = 0; k < complex.size(); ++k) {
				vec.push_back(elements[complex[k]]);
			}
			simplex_simplices[i].insert(simplex_simplices[i].begin(), vec);
			vec.clear();
			
			complex.clear();
		}
	}

	return SUCCESS;
}


Vector2D SimplicialComplex::circumcenter(int dim) {
	
	VectorD center;
	double radius = -1;
	Vector2D pts;
	Vector2D circumcenters; 

	for (int i = 0; i < num_simplices[dim]; ++i) {
		pts.clear();
		for (int j = 0; j < simplices[dim][i].size(); ++j) {
			pts.push_back(vertices[simplices[dim][i][j]]);
			
		}
		get_circumcenter(center,
					     radius,
					     pts);
		circumcenters.push_back(center);
		center.clear();
		radius = -1;
	}

	return circumcenters;
}


int SimplicialComplex::read_files(int &rows, VectorS &files) {
	if(rows == 1) {
		std::string off_file = files[0];

		std::ifstream f1(off_file);
		std::string line;

		getline(f1, line);
		getline(f1, line);
		getline(f1, line);

		int v_cols;
		count_columns(line, v_cols);

	    f1.close();

		int num_vertices;
		int num_triangles;
		int num_edges;

		std::fstream file;
		std::fstream f2;
		std::string f;
		int s_cols;
		std::vector<double> row_vector1(v_cols);
		std::vector<int> row_vector2;
		int row = 0;

		file.open(off_file, std::ios::in);
		f2.open(off_file, std::ios::in);
		if (file.is_open()) {
			std::string temp;
			getline(file, temp);
			file >> num_vertices;
			file >> num_triangles;
			file >> num_edges;

			getline(f2, f);
			getline(f2, f);

			while (file.good()) {
				if (row == num_vertices) {
					getline(f2, f);
					count_columns(f, s_cols);
					row_vector2.resize(s_cols);
				}

				if (row < num_vertices) {
					getline(file, line);
					std::istringstream iss(line);
					vertices.push_back(row_vector1);
					getline(f2, f);
					for (int col = 0; col < v_cols; ++col) 
						iss >> vertices[row][col];
				}
				else if (row >= num_vertices &&
	                     row < num_vertices + num_triangles) {
					getline(file, line);
					std::istringstream iss(line);
					simplex.push_back(row_vector2);
					int r = row - num_vertices;
					for (int col = 0; col < s_cols; ++col) {
						iss >> simplex[r][col];
					}
				}
				else
					break;
	            
				row += 1;
			}
		}
		else
			throw "[INPUT ERROR] Unable to open file.";

		file.close();
		f2.close();
	}

	else if(rows == 2) {
		std::string v_file = files[0];
		std::string s_file = files[1];

		std::ifstream f1(v_file);
		std::string line_v;
		getline(f1, line_v);
		int v_cols;
		count_columns(line_v, v_cols);
		f1.close();

		std::fstream file;
		VectorD row_vector1(v_cols);
		int row = 0;

		file.open(v_file, std::ios::in);
		if (file.is_open()) {
			while (file.good()) {
				std::string line;
				getline(file, line);
				
				if (not line.empty()) {
					std::istringstream iss(line);
					vertices.push_back(row_vector1);
					
					for (int col = 0; col < v_cols; ++col) {
						double n;
						iss >> n;
						vertices[row][col] = n;
					}

					row += 1;
				}
			}
		}    
		else
			throw "[INPUT ERROR] Unable to open file.";
	    
		file.close();

	    std::ifstream f2(s_file);
		std::string line_t;
		getline(f2, line_t);
		int t_cols;
		count_columns(line_t, t_cols);
		f2.close();

	    VectorI row_vector2(t_cols);
		row = 0;

		file.open(s_file, std::ios::in);
		if (file.is_open()) {
			while (file.good()) {
				std::string line;
				getline(file, line);
				
				if(not line.empty()) {
					std::istringstream iss(line);
					simplex.push_back(row_vector2);
					
					for (int col = 0; col < t_cols; ++col) {
						int n;
						iss >> n;
						simplex[row][col] = n;
					}

					row += 1;
				}
			}
		}
		else
			std::cout << "[INPUT ERROR] Unable to open file." << std::endl
	                  << std::flush;

		file.close();
	}

	else {
		throw std::invalid_argument("\n[INPUT ERROR] Input file format not correct.\n");
	}

	return SUCCESS;
}


SimplicialComplex::SimplicialComplex() {
	const std::string input_file = "input.txt";

	std::fstream inp_file;
	VectorS files;
	int rows = 0;

	inp_file.open(input_file, std::ios::in);
	if (inp_file.is_open()) {
		while (inp_file.good()) {
			std::string line;
			getline(inp_file, line);
			if (not line.empty()) {
				files.push_back(line);
				rows += 1;
			}
		}
	}

	read_files(rows, files);
    
	inp_file.close();
}


SimplicialComplex::SimplicialComplex(Vector2D &v, 
									 Vector2I &f) {
	
	vertices = v;
	simplex = f;
}


SimplicialComplex::~SimplicialComplex() {}







