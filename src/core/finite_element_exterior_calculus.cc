#include "definitions.h"
#include "simplicial_complex.h"
#include "core_utils.h"
#include "discrete_exterior_calculus.h"
#include "geometry.h"
#include "finite_element_exterior_calculus.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <Eigen/Eigen>
#include <cmath>
#include <limits>
#include <iomanip>
#include <numeric>

#ifdef PYTHON
	#include <pybind11/pybind11.h>
	#include <pybind11/stl.h>
	#include <pybind11/numpy.h>
#endif


int set_increasing_ordering(Vector2I &sets) {

	size_t s1 = sets.size();
	size_t s2 = sets[0].size();

	std::vector < std::pair <int, VectorI> > sorted_pair; 
	for(size_t i = 0; i < s1; ++i) {
		int k = 0;
		for(size_t j = 0; j < s2; ++j) {
			if(sets[i][j] != 0) {
				k = 10 * k + sets[i][j];
			}
		}
		sorted_pair.push_back(std::make_pair(k, sets[i]));
	}

	std::sort(sorted_pair.begin(), sorted_pair.end()); 
	
	for(size_t i = 0; i < s1; ++i) {
		sets[i] = sorted_pair[i].second;
	}

	return SUCCESS;
}


int get_sets_sum_to_n(Vector2I &sets,
					  int &n,
					  int &dim,
					  int d = 3) {
	VectorI p(d+1);
    int k = 0;
    p[k] = n;

    while (true) { 
    	int s = std::accumulate(p.begin(), p.end(), 0);
        if (k + 1 == dim || (dim == 0 && s == n)) {
        	sets.push_back(p);
        }
  
        int rem_val = 0; 
        while (k >= 0 && p[k] == 1) { 
            rem_val += p[k]; 
            --k; 
        } 
  
        if (k < 0) {
        	break;
        }
  
        --p[k]; 
        ++rem_val; 
  
        while (rem_val > p[k]) { 
            p[k+1] = p[k]; 
            rem_val = rem_val - p[k]; 
            ++k; 
        } 
  
        p[k+1] = rem_val; 
        ++k; 
    }
}

int FiniteElementExteriorCalculus::omega_ij(double &omega,
											VectorD &bary_coords,
											VectorD &grad_bary_coords) {
	
	if(bary_coords.size() != 2 || grad_bary_coords.size() != 2) {
		return FAILURE;
	}

	omega = bary_coords[0] * grad_bary_coords[1] - bary_coords[1] * grad_bary_coords[0];

	return SUCCESS;
}


int FiniteElementExteriorCalculus::bb_basis(double &bernstein_poly,
											 VectorI &alpha,
											 int n,
											 VectorD &bary_coords) {

	size_t size = alpha.size();

	if(size != bary_coords.size()) {
		return FAILURE;
	}

	double temp_value;

	for(size_t i = 0; i < size; ++i) {
		int fact = alpha[i];
		factorial(fact);
		temp_value = pow(bary_coords[i], alpha[i])/fact;
	}

	factorial(n);
	bernstein_poly = n * temp_value;

	return SUCCESS;
}


int FiniteElementExteriorCalculus::compute_index_sets_p(Vector2I &sets,
												   		 int n,
												   		 int dim,
												   		 int d) {

	if(dim == 3) {
		int num_faces;
		binomialCoeff(num_faces,
					  d+1,
					  dim);

		for(int i = num_faces - 1; i >= 0; --i) {
			Vector2I temp_sets;
			compute_index_sets_o(temp_sets,
						   		 n,
						   		 0,
						   		 dim - 1);
			size_t temp_sets_size = temp_sets.size();
			for(size_t j = 0; j < temp_sets_size - 1; ++j) {
				temp_sets[j].insert(temp_sets[j].begin() + i, 0);
				sets.push_back(temp_sets[j]);
			}
		}
	}

	else if(dim == 4) {
		compute_index_sets_o(sets,
					   		 n,
					   		 0);
		sets.pop_back();
	}

	else {
		return FAILURE;
	}

	return SUCCESS;

}


int FiniteElementExteriorCalculus::compute_index_sets_t(Vector2I &sets,
												   		 int n,
												   		 int dim,
												   		 int d) {

	if(dim != 3) {
		return FAILURE;
	}

	int num_faces;
	binomialCoeff(num_faces,
				  d+1,
				  dim);

	for(int i = num_faces - 1; i >= 0; --i) {
		Vector2I temp_sets;
		Vector2I vec;
		compute_index_sets_o(temp_sets,
					   		 n,
					   		 0,
					   		 dim - 1);
		size_t temp_sets_size = temp_sets.size();
		for(size_t j = 0; j < temp_sets_size; ++j) {
			temp_sets[j].insert(temp_sets[j].begin() + i, 0);
			vec.push_back(temp_sets[j]);
		}
		set_increasing_ordering(vec);
		for(size_t j = 0; j < temp_sets_size; ++j) {
			sets.push_back(vec[j]);
		}
	}

	return SUCCESS;

}

int FiniteElementExteriorCalculus::compute_index_sets_o(Vector2I &sets,
												   		 int n,
												   		 int dim,
												   		 int d) {

  	Vector2I temp_sets;
  	Vector2I face_indices;

  	if(dim == 0) {
  		get_sets_sum_to_n(temp_sets,
  					  n,
  					  dim,
  					  d);

  		get_permutations(sets,
    				 temp_sets);

  		set_increasing_ordering(sets);
  	}

  	else {
  		Vector2I temp;
  		get_sets_sum_to_n(temp,
	  					  dim,
	  					  dim,
	  					  d);

  		get_permutations(face_indices,
    				 	 temp);

  		size_t face_indices_size = face_indices.size();

  		temp.clear();
		get_sets_sum_to_n(temp,
				   		  n,
				   		  dim,
				   		  dim - 1);
		get_permutations(temp_sets,
				 	 	 temp);

		size_t temp_sets_size = temp_sets.size();

  		for(size_t i = 0; i < face_indices_size; ++i) {
  			Vector2I vec2;
			for(size_t k = 0; k < temp_sets_size; ++k) {
				VectorI vec1;
				int index = 0;
				
				for(size_t j = 0; j < d + 1; ++j) {
					if(face_indices[i][j] == 0) {
						vec1.push_back(0);
					}
					else {
						vec1.push_back(temp_sets[k][index]);
						++index;
					}
				}
				
				vec2.push_back(vec1);
			}
			set_increasing_ordering(vec2);
			for(size_t j = 0; j < temp_sets_size; ++j) {
				sets.push_back(vec2[j]);
			}
		}
  	}
  
    // set_ordering(dim, 
    // 			 sets);

    return SUCCESS;
}


int get_triplets(Vector3I &simplex_simplices,
				 int &i,
				 int &k,
				 EigVectorD &vals,
				 VectorTripletD &triplet) {

	int n = simplex_simplices[i][k].size();
	int count = 0;
	for (int v = 0; v < n; ++v) {
		for (int w = 0; w < n; ++w) {
			#ifdef MULTICORE
				#pragma omp critical
			#endif
			triplet.push_back(TripletD(simplex_simplices[i][k][v], simplex_simplices[i][k][w], vals.coeffRef(count)));
			++count;
		}
	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::set_hodge_stars_to_null() {

    for (int i = 0; i < num_simplices[complex_dimension]; ++i) {
    	auto x = std::numeric_limits<SpMatD>::quiet_NaN();
    	hodge_stars.push_back(x);
    }

    return SUCCESS;
}


int FiniteElementExteriorCalculus::compute_hodge_star_k(int &k) {
	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif
	
	Vector3I temp_simplices;
	temp_simplices.push_back(simplex_sorted);
	for (int i = 1; i < complex_dimension + 1; ++i) {
		temp_simplices.insert(temp_simplices.begin(), simplices[i]);
	}


	size_t num_k_simplices = num_simplices[k];
	VectorTripletD triplet;


	if (k < 0 || k > temp_simplices.size()) {
		throw "[INPUT ERROR] Value of k is not in the correct range.";
	}

	int f = k;
	factorial(f);

	double scale_integration = (double)(f*f)/((complex_dimension + 2)*(complex_dimension + 1));

	VectorI nodes;
	nodes.reserve(complex_dimension + 1);
	for (int i = 0; i < complex_dimension + 1; ++i) {
		nodes.push_back(i);
	}

	Vector2I k_forms;
	get_combinations_simplex(nodes,
						     k_forms,
						     k);

	Vector2I k_faces;
	get_combinations_simplex(nodes,
						     k_faces,
						     k+1);

	nodes.clear();
	size_t num_k_faces = k_faces.size();

	Vector3I k_form_pairs;
	get_combinations_simplex(k_forms,
						     k_form_pairs,
						     2);

	size_t num_k_forms = k_forms.size();
	size_t num_k_form_pairs = k_form_pairs.size();

	for (int i = 0; i < num_k_forms; ++i) {
		k_form_pairs.push_back(Vector2I());
		k_form_pairs[num_k_form_pairs + i].push_back(k_forms[i]);
		k_form_pairs[num_k_form_pairs + i].push_back(k_forms[i]);
	}

	num_k_form_pairs = k_form_pairs.size();
	MapVector2I k_form_pairs_to_index;
	
	int count = 0;
	Vector2I temp;
	for (int i = 0; i < num_k_form_pairs; ++i) {
		temp = k_form_pairs[i];
		k_form_pairs_to_index[temp] = count;
		std::reverse(temp.begin(), temp.end());
		k_form_pairs_to_index[temp] = count;
		temp.clear();

		++count;
	}

	size_t num_k_face_pairs = num_k_faces * num_k_faces;

	SpMatD dets_to_vals(num_k_face_pairs, num_k_form_pairs);
	
	Vector3I k_face_pairs;
	Vector2I temp_vec;
	int x = 0;
	for (int i = 0; i < num_k_faces; ++i) {
		for (int j = 0; j < num_k_faces; ++j) {
			
			int row_index = k_face_pairs.size();

			k_face_pairs.push_back(Vector2I());
			k_face_pairs[x].push_back(k_faces[i]);
			k_face_pairs[x].push_back(k_faces[j]);
			++x;

			for (int n = 0; n < k+1; ++n) {
				for (int m = 0; m < k+1; ++m) {
					temp_vec.reserve(2);
					temp_vec.push_back(k_faces[i]);
					temp_vec.push_back(k_faces[j]);

					temp_vec[0].erase(temp_vec[0].begin()+n);
					temp_vec[1].erase(temp_vec[1].begin()+m);

					int col_index = k_form_pairs_to_index[temp_vec];
					temp_vec.clear();

					int equal;
					if (k_faces[i][n] == k_faces[j][m]) {
						equal = 1;
					}
					else {
						equal = 0;
					}

					dets_to_vals.coeffRef(row_index, col_index) += std::pow(-1, n+m) * (equal + 1);
				}
			}
		}
	}

	MapVector2I k_face_pairs_to_index;
	for (int i = 0; i < num_k_faces * num_k_faces; ++i) {
		k_face_pairs_to_index[k_face_pairs[i]] = i;
	}
	
	dets_to_vals.makeCompressed();

	EigVectorD dets;
	size_t num_n_simplex = num_simplices[complex_dimension];
	DenMatD A;
	DenMatD B;
	EigVectorD vals;
	DenMatD d_lambda;
	Vector2D pts;

	#ifdef MULTICORE
		#pragma omp parallel for private(A, B, pts, d_lambda, dets, vals)
	#endif
	for (int i = 0; i < num_n_simplex; ++i) {
		if (k > 0) {
			dets.resize(num_k_form_pairs);
			size_t size = temp_simplices[complex_dimension][i].size();
			for (int j = 0; j < size; ++j) {
				size_t index = temp_simplices[complex_dimension][i][j];
				pts.push_back(vertices[index]);
			}
			barycentric_gradients(pts, d_lambda);

			for (int j = 0; j < num_k_form_pairs; ++j) {

				size_t num_form1 = k_form_pairs[j][0].size();
				size_t num_form2 = k_form_pairs[j][1].size();

				A.resize(num_form1, d_lambda.cols());
				B.resize(num_form2, d_lambda.cols());

				for (int k = 0; k < num_form1; ++k) {
					A.row(k) = d_lambda.row(k_form_pairs[j][0][k]);
				}
				for (int k = 0; k < num_form2; ++k) {
					B.row(k) = d_lambda.row(k_form_pairs[j][1][k]);
				}
				
				dets.coeffRef(j) = (A * B.transpose()).determinant();
				A.resize(0, 0);
				B.resize(0, 0);
			}
			d_lambda.resize(0, 0);
			pts.clear();
		}
		else {
			dets = EigVectorD::Ones(num_k_form_pairs);
		}

		int dim = complex_dimension;

		compute_primal_volume_k(dim,
								i);

		vals = dets_to_vals * dets * primal_volume[dim][i] * scale_integration;

		get_triplets(simplex_sub_simplices, i, k, vals, triplet);
	}
	
	SpMatD mass_matrix(num_k_simplices, num_k_simplices);
	mass_matrix.setFromTriplets(triplet.begin(), triplet.end());

	mass_matrix.makeCompressed();
	hodge_stars[k] = mass_matrix;


	return SUCCESS;
}


int FiniteElementExteriorCalculus::compute_hodge_stars() {
	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif
	
	Vector3I temp_simplices;
	temp_simplices.push_back(simplex_sorted);
	for (int i = 1; i < complex_dimension + 1; ++i) {
		temp_simplices.insert(temp_simplices.begin(), simplices[i]);
	}


	for (int k = 0; k < complex_dimension + 1; ++k) {

		size_t num_k_simplices = num_simplices[k];
		VectorTripletD triplet;


		if (k < 0 || k > temp_simplices.size()) {
			throw "[INPUT ERROR] Value of k is not in the correct range.";
		}

		int f = k;
		factorial(f);

		double scale_integration = (double)(f*f)/((complex_dimension + 2)*(complex_dimension + 1));

		VectorI nodes;
		nodes.reserve(complex_dimension + 1);
		for (int i = 0; i < complex_dimension + 1; ++i) {
			nodes.push_back(i);
		}

		Vector2I k_forms;
		get_combinations_simplex(nodes,
							     k_forms,
							     k);

		Vector2I k_faces;
		get_combinations_simplex(nodes,
							     k_faces,
							     k+1);

		nodes.clear();
		size_t num_k_faces = k_faces.size();

		Vector3I k_form_pairs;
		get_combinations_simplex(k_forms,
							     k_form_pairs,
							     2);

		size_t num_k_forms = k_forms.size();
		size_t num_k_form_pairs = k_form_pairs.size();

		for (int i = 0; i < num_k_forms; ++i) {
			k_form_pairs.push_back(Vector2I());
			k_form_pairs[num_k_form_pairs + i].push_back(k_forms[i]);
			k_form_pairs[num_k_form_pairs + i].push_back(k_forms[i]);
		}

		num_k_form_pairs = k_form_pairs.size();
		MapVector2I k_form_pairs_to_index;
		
		int count = 0;
		Vector2I temp;
		for (int i = 0; i < num_k_form_pairs; ++i) {
			temp = k_form_pairs[i];
			k_form_pairs_to_index[temp] = count;
			std::reverse(temp.begin(), temp.end());
			k_form_pairs_to_index[temp] = count;
			temp.clear();

			++count;
		}

		size_t num_k_face_pairs = num_k_faces * num_k_faces;

		SpMatD dets_to_vals(num_k_face_pairs, num_k_form_pairs);
		
		Vector3I k_face_pairs;
		Vector2I temp_vec;
		int x = 0;
		for (int i = 0; i < num_k_faces; ++i) {
			for (int j = 0; j < num_k_faces; ++j) {
				
				int row_index = k_face_pairs.size();

				k_face_pairs.push_back(Vector2I());
				k_face_pairs[x].push_back(k_faces[i]);
				k_face_pairs[x].push_back(k_faces[j]);
				++x;

				for (int n = 0; n < k+1; ++n) {
					for (int m = 0; m < k+1; ++m) {
						temp_vec.reserve(2);
						temp_vec.push_back(k_faces[i]);
						temp_vec.push_back(k_faces[j]);

						temp_vec[0].erase(temp_vec[0].begin()+n);
						temp_vec[1].erase(temp_vec[1].begin()+m);

						int col_index = k_form_pairs_to_index[temp_vec];
						temp_vec.clear();

						int equal;
						if (k_faces[i][n] == k_faces[j][m]) {
							equal = 1;
						}
						else {
							equal = 0;
						}

						dets_to_vals.coeffRef(row_index, col_index) += std::pow(-1, n+m) * (equal + 1);
					}
				}
			}
		}

		MapVector2I k_face_pairs_to_index;
		for (int i = 0; i < num_k_faces * num_k_faces; ++i) {
			k_face_pairs_to_index[k_face_pairs[i]] = i;
		}
		
		dets_to_vals.makeCompressed();

		EigVectorD dets;
		size_t num_n_simplex = num_simplices[complex_dimension];
		DenMatD A;
		DenMatD B;
		EigVectorD vals;
		DenMatD d_lambda;
		Vector2D pts;

		#ifdef MULTICORE
			#pragma omp parallel for private(A, B, pts, d_lambda, dets, vals)
		#endif
		for (int i = 0; i < num_n_simplex; ++i) {
			if (k > 0) {
				dets.resize(num_k_form_pairs);
				size_t size = temp_simplices[complex_dimension][i].size();
				for (int j = 0; j < size; ++j) {
					size_t index = temp_simplices[complex_dimension][i][j];
					pts.push_back(vertices[index]);
				}
				barycentric_gradients(pts, d_lambda);

				for (int j = 0; j < num_k_form_pairs; ++j) {

					size_t num_form1 = k_form_pairs[j][0].size();
					size_t num_form2 = k_form_pairs[j][1].size();

					A.resize(num_form1, d_lambda.cols());
					B.resize(num_form2, d_lambda.cols());

					for (int k = 0; k < num_form1; ++k) {
						A.row(k) = d_lambda.row(k_form_pairs[j][0][k]);
					}
					for (int k = 0; k < num_form2; ++k) {
						B.row(k) = d_lambda.row(k_form_pairs[j][1][k]);
					}
					
					dets.coeffRef(j) = (A * B.transpose()).determinant();
					A.resize(0, 0);
					B.resize(0, 0);
				}
				d_lambda.resize(0, 0);
				pts.clear();
			}
			else {
				dets = EigVectorD::Ones(num_k_form_pairs);
			}

			int dim = complex_dimension;

			compute_primal_volume_k(dim,
									i);

			vals = dets_to_vals * dets * primal_volume[dim][i] * scale_integration;

			get_triplets(simplex_sub_simplices, i, k, vals, triplet);
		}
		
		SpMatD mass_matrix(num_k_simplices, num_k_simplices);
		mass_matrix.setFromTriplets(triplet.begin(), triplet.end());

		mass_matrix.makeCompressed();
		hodge_stars[k] = mass_matrix;

	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::mass_matrix_bb_0(DenMatD &mass_matrix,
													int n,
													int m,
													int d) {

	int max = std::max(n, m);

	if (max < 1) {
		return FAILURE;
	}

	Vector2I index_sets;
	Vector2I temp_index_sets;
	compute_index_sets_o(index_sets,
						 1,
						 1,
						 d);
	for (int i = 2; i <= std::min(max, d+1); ++i) {
		temp_index_sets.clear();
		compute_index_sets_o(temp_index_sets,
							  max,
							  i,
							  d);
		index_sets.insert(index_sets.end(), temp_index_sets.begin(), temp_index_sets.end());
	}

	size_t size = index_sets.size();
	mass_matrix.resize(size, size);

	for (size_t i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {

			int nCk_1;
			int nCk_2;
			int s1 = std::accumulate(index_sets[i].begin(), index_sets[i].end(), 0);
			int s2 = std::accumulate(index_sets[j].begin(), index_sets[j].end(), 0);
			binomialCoeff(nCk_1,
						  s1 + s2,
						  s2);
			binomialCoeff(nCk_2,
						  s1 + s2 + 3,
						  3);

			double den = nCk_1 * nCk_2;
			double num = 1.0;

			for (int k = 0; k < d + 1; ++k) {
				int nCk;
				binomialCoeff(nCk,
							  index_sets[i][k] + index_sets[j][k],
							  index_sets[i][k]);
				num *= nCk;
			}

			mass_matrix.coeffRef(i, j) = num/den;
		}
	}

	return SUCCESS;
	
}

FiniteElementExteriorCalculus::FiniteElementExteriorCalculus() {

    set_hodge_stars_to_null();
    
}


FiniteElementExteriorCalculus::FiniteElementExteriorCalculus(SimplicialComplex sc) : GeometryComplex(sc) {

    set_hodge_stars_to_null();
}

FiniteElementExteriorCalculus::~FiniteElementExteriorCalculus() {}