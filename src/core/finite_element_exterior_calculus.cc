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
					  int n,
					  int dim,
					  int d = 3) {

	VectorI p(std::max(n, d+1));
    int k = 0;
    p[k] = n;

    while (true) { 
    	int s = std::accumulate(p.begin(), p.begin() + d + 1, 0);

        if ((k + 1 == dim && s == n) || (dim == 0 && s == n)) {
        	VectorI sliced = VectorI(p.begin(), p.begin() + d + 1);
        	sets.push_back(sliced);
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

    return SUCCESS;
}


int FiniteElementExteriorCalculus::phi_FT(double &phi,
										  VectorI &alpha,
										  int n,
										  VectorD &bary_coords,
										  VectorD &grad_bary_coords,
										  VectorI &local_indices) {

	size_t size = alpha.size();

	if(size != bary_coords.size()) {
		return FAILURE;
	}

	// double temp_value = 1;
	// double bernstein_poly;

	// for(size_t i = 0; i < size; ++i) {
	// 	int fact = alpha[i];
	// 	factorial(fact);
	// 	temp_value *= pow(bary_coords[i], alpha[i])/fact;
	// }

	// int fact = n;
	// factorial(fact);
	// bernstein_poly = fact * temp_value;

	// temp_value = 0;
	// size_t local_indices_size = local_indices.size();

	// for(size_t i = 0; i < local_indices_size; ++i) {
	// 	double omega = 0;
	// 	omega = bary_coords[(i+1)%local_indices_size] * grad_bary_coords[(i+2)%local_indices_size] - bary_coords[(i+2)%local_indices_size] * grad_bary_coords[(i+1)%local_indices_size];
	// 	temp_value += alpha[i]*omega;
	// }

	// phi = (n+1) * bernstein_poly * temp_value;

	Vector2I SIGMA;
	Vector2I temp_local_indices;

	temp_local_indices.push_back(local_indices);
	get_permutations(SIGMA,
					 temp_local_indices);

	size_t SIGMA_size = SIGMA.size();

	phi = 0;

	for (size_t i = 0; i < SIGMA_size; ++i) {
		VectorI sigma = SIGMA[i];

		double c = (alpha[sigma[0]] + 1) * (alpha[sigma[2]] * grad_bary_coords[1] - alpha[sigma[1]] * grad_bary_coords[2]);

		double temp_phi = 0;
		VectorI temp_alpha = alpha;
		temp_alpha[sigma[0]] += 1;

		bernstein(temp_phi,
				 temp_alpha,
				 n + 1,
				 bary_coords);

		phi += c * temp_phi;

	}

	return SUCCESS;
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


int FiniteElementExteriorCalculus::bernstein(double &bernstein_poly,
											 VectorI &alpha,
											 int n,
											 VectorD &bary_coords) {

	size_t size = alpha.size();

	if(size != bary_coords.size()) {
		return FAILURE;
	}

	double temp_value = 1;

	for(size_t i = 0; i < size; ++i) {
		long long fact = alpha[i];
		factorial(fact);
		temp_value *= pow(bary_coords[i], alpha[i])/fact;
	}

	int sum = std::accumulate(alpha.begin(), alpha.end(), 0);
	
	long long num = n;
	factorial(num);
	long long den = n - sum;
	factorial(den);

	bernstein_poly = (num/den) * temp_value;

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
  	sets.clear();

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
  		if(n < dim) {
			return FAILURE;
		}

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
			barycentric_gradients(d_lambda, pts);

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
				barycentric_gradients(d_lambda, pts);

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


int FiniteElementExteriorCalculus::M_alpha_beta(double &M,
										  		VectorI &alpha,
										  		VectorI &beta) {

	int n = std::accumulate(alpha.begin(), alpha.end(), 0);
	int m = std::accumulate(beta.begin(), beta.end(), 0);

	size_t d = alpha.size();

	long long nCk_1;
	long long nCk_2;

	binomialCoeff(nCk_1,
				  n + m,
				  n);
	binomialCoeff(nCk_2,
				  n + m + 3,
				  3);

	double den = nCk_1 * nCk_2;
	double num = 1.0;

	for (size_t i = 0; i < d; ++i) {
		long long nCk;
		binomialCoeff(nCk,
					  alpha[i] + beta[i],
					  alpha[i]);
		num *= nCk;
	}

	M = num/den;

	return SUCCESS;
}


int FiniteElementExteriorCalculus::S_n(DenMatD &S,
									   Vector2D &pts,
									   int n) {

	if (n == 0) {
		DenMatD gradients;
		barycentric_gradients(gradients,
							  pts);

		size_t gradients_size = gradients.rows();
		S.resize(gradients_size, gradients_size);

		for (size_t i = 0; i < gradients_size; ++i) {
			for (size_t j = 0; j <=i; ++j) {
				S.coeffRef(i,j) = gradients.row(i).dot(gradients.row(j));
				if (i != j) {
					S.coeffRef(j, i) = S.coeffRef(i, j);
				}
			}
		}
	}
		
	return SUCCESS;

}


int FiniteElementExteriorCalculus::mass_matrix_bb_curl(DenMatD &mass_matrix,
													   Vector2D &pts,
													   int n,
													   int m) {

	int max = std::max(n, m);

	if (max < 0) {
		return FAILURE;
	}

	Vector2I alpha;
	Vector2I temp_alpha;
	VectorI ordered_basis_sizes;
	Vector2I e;

	compute_index_sets_o(e,
						 1,
						 1);

	compute_index_sets_o(alpha,
						 2,
						 2);
	ordered_basis_sizes.push_back(alpha.size());

	size_t total = alpha.size();
	for (int i = 2; i <= 4; ++i) {
		temp_alpha.clear();
		compute_index_sets_o(temp_alpha,
							  max + 1,
							  i);
		
		size_t temp_alpha_size = temp_alpha.size();
		if (temp_alpha_size == 0) {
			continue;
		}
		ordered_basis_sizes.push_back(total + temp_alpha_size);
		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
		total += temp_alpha_size;
	}

	size_t size = alpha.size();
	mass_matrix.resize(size, size);

	DenMatD S_0;
	S_n(S_0,
		pts,
		0);

	for (size_t i = 0; i < size; ++i) {
		for (size_t j = i; j < size; ++j) {

			double x = 0;
			
			// 33a
			if (i < ordered_basis_sizes[0] &&
				j < ordered_basis_sizes[0]) {

				VectorI index_p;
				VectorI index_q;

				for (int k = 0; k < 4; ++k) {
					if (alpha[i][k] > 0) {
						index_p.push_back(k);
					}
					if (alpha[j][k] > 0) {
						index_q.push_back(k);
					}
				}

				double M_ii;
				M_alpha_beta(M_ii,
							 e[index_p[0]],
							 e[index_q[0]]);
				double M_ij;
				M_alpha_beta(M_ij,
							 e[index_p[0]],
							 e[index_q[1]]);
				
				double M_ji = M_ij;

				double M_jj;
				M_alpha_beta(M_jj,
							 e[index_p[1]],
							 e[index_q[1]]);
				
				double S_ii = S_0(index_p[0], index_q[0]);
				double S_ij = S_0(index_p[0], index_q[1]);
				double S_ji = S_0(index_p[1], index_q[0]);
				double S_jj = S_0(index_p[1], index_q[1]);
				

				x = (M_ii * S_jj) - (M_ij * S_ji) - (M_ji * S_ij) + (M_jj * S_ii);
			}

			// 33b
			else if(i < ordered_basis_sizes[0] && 
				    j >= ordered_basis_sizes[0] && 
				    j < ordered_basis_sizes[1]) {
				
				x = 0;
				for (int s = 0; s < 4; ++s) {
					VectorI vec;
					std::transform(alpha[j].begin(), alpha[j].end(), e[s].begin(),
									std::back_inserter(vec),
									std::minus<int>());
					
					bool flag = false;
					for (int k = 0; k < 4; ++k) {
						if (vec[k] < 0) {
							flag = true;
							break;
						}
					}

					if (flag) {
						continue;
					}

					VectorI index_p;

					for (int k = 0; k < 4; ++k) {
						if (alpha[i][k] > 0) {
							index_p.push_back(k);
						}
					}

					double M_i;
					M_alpha_beta(M_i,
								 e[index_p[0]],
								 vec);

					double M_j;
					M_alpha_beta(M_j,
								 e[index_p[1]],
								 vec);

					double S_i = S_0.coeffRef(index_p[0], s);
					double S_j = S_0.coeffRef(index_p[1], s);

					x += (M_i * S_j) - (M_j * S_i);
				}

				x = x * (max + 1);
			}

			// 33c
			else if(i < ordered_basis_sizes[0] && 
				    j >= ordered_basis_sizes[1] && 
				    j < ordered_basis_sizes[2]) {
				
				x = 0;
				
				Vector2I SIGMA;
				Vector2I temp_local_indices;
				VectorI local_indices;

				for (int k = 0; k < 4; ++k) {
					if (alpha[j][k] > 0) {
						local_indices.push_back(k);
					}
				}

				temp_local_indices.push_back(local_indices);
				get_permutations(SIGMA,
								 temp_local_indices);

				size_t SIGMA_size = SIGMA.size();

				VectorI index_p;

				for (int k = 0; k < 4; ++k) {
					if (alpha[i][k] > 0) {
						index_p.push_back(k);
					}
				}

				for (int k = 0; k < SIGMA_size; ++k) {
					VectorI vec;
					std::transform(alpha[j].begin(), alpha[j].end(), e[SIGMA[k][0]].begin(),
										std::back_inserter(vec),
										std::plus<int>());

					double M_i;
					M_alpha_beta(M_i,
								 e[index_p[0]],
								 vec);

					double M_j;
					M_alpha_beta(M_j,
								 e[index_p[1]],
								 vec);

					double sigma_1 = SIGMA[k][0];
					double sigma_2 = SIGMA[k][1];
					double sigma_3 = SIGMA[k][2];

					x += (M_i * (alpha[j][sigma_1] + 1) * (alpha[j][sigma_3] * S_0.coeffRef(index_p[1], sigma_2) - alpha[j][sigma_2] * S_0.coeffRef(index_p[1], sigma_3)))
						 - (M_j * (alpha[j][sigma_1] + 1) * (alpha[j][sigma_3] * S_0.coeffRef(index_p[0], sigma_2) - alpha[j][sigma_2] * S_0.coeffRef(index_p[0], sigma_3)));
				}
			}

			// 33e
			else if(i >= ordered_basis_sizes[0] &&
					i < ordered_basis_sizes[1] && 
				    j >= ordered_basis_sizes[0] && 
				    j < ordered_basis_sizes[1]) {
				
				x = 0;
				for (int t = 0; t < 4; ++t) {
					VectorI vec1;
					std::transform(alpha[i].begin(), alpha[i].end(), e[t].begin(),
										std::back_inserter(vec1),
										std::minus<int>());

					for (int s = 0; s < 4; ++s) {
						VectorI vec2;
						std::transform(alpha[j].begin(), alpha[j].end(), e[s].begin(),
										std::back_inserter(vec2),
										std::minus<int>());
						
						bool flag = false;
						for (int k = 0; k < 4; ++k) {
							if (vec1[k] < 0 || vec2[k] < 0) {
								flag = true;
								break;
							}
						}

						if (flag) {
							continue;
						}

						double M;
						M_alpha_beta(M,
									 vec1,
									 vec2);

						double S = S_0.coeffRef(t, s);

						x += (M * S);
					}
				}

				x = x * pow(max + 1, 2);
			}

			// 33f
			else if(i >= ordered_basis_sizes[0] && 
					i < ordered_basis_sizes[1] &&
				    j >= ordered_basis_sizes[1] && 
				    j < ordered_basis_sizes[2]) {
				
				x = 0;
				
				Vector2I SIGMA;
				Vector2I temp_local_indices;
				VectorI local_indices;

				for (int k = 0; k < 4; ++k) {
					if (alpha[j][k] > 0) {
						local_indices.push_back(k);
					}
				}

				temp_local_indices.push_back(local_indices);
				get_permutations(SIGMA,
								 temp_local_indices);

				size_t SIGMA_size = SIGMA.size();

				for (int t = 0; t < 4; ++t) {
					VectorI vec1;
					std::transform(alpha[i].begin(), alpha[i].end(), e[t].begin(),
										std::back_inserter(vec1),
										std::minus<int>());

					for (int k = 0; k < SIGMA_size; ++k) {
						VectorI vec2;
						std::transform(alpha[j].begin(), alpha[j].end(), e[SIGMA[k][0]].begin(),
											std::back_inserter(vec2),
											std::plus<int>());

						double M;
						M_alpha_beta(M,
									 vec1,
									 vec2);

						double sigma_1 = SIGMA[k][0];
						double sigma_2 = SIGMA[k][1];
						double sigma_3 = SIGMA[k][2];

						x += (M * (alpha[j][sigma_1] + 1) * (alpha[j][sigma_3] * S_0.coeffRef(t, sigma_2) - alpha[j][sigma_2] * S_0.coeffRef(t, sigma_3)));
					}
				}

				x = (max + 1) * x;
			}

			// 33h
			else if(i >= ordered_basis_sizes[1] && 
					i < ordered_basis_sizes[2] &&
				    j >= ordered_basis_sizes[1] && 
				    j < ordered_basis_sizes[2]) {
				
				x = 0;
				
				Vector2I SIGMA_p;
				Vector2I SIGMA_q;
				Vector2I temp_local_indices;
				VectorI local_indices;

				for (int k = 0; k < 4; ++k) {
					if (alpha[i][k] > 0) {
						local_indices.push_back(k);
					}
				}

				temp_local_indices.push_back(local_indices);
				get_permutations(SIGMA_p,
								 temp_local_indices);

				local_indices.clear();
				temp_local_indices.clear();

				for (int k = 0; k < 4; ++k) {
					if (alpha[j][k] > 0) {
						local_indices.push_back(k);
					}
				}

				temp_local_indices.push_back(local_indices);
				get_permutations(SIGMA_q,
								 temp_local_indices);

				size_t SIGMA_p_size = SIGMA_p.size();
				size_t SIGMA_q_size = SIGMA_q.size();

				for (int t = 0; t < SIGMA_p_size; ++t) {
					VectorI vec1;
					std::transform(alpha[i].begin(), alpha[i].end(), e[SIGMA_p[t][0]].begin(),
										std::back_inserter(vec1),
										std::plus<int>());

					for (int s = 0; s < SIGMA_q_size; ++s) {
						VectorI vec2;
						std::transform(alpha[j].begin(), alpha[j].end(), e[SIGMA_q[s][0]].begin(),
											std::back_inserter(vec2),
											std::plus<int>());

						double M;
						M_alpha_beta(M,
									 vec1,
									 vec2);

						x += M * (alpha[i][SIGMA_p[t][0]] + 1) * (alpha[j][SIGMA_q[s][0]] + 1)  
						     * (alpha[i][SIGMA_p[t][2]] * alpha[j][SIGMA_q[s][2]] * S_0.coeffRef(SIGMA_p[t][1], SIGMA_q[s][1])
						     	- alpha[i][SIGMA_p[t][2]] * alpha[j][SIGMA_q[s][1]] * S_0.coeffRef(SIGMA_p[t][1], SIGMA_q[s][2])
						     	- alpha[i][SIGMA_p[t][1]] * alpha[j][SIGMA_q[s][2]] * S_0.coeffRef(SIGMA_p[t][2], SIGMA_q[s][1])
						     	+ alpha[i][SIGMA_p[t][1]] * alpha[j][SIGMA_q[s][1]] * S_0.coeffRef(SIGMA_p[t][2], SIGMA_q[s][2]));
					}
				}
			}

			mass_matrix.coeffRef(i, j) = x;
			if (i != j) {
				mass_matrix.coeffRef(j, i) = x;
			}
		}
	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::mass_matrix_bb_1(DenMatD &mass_matrix,
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
						 1);

	for (int i = 2; i <= std::min(max, d+1); ++i) {
		temp_index_sets.clear();
		compute_index_sets_o(temp_index_sets,
							  max,
							  i);
		index_sets.insert(index_sets.end(), temp_index_sets.begin(), temp_index_sets.end());
	}

	size_t size = index_sets.size();
	mass_matrix.resize(size, size);

	for (size_t i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {

			long long nCk_1;
			long long nCk_2;
			
			int s1 = std::accumulate(index_sets[i].begin(), index_sets[i].end(), 0);
			int s2 = std::accumulate(index_sets[j].begin(), index_sets[j].end(), 0);

			binomialCoeff(nCk_1,
						  s1 + s2,
						  s1);
			binomialCoeff(nCk_2,
						  s1 + s2 + 3,
						  3);

			double den = nCk_1 * nCk_2;
			double num = 1.0;

			for (int k = 0; k < d + 1; ++k) {
				long long nCk;
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

double FiniteElementExteriorCalculus::bb_error(int n,
											    Vector3I &simplices,
												Vector2D &vertices,
												VectorI &num_simplices,
												int q_order) {

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

	Vector2I alpha;
	Vector2I temp_alpha;
	Vector2D basis_vector;

	compute_index_sets_o(alpha,
					 	 1,
					 	 1);

	int d = alpha[0].size();

	for(int i = 2; i <= std::min(n, d); ++i) {
		temp_alpha.clear();
		compute_index_sets_o(temp_alpha,
						 	 n,
						 	 i);

		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	}
	size_t alpha_size = alpha.size();

	DenMatD basis_elements(alpha_size, nodes_size);
	for(size_t i = 0; i < alpha_size; ++i) {
		EigVectorD temp_basis_elements(nodes_size);
		int sum;
		get_sum(sum,
				alpha[i]);
		
		for(size_t j = 0; j < nodes_size; ++j) {
			double bernstein_poly = 0;
			bernstein(bernstein_poly,
					 alpha[i],
					 sum,
					 nodes[j]);
			temp_basis_elements(j) = bernstein_poly;
		}

		basis_elements.row(i) = temp_basis_elements;
	}

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t i = 0; i < num_simplices[N-1]; ++i) {
		double e = 0;

		Vector2D pts;
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplices[N-1][i][k]]);
		}
		double vol = get_simplex_volume(pts);

		DenMatD M;
		mass_matrix_bb_1(M,
						 n,
						 n);
		M = M * vol;

		double sum_weight = 0.0;

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {
			EigVectorD b(alpha_size);

			for(size_t j = 0; j < alpha_size; ++j) {
				double inner_product = 0.0;
				double sum_weights = 0.0;

				for(size_t k = 0; k < nodes_size; ++k) {
					VectorD vec(embed_dim, 0.0);

					for(size_t v = 0; v < N; ++v) {
						for(size_t l = 0; l < embed_dim; ++l) {
							vec[l] += pts[v][l] * nodes[k][v];
						}
					}

					sum_weights += weights[k];
					inner_product += vol * weights[k] * get_analytical_soln(vec) * basis_elements.coeffRef(j, k);
				}

				b.coeffRef(j) = inner_product/sum_weights;
			}

			EigVectorD coeffs = M.llt().solve(b);

			double f_dash = coeffs.dot(basis_elements.col(node_index));
			
			VectorD points(embed_dim, 0.0);
			for(size_t v = 0; v < N; ++v) {
				for(size_t l = 0; l < embed_dim; ++l) {
					points[l] += pts[v][l] * nodes[node_index][v];
				}
			}

			e += weights[node_index] * pow(get_analytical_soln(points) - f_dash, 2);
			sum_weight += weights[node_index];
		}

		#ifdef MULTICORE
			#pragma omp critical
		#endif
		E += sqrt(vol*e/sum_weight);
	}

	return E;
}


FiniteElementExteriorCalculus::FiniteElementExteriorCalculus() {

    set_hodge_stars_to_null();
    
}


FiniteElementExteriorCalculus::FiniteElementExteriorCalculus(SimplicialComplex sc) : GeometryComplex(sc) {

    set_hodge_stars_to_null();
}

FiniteElementExteriorCalculus::~FiniteElementExteriorCalculus() {}