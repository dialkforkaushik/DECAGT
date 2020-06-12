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


int FiniteElementExteriorCalculus::psi_T(EigVectorD &psi,
										  VectorI &alpha,
										  int n,
										  int l,
										  VectorD &bary_coords,
										  DenMatD &grad_bary_coords) {

	psi = EigVectorD::Zero(grad_bary_coords.cols());

	for (int i = 0; i < 4; ++i) {
		EigVectorD c;
		int delta = 1;
		double B = 0;

		VectorI temp_alpha = alpha;
		temp_alpha[i] -= 1;

		if (temp_alpha[i] < 0) {
			continue;
		}

		if (i != l) {
			delta = 0;
		}

		c = (delta * (n + 1) - alpha[i]) * grad_bary_coords.row(i);

		bernstein(B,
				  temp_alpha,
				  n,
				  bary_coords);

		psi += c * B;
	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::phi_FT(EigVectorD &phi,
										  VectorI &alpha,
										  int n,
										  VectorD &bary_coords,
										  DenMatD &grad_bary_coords,
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

	for (int k = 0; k < 3; ++k) {
		VectorI vec {local_indices[k%3], local_indices[(k+1)%3], local_indices[(k+2)%3]};
		SIGMA.push_back(vec);
	}

	size_t SIGMA_size = SIGMA.size();

	phi = EigVectorD::Zero(grad_bary_coords.cols());

	for (size_t i = 0; i < SIGMA_size; ++i) {
		VectorI sigma = SIGMA[i];

		EigVectorD c = (alpha[sigma[0]] + 1) * (alpha[sigma[2]] * grad_bary_coords.row(sigma[1]) - alpha[sigma[1]] * grad_bary_coords.row(sigma[2]));

		double B = 0;
		VectorI temp_alpha = alpha;
		temp_alpha[sigma[0]] += 1;

		bernstein(B,
				  temp_alpha,
				  n + 1,
				  bary_coords);

		phi += c * B;

	}

	return SUCCESS;
}

int FiniteElementExteriorCalculus::grad_B(EigVectorD &grad_b,
										  VectorI &alpha,
										  int n,
										  VectorD &bary_coords,
										  DenMatD &grad_bary_coords) {

	grad_b = EigVectorD::Zero(grad_bary_coords.cols());

	for (int i = 0; i < 4; ++i) {
		EigVectorD c;
		double B = 0;

		VectorI temp_alpha = alpha;
		temp_alpha[i] -= 1;

		if (temp_alpha[i] < 0) {
			continue;
		}

		c = n * grad_bary_coords.row(i);

		bernstein(B,
				  temp_alpha,
				  n - 1,
				  bary_coords);

		grad_b += c * B;
	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::omega_ij(EigVectorD &omega,
											VectorD &bary_coords,
											DenMatD &grad_bary_coords) {
	
	if(bary_coords.size() != 2 || grad_bary_coords.rows() != 2) {
		return FAILURE;
	}

	omega = bary_coords[0] * grad_bary_coords.row(1) - bary_coords[1] * grad_bary_coords.row(0);

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

	if (n == 0) {
		return FAILURE;
	}

	if (dim == 3) {
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

	else if (dim == 4) {
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

	if (n == 0) {
		return FAILURE;
	}

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

	S = S * get_simplex_volume(pts);
		
	return SUCCESS;
}


int FiniteElementExteriorCalculus::bb_mass_matrix_H_curl(DenMatD &mass_matrix,
													    Vector2D &pts,
													    int n,
													    Vector2I &alpha,
													    VectorI &ordered_basis_sizes) {

	if (n < 0) {
		return FAILURE;
	}

	Vector2I e;
	compute_index_sets_o(e,
						 1,
						 1);

	// compute_index_sets_o(alpha,
	// 					 2,
	// 					 2);
	// ordered_basis_sizes.push_back(alpha.size());

	// size_t total = alpha.size();
	// for (int i = 2; i <= 4; ++i) {
	// 	temp_alpha.clear();
	// 	compute_index_sets_o(temp_alpha,
	// 						 n + 1,
	// 						 i);
		
	// 	size_t temp_alpha_size = temp_alpha.size();
	// 	if (temp_alpha_size != 0) {
	// 		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	// 	}
	// 	ordered_basis_sizes.push_back(total + temp_alpha_size);
	// 	total += temp_alpha_size;
			
	// 	if (i == 3) {
	// 		temp_alpha.clear();
	// 		compute_index_sets_p(temp_alpha,
	// 							 n,
	// 							 i);
			
	// 		size_t temp_alpha_size = temp_alpha.size();
	// 		if (temp_alpha_size != 0) {
	// 			alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	// 		}
	// 		ordered_basis_sizes.push_back(total + temp_alpha_size);
	// 		total += temp_alpha_size;
	// 	}

	// 	else if (i == 4) {
	// 		for (int l = 0; l < 2; ++l) {	
	// 			temp_alpha.clear();
	// 			compute_index_sets_o(temp_alpha,
	// 								 n + 2,
	// 								 i);
				
	// 			size_t temp_alpha_size = temp_alpha.size();
	// 			if (temp_alpha_size != 0) {
	// 				alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	// 			}
	// 			ordered_basis_sizes.push_back(total + temp_alpha_size);
	// 			total += temp_alpha_size;
	// 		}

	// 		temp_alpha.clear();
	// 		compute_index_sets_o(temp_alpha,
	// 							 n + 2,
	// 							 i);
			
	// 		size_t temp_alpha_size = temp_alpha.size();

	// 		int counter = 0; 
	// 		for (size_t j = 0; j < temp_alpha_size; ++j) {
	// 			if (temp_alpha[j][2] == 1) {
	// 				alpha.push_back(temp_alpha[j]);
	// 				++counter;
	// 			}
	// 		}

	// 		ordered_basis_sizes.push_back(total + counter);	
	// 	}
	// }

	size_t size = alpha.size();
	mass_matrix.resize(size, size);

	DenMatD S_0;
	S_n(S_0,
		pts,
		0);

	Vector2I all_faces;
	compute_index_sets_o(all_faces,
						 3,
						 3);

	for (size_t i = 0; i < size; ++i) {
		for (size_t j = i; j < size; ++j) {

			double x = 0;
			
			// 33a
			if (i < ordered_basis_sizes[0] && j < ordered_basis_sizes[0]) {

				x = 0;

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
				
				double M_ji;
				M_alpha_beta(M_ji,
							 e[index_p[1]],
							 e[index_q[0]]);

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
				   ((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) ||
				   (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) || 
				   (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]))) {
				
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

				x = x * (n + 1);
			}

			// 33c
			else if(i < ordered_basis_sizes[0] && 
				    j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) {
				
				x = 0;

				int E_nF_size = ordered_basis_sizes[3] - ordered_basis_sizes[2];
				int face_index = std::floor((j - ordered_basis_sizes[2])/(E_nF_size/4));
				VectorI face = all_faces[face_index];
				
				Vector2I SIGMA;
				VectorI local_indices;

				for (int k = 0; k < 4; ++k) {
					if (face[k] > 0) {
						local_indices.push_back(k);
					}
				}

				for (int k = 0; k < 3; ++k) {
					VectorI vec {local_indices[k%3], local_indices[(k+1)%3], local_indices[(k+2)%3]};
					SIGMA.push_back(vec);
				}

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

			//33d
			else if (i < ordered_basis_sizes[0] &&
					 j >= ordered_basis_sizes[4]) {

				x = 0;

				int l_q;
				if (j >= ordered_basis_sizes[4] && j < ordered_basis_sizes[5]) {
					l_q = 0;
				}
				else if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
					l_q = 1;
				}
				if (j >= ordered_basis_sizes[6] && j < ordered_basis_sizes[7]) {
					l_q = 2;
				}

				for (int s = 0; s < 4; ++s) {
					int delta = 0;
					if (s == l_q) {
						delta = 1;
					}

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

					x += (delta * (n + 2) - alpha[j][s]) * ((M_i * S_j) - (M_j * S_i));
				}
			}

			// 33e
			else if(((i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) ||
				   (i >= ordered_basis_sizes[1] && i < ordered_basis_sizes[2]) || 
				   (i >= ordered_basis_sizes[3] && i < ordered_basis_sizes[4])) && 
				   ((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) ||
				   (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) || 
				   (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]))) {
				
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

				x = x * pow(n + 1, 2);
			}

			// 33f
			else if(((i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) ||
				   (i >= ordered_basis_sizes[1] && i < ordered_basis_sizes[2]) || 
				   (i >= ordered_basis_sizes[3] && i < ordered_basis_sizes[4])) && 
				   (j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3])) {
				
				x = 0;
				
				int E_nF_size = ordered_basis_sizes[3] - ordered_basis_sizes[2];
				int face_index = std::floor((j - ordered_basis_sizes[2])/(E_nF_size/4));
				VectorI face = all_faces[face_index];
				
				Vector2I SIGMA;
				VectorI local_indices;

				for (int k = 0; k < 4; ++k) {
					if (face[k] > 0) {
						local_indices.push_back(k);
					}
				}

				for (int k = 0; k < 3; ++k) {
					VectorI vec {local_indices[k%3], local_indices[(k+1)%3], local_indices[(k+2)%3]};
					SIGMA.push_back(vec);
				}

				size_t SIGMA_size = SIGMA.size();

				for (int t = 0; t < 4; ++t) {
					VectorI vec1;
					std::transform(alpha[i].begin(), alpha[i].end(), e[t].begin(),
										std::back_inserter(vec1),
										std::minus<int>());

					bool flag = false;
					for (int k = 0; k < 4; ++k) {
						if (vec1[k] < 0) {
							flag = true;
							break;
						}
					}

					if (flag) {
						continue;
					}

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

				x = (n + 1) * x;
			}

			// 33g
			else if(((i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) ||
				   (i >= ordered_basis_sizes[1] && i < ordered_basis_sizes[2]) || 
				   (i >= ordered_basis_sizes[3] && i < ordered_basis_sizes[4])) &&
				   (j >= ordered_basis_sizes[4])) {

				x = 0;

				int l_q;
				if (j >= ordered_basis_sizes[4] && j < ordered_basis_sizes[5]) {
					l_q = 0;
				}
				else if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
					l_q = 1;
				}
				else if (j >= ordered_basis_sizes[6] && j < ordered_basis_sizes[7]) {
					l_q = 2;
				}

				for (int t = 0; t < 4; ++t) {
					VectorI vec1;
					std::transform(alpha[i].begin(), alpha[i].end(), e[t].begin(),
										std::back_inserter(vec1),
										std::minus<int>());

					for (int s = 0; s < 4; ++s) {
						int delta = 0;
						if (s == l_q) {
							delta = 1;
						}

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

						x += (delta * (n + 2) - alpha[j][s]) * (M * S);
					}
				}

				x = x * (n + 1);
			}

			// 33h
			else if((i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) &&
				    (j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3])) {
				
				x = 0;
				
				Vector2I SIGMA_p;
				Vector2I SIGMA_q;
				VectorI local_indices_p;
				VectorI local_indices_q;

				int E_nF_size = ordered_basis_sizes[3] - ordered_basis_sizes[2];
				int face_index_p = std::floor((i - ordered_basis_sizes[2])/(E_nF_size/4));
				int face_index_q = std::floor((j - ordered_basis_sizes[2])/(E_nF_size/4));
				VectorI face_p = all_faces[face_index_p];
				VectorI face_q = all_faces[face_index_q];

				for (int k = 0; k < 4; ++k) {
					if (face_p[k] > 0) {
						local_indices_p.push_back(k);
					}
					if (face_q[k] > 0) {
						local_indices_q.push_back(k);
					}
				}

				for (int k = 0; k < 3; ++k) {
					VectorI vec_p {local_indices_p[k%3], local_indices_p[(k+1)%3], local_indices_p[(k+2)%3]};
					VectorI vec_q {local_indices_q[k%3], local_indices_q[(k+1)%3], local_indices_q[(k+2)%3]};
					SIGMA_p.push_back(vec_p);
					SIGMA_q.push_back(vec_q);
				}

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

			// 33i
			else if((i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) &&
				    (j >= ordered_basis_sizes[4])) {

				x = 0;

				int l_q;
				if (j >= ordered_basis_sizes[4] && j < ordered_basis_sizes[5]) {
					l_q = 0;
				}
				else if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
					l_q = 1;
				}
				else if (j >= ordered_basis_sizes[6] && j < ordered_basis_sizes[7]) {
					l_q = 2;
				}
				
				Vector2I SIGMA_p;
				VectorI local_indices_p;

				int E_nF_size = ordered_basis_sizes[3] - ordered_basis_sizes[2];
				int face_index_p = std::floor((i - ordered_basis_sizes[2])/(E_nF_size/4));
				VectorI face_p = all_faces[face_index_p];

				for (int k = 0; k < 4; ++k) {
					if (face_p[k] > 0) {
						local_indices_p.push_back(k);
					}
				}

				for (int k = 0; k < 3; ++k) {
					VectorI vec_p {local_indices_p[k%3], local_indices_p[(k+1)%3], local_indices_p[(k+2)%3]};
					SIGMA_p.push_back(vec_p);
				}

				size_t SIGMA_p_size = SIGMA_p.size();

				for (int t = 0; t < SIGMA_p_size; ++t) {
					VectorI vec1;
					std::transform(alpha[i].begin(), alpha[i].end(), e[SIGMA_p[t][0]].begin(),
										std::back_inserter(vec1),
										std::plus<int>());

					for (int s = 0; s < 4; ++s) {
						int delta = 0;
						if (s == l_q) {
							delta = 1;
						}

						VectorI vec2;
						std::transform(alpha[j].begin(), alpha[j].end(), e[s].begin(),
											std::back_inserter(vec2),
											std::minus<int>());

						double M;
						M_alpha_beta(M,
									 vec1,
									 vec2);

						x += M * (delta * (n + 2) - alpha[j][s]) * (alpha[i][SIGMA_p[t][0]] + 1)
						       * (alpha[i][SIGMA_p[t][2]] * S_0.coeffRef(SIGMA_p[t][1], s)
						       	  - alpha[i][SIGMA_p[t][1]] * S_0.coeffRef(SIGMA_p[t][2], s));
					}
				}
			}

			// 33j
			else if(i >= ordered_basis_sizes[4] &&
				    j >= ordered_basis_sizes[4]) {

				x = 0;

				int l_q;
				int l_p;
				if (i >= ordered_basis_sizes[4] && i < ordered_basis_sizes[5]) {
					l_p = 0;
				}
				else if (i >= ordered_basis_sizes[5] && i < ordered_basis_sizes[6]) {
					l_p = 1;
				}
				else if (i >= ordered_basis_sizes[6] && i < ordered_basis_sizes[7]) {
					l_p = 2;
				}
				if (j >= ordered_basis_sizes[4] && j < ordered_basis_sizes[5]) {
					l_q = 0;
				}
				else if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
					l_q = 1;
				}
				else if (j >= ordered_basis_sizes[6] && j < ordered_basis_sizes[7]) {
					l_q = 2;
				}

				for (int t = 0; t < 4; ++t) {
					int delta_t = 0;
					if (l_p == t) {
						delta_t = 1;
					}

					VectorI vec1;
					std::transform(alpha[i].begin(), alpha[i].end(), e[t].begin(),
										std::back_inserter(vec1),
										std::minus<int>());

					for (int s = 0; s < 4; ++s) {
						int delta_s = 0;
						if (l_q == s) {
							delta_s = 1;
						}
						
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

						x += M * (delta_t * (n + 2) - alpha[i][t]) * (delta_s * (n + 2) - alpha[j][s]) * S;
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



int FiniteElementExteriorCalculus::bb_mass_matrix_H_1(DenMatD &mass_matrix,
													  int n,
													  int m,
													  Vector2I &index_sets,
													  int d) {

	int max = std::max(n, m);

	if (max < 0) {
		return FAILURE;
	}

	// Vector2I index_sets;
	// Vector2I temp_index_sets;
	// compute_index_sets_o(index_sets,
	// 					 1,
	// 					 1);

	// for (int i = 2; i <= std::min(max, d+1); ++i) {
	// 	temp_index_sets.clear();
	// 	compute_index_sets_o(temp_index_sets,
	// 						  max,
	// 						  i);
	// 	index_sets.insert(index_sets.end(), temp_index_sets.begin(), temp_index_sets.end());
	// }

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

double FiniteElementExteriorCalculus::bb_error_H_curl_1d_quad(int n,
														      Vector3I &simplices,
															  Vector2D &vertices,
															  VectorI &num_simplices,
															  int q_order) {

	size_t N = num_simplices.size();
	double E = 0.0;
	size_t embed_dim = vertices[0].size();

	Vector2D nodes;
	VectorD weights;
	Vector2D nodes_1d;
	VectorD weights_1d;

	std::string data = "./data/quadrature/d" + std::to_string(N-1) + "o" + std::to_string(q_order) + ".txt";
	std::string data_1d = "./data/quadrature/d1o" + std::to_string(q_order) + ".txt";
	
	read_quadratures(nodes,
					 weights,
					 data);
	size_t nodes_size = nodes.size();

	read_quadratures(nodes_1d,
					 weights_1d,
					 data_1d);
	size_t nodes_1d_size = nodes_1d.size();

	Vector2I alpha;
	Vector2I temp_alpha;

	compute_index_sets_o(alpha,
					 	 2,
					 	 2);

	int d = alpha[0].size();

	// for(int i = 2; i <= std::min(n, d); ++i) {
	// 	temp_alpha.clear();
	// 	compute_index_sets_o(temp_alpha,
	// 					 	 n,
	// 					 	 i);

	// 	alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	// }
	size_t alpha_size = alpha.size();

	// #ifdef MULTICORE
	// 	#pragma omp parallel for
	// #endif
	// for(size_t i = 0; i < num_simplices[N-1]; ++i) {
	for(size_t i = 0; i < 1; ++i) {
		double e = 0;

		Vector2D pts;
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplices[N-1][i][k]]);
		}
		double vol = get_simplex_volume(pts);

		Vector3D edges_pts;
		VectorI edge_indices = simplex_sub_simplices[i][1];
		size_t edge_indices_size = edge_indices.size();

		for(size_t j = 0; j < edge_indices_size; ++j) {
			Vector2D temp_edges_pts;
			for(size_t k = 0; k < 2; ++k) {
				temp_edges_pts.push_back(vertices[simplices[1][edge_indices[j]][k]]);
			}
			edges_pts.push_back(temp_edges_pts);
		}

		DenMatD grad_bary_coords;
		barycentric_gradients(grad_bary_coords,
							  pts);

		double sum_weight = 0.0;

		EigVectorD b(alpha_size);

		for(size_t j = 0; j < edge_indices_size; ++j) {
			double trace = 0;

			EigVectorD v0(embed_dim);
			EigVectorD v1(embed_dim);
			EigVectorD func(embed_dim);

			for (size_t v = 0; v < embed_dim; ++v) {
				v0.coeffRef(v) = edges_pts[j][0][v];
				v1.coeffRef(v) = edges_pts[j][1][v];
			}

			double sum_weights = 0;

			for (size_t node_1d = 0; node_1d < nodes_1d_size; ++node_1d) {
				EigVectorD vec1;
				VectorD vec2;
				VectorD temp_f;

				vec1 = (v0+v1)/2 + nodes_1d[node_1d][0] * (v1 - v0)/2;

				for (size_t v = 0; v < embed_dim; ++v) {
					vec2.push_back(vec1.coeffRef(v));
				}
				
				get_analytical_soln_vec(temp_f,
										vec2);

				for (size_t v = 0; v < embed_dim; ++v) {
					func.coeffRef(v) = temp_f[v];
				}

				trace += weights_1d[node_1d] * func.dot(v1 - v0);
				sum_weights += weights_1d[node_1d];
			}
			
			b(j) = trace/sum_weights;
		}

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {
			DenMatD omega(alpha_size, embed_dim);
			EigVectorD f(embed_dim);
			EigVectorD f_dash(embed_dim);

			for (size_t k = 0; k < alpha_size; ++k) {
				EigVectorD temp_omega;
				VectorD temp_bary_coords;
				DenMatD temp_grad_bary_coords(2, embed_dim);
				int c = 0;
				for (size_t j = 0; j < d; ++j) {
					if (alpha[k][j] > 0) {
						temp_bary_coords.push_back(nodes[node_index][j]);
						temp_grad_bary_coords.row(c) = grad_bary_coords.row(j);
						++c;
					}
				}

				omega_ij(temp_omega,
						 temp_bary_coords,
						 temp_grad_bary_coords);

				omega.row(k) = temp_omega;
			}

			VectorD points(embed_dim, 0.0);
			for(size_t v = 0; v < N; ++v) {
				for(size_t l = 0; l < embed_dim; ++l) {
					points[l] += pts[v][l] * nodes[node_index][v];
				}
			}

			f_dash = b.transpose() * omega;

			VectorD F;
			get_analytical_soln_vec(F,
									points);

			for (size_t j = 0; j < embed_dim; ++j) {
				f.coeffRef(j) = F[j];
			}

			e += weights[node_index] * pow((f-f_dash).norm(), 2);
			sum_weight += weights[node_index];
		}

		#ifdef MULTICORE
			#pragma omp critical
		#endif
		E += vol*e/sum_weight;
	}

	return E;
}


double FiniteElementExteriorCalculus::bb_error_H_1(int n,
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

	double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);

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

	DenMatD M;
	bb_mass_matrix_H_1(M,
					   n,
					   n,
					   alpha);

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

		DenMatD temp_M = M * vol;

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {
			EigVectorD b(alpha_size);

			for(size_t j = 0; j < alpha_size; ++j) {
				double inner_product = 0.0;

				for(size_t k = 0; k < nodes_size; ++k) {
					VectorD vec(embed_dim, 0.0);

					for(size_t v = 0; v < N; ++v) {
						for(size_t l = 0; l < embed_dim; ++l) {
							vec[l] += pts[v][l] * nodes[k][v];
						}
					}

					inner_product += vol * weights[k] * get_analytical_soln(vec) * basis_elements.coeffRef(j, k);
				}

				b.coeffRef(j) = inner_product/sum_weights;
			}

			EigVectorD coeffs = temp_M.colPivHouseholderQr().solve(b);

			double f_dash = coeffs.dot(basis_elements.col(node_index));
			
			VectorD points(embed_dim, 0.0);
			for(size_t v = 0; v < N; ++v) {
				for(size_t l = 0; l < embed_dim; ++l) {
					points[l] += pts[v][l] * nodes[node_index][v];
				}
			}

			e += weights[node_index] * pow(get_analytical_soln(points) - f_dash, 2);
		}

		#ifdef MULTICORE
			#pragma omp critical
		#endif
		E += vol*e/sum_weights;
	}

	E = sqrt(E);
	return E;
}


double FiniteElementExteriorCalculus::bb_error_H_curl(int n,
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

	double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);

	Vector2I alpha;
	Vector2I temp_alpha;
	VectorI ordered_basis_sizes;
	Vector2D basis_vector;

	compute_index_sets_o(alpha,
						 2,
						 2);
	ordered_basis_sizes.push_back(alpha.size());

	size_t total = alpha.size();
	for (int i = 2; i <= 4; ++i) {
		temp_alpha.clear();
		compute_index_sets_o(temp_alpha,
							 n + 1,
							 i);
		
		size_t temp_alpha_size = temp_alpha.size();
		if (temp_alpha_size != 0) {
			alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
		}
		ordered_basis_sizes.push_back(total + temp_alpha_size);
		total += temp_alpha_size;
			
		if (i == 3) {
			temp_alpha.clear();
			compute_index_sets_p(temp_alpha,
								 n,
								 i);
			
			size_t temp_alpha_size = temp_alpha.size();
			if (temp_alpha_size != 0) {
				alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
			}
			ordered_basis_sizes.push_back(total + temp_alpha_size);
			total += temp_alpha_size;
		}

		else if (i == 4) {
			for (int l = 0; l < 2; ++l) {	
				temp_alpha.clear();
				compute_index_sets_o(temp_alpha,
									 n + 2,
									 i);
				
				size_t temp_alpha_size = temp_alpha.size();
				if (temp_alpha_size != 0) {
					alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
				}
				ordered_basis_sizes.push_back(total + temp_alpha_size);
				total += temp_alpha_size;
			}

			temp_alpha.clear();
			compute_index_sets_o(temp_alpha,
								 n + 2,
								 i);
			
			size_t temp_alpha_size = temp_alpha.size();

			size_t counter = 0; 
			for (size_t j = 0; j < temp_alpha_size; ++j) {
				if (temp_alpha[j][2] == 1) {
					alpha.push_back(temp_alpha[j]);
					++counter;
				}
			}

			ordered_basis_sizes.push_back(total + counter);	
		}
	}

	size_t alpha_size = alpha.size();

	Vector2I all_faces;
	compute_index_sets_o(all_faces,
						 3,
						 3);

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

		DenMatD grad_bary_coords;
		barycentric_gradients(grad_bary_coords,
							  pts);

		VectorDenMatD basis_elements;
		for(size_t i = 0; i < alpha_size; ++i) {
			DenMatD temp_basis_elements(nodes_size, embed_dim);

			VectorI local_indices;
			for (int j = 0; j < 4; ++j) {
				if (alpha[i][j] > 0) {
					local_indices.push_back(j);
				}
			}
			size_t local_indices_size = local_indices.size();
			
			for(size_t j = 0; j < nodes_size; ++j) {
				if (i < ordered_basis_sizes[0]) {
					EigVectorD omega;

					DenMatD temp_grad_bary_coords(local_indices_size, embed_dim);
					VectorD temp_bary_coords;

					for (size_t k = 0; k < local_indices_size; ++k) {
						temp_bary_coords.push_back(nodes[j][local_indices[k]]);
						temp_grad_bary_coords.row(k) = grad_bary_coords.row(local_indices[k]);
					}
					
					omega_ij(omega,
							 temp_bary_coords,
							 temp_grad_bary_coords);
					
					temp_basis_elements.row(j) = omega;
				}
				else if ((i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) ||
						 (i >= ordered_basis_sizes[1] && i < ordered_basis_sizes[2]) ||
						 (i >= ordered_basis_sizes[3] && i < ordered_basis_sizes[4])) {
					EigVectorD grad_b;

					grad_B(grad_b,
						   alpha[i],
						   n + 1,
						   nodes[j],
						   grad_bary_coords);

					temp_basis_elements.row(j) = grad_b;
				}
				else if (i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) {
					EigVectorD phi;

					int E_nF_size = ordered_basis_sizes[3] - ordered_basis_sizes[2];
					int face_index = std::floor((i - ordered_basis_sizes[2])/(E_nF_size/4));
					VectorI face = all_faces[face_index];
					
					VectorI temp_local_indices;

					for (int k = 0; k < 4; ++k) {
						if (face[k] > 0) {
							temp_local_indices.push_back(k);
						}
					}

					phi_FT(phi,
						   alpha[i],
						   n,
						   nodes[j],
						   grad_bary_coords,
						   temp_local_indices);

					temp_basis_elements.row(j) = phi;

				}
				else if (i >= ordered_basis_sizes[4] && i < ordered_basis_sizes[5]) {
					EigVectorD psi;

					psi_T(psi,
						   alpha[i],
						   n + 1,
						   0,
						   nodes[j],
						   grad_bary_coords);

					temp_basis_elements.row(j) = psi;
				}
				else if (i >= ordered_basis_sizes[5] && i < ordered_basis_sizes[6]) {
					EigVectorD psi;

					psi_T(psi,
						   alpha[i],
						   n + 1,
						   1,
						   nodes[j],
						   grad_bary_coords);

					temp_basis_elements.row(j) = psi;
				}
				else if (i >= ordered_basis_sizes[6] && i < ordered_basis_sizes[7]) {
					EigVectorD psi;

					psi_T(psi,
						   alpha[i],
						   n + 1,
						   2,
						   nodes[j],
						   grad_bary_coords);

					temp_basis_elements.row(j) = psi;
				}
			}

			basis_elements.push_back(temp_basis_elements);
		}

		DenMatD M;
		bb_mass_matrix_H_curl(M,
							  pts,
						 	  n,
						 	  alpha,
						 	  ordered_basis_sizes);

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {
			EigVectorD b(alpha_size);

			for(size_t j = 0; j < alpha_size; ++j) {
				double inner_product = 0.0;

				for(size_t k = 0; k < nodes_size; ++k) {
					VectorD vec(embed_dim, 0.0);

					for(size_t v = 0; v < N; ++v) {
						for(size_t l = 0; l < embed_dim; ++l) {
							vec[l] += pts[v][l] * nodes[k][v];
						}
					}

					VectorD temp_vec;
					get_analytical_soln_vec(temp_vec,
											vec);
					EigVectorD f(embed_dim);
					for (size_t v = 0; v < embed_dim; ++v) {
						f.coeffRef(v) = temp_vec[v];
					}
					inner_product += vol * weights[k] * f.dot(basis_elements[j].row(k));
				}

				b.coeffRef(j) = inner_product/sum_weights;
			}

			// std::cout<<"b\n"<<b;

			EigVectorD coeffs = M.colPivHouseholderQr().solve(b);

			// std::cout<<"coeffs\n"<<coeffs;

			EigVectorD f_dash = EigVectorD::Zero(embed_dim);
			for (size_t j = 0; j < alpha_size; ++j) {
				f_dash += coeffs.coeffRef(j) * basis_elements[j].row(node_index);
			}

			// std::cout<<"f_dash\n"<<f_dash;

			VectorD points(embed_dim, 0.0);
			for(size_t v = 0; v < N; ++v) {
				for(size_t l = 0; l < embed_dim; ++l) {
					points[l] += pts[v][l] * nodes[node_index][v];
				}
			}

			VectorD temp_vec;
			get_analytical_soln_vec(temp_vec,
									points);
			EigVectorD f(embed_dim);
			for (size_t v = 0; v < embed_dim; ++v) {
				f.coeffRef(v) = temp_vec[v];
			}

			e += weights[node_index] * pow((f - f_dash).norm(), 2);
		}

		#ifdef MULTICORE
			#pragma omp critical
		#endif
		E += vol*e/sum_weights;
	}

	E = sqrt(E);
	return E;
}


FiniteElementExteriorCalculus::FiniteElementExteriorCalculus() {

    set_hodge_stars_to_null();
    
}


FiniteElementExteriorCalculus::FiniteElementExteriorCalculus(SimplicialComplex sc) : GeometryComplex(sc) {

    set_hodge_stars_to_null();
}

FiniteElementExteriorCalculus::~FiniteElementExteriorCalculus() {}