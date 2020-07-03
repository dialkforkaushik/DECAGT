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
#include<Eigen/SparseCholesky>
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


int set_increasing_ordering1(Vector2I &sets) {

	size_t s1 = sets.size();
	size_t s2 = sets[0].size();

	std::vector < std::pair <int, VectorI> > sorted_pair; 
	for(size_t i = 0; i < s1; ++i) {
		int k = 0;
		for(size_t j = 0; j < s2; ++j) {
			k = 10 * k + sets[i][j];
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
    	VectorI temp = VectorI(p.begin(), p.begin() + std::min(k + 1, d + 1));
    	int s = std::accumulate(temp.begin(), temp.end(), 0);

        if ((k + 1 == dim && s == n) || (dim == 0 && s == n)) {
        	VectorI sliced = VectorI(p.begin(), p.begin() + std::min(k + 1, d + 1));
        	for (int i = std::min(k + 1, d + 1); i < d + 1; ++i) {
        		sliced.push_back(0);
        	}
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


int t_ij(EigVectorD &t,
		int i,
		int j,
		DenMatD &grad_bary_coords) {

	EigVectorD v1 = grad_bary_coords.row(i);
	EigVectorD v2 = grad_bary_coords.row(j);

	cross_product(t,
				  v1,
				  v2);

	return SUCCESS;
}


int e_ijk(double &e,
		  int i,
		  int j, 
		  int k, 
		  DenMatD &grad_bary_coords) {

	EigVectorD t;
	t_ij(t,
		 i,
		 j,
		 grad_bary_coords);

	e = t.dot(grad_bary_coords.row(k));

	return SUCCESS;
}


int S_ij_kl(double &S,
		   int i,
		   int j,
		   int k,
		   int l,
		   DenMatD &S_1,
		   Vector2I &all_edges) {

	if (i == j || k == l) {
		S = 0;
		return SUCCESS;
	}

	size_t all_edges_size = all_edges.size();

	int index1 = -1;
	int index2 = -1;
	int sign = copysign(1, (j - i) * (l - k));

	for (size_t e = 0; e < all_edges_size; ++e) {
		if (all_edges[e][i] > 0 && all_edges[e][j] > 0) {
			index1 = e;
		}
		if (all_edges[e][k] > 0 && all_edges[e][l] > 0) {
			index2 = e;
		}
	}

	S = S_1.coeffRef(index1, index2) * sign;

	return SUCCESS;
}


int get_global_index(int &index,
					 int k,
					 int i,
					 int j,
					 VectorI &ordered_basis_sizes,
					 VectorI &ndofs,
					 VectorI &sizes,
					 Vector3I &simplex_sub_simplices,
					 int complex_dimension) {

	int num_edges;
	binomialCoeff(num_edges,
				  complex_dimension + 1,
				  2);
	int num_faces;
	binomialCoeff(num_faces,
				  complex_dimension + 1,
				  3);

	if (k == 0) {
		if (j < ordered_basis_sizes[0]) {
			index = simplex_sub_simplices[i][0][j];
		}
		else if (j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) {
			int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
			int temp_index1 = std::floor((j - ordered_basis_sizes[0])/(temp/num_edges));
			index = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index1] + (j - ordered_basis_sizes[0])%ndofs[1];
		}
		else if (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) {
			int temp1 = ordered_basis_sizes[2] - ordered_basis_sizes[1];
			int temp_index1 = std::floor((j - ordered_basis_sizes[1])/(temp1/num_faces));
			index = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[1])%ndofs[2];
		}
		else if (j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) {
			int temp1 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
			int temp_index1 = std::floor((j - ordered_basis_sizes[2])/(temp1));
			index = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[2])%ndofs[3];
		}
	}

	else if (k == 1) {
		if (j < ordered_basis_sizes[0]) {
			index = simplex_sub_simplices[i][1][j];
		}
		
		else if (((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) ||
				  (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) ||
				  (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]))) {

			if (j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) {
				int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
				int temp_index = std::floor((j - ordered_basis_sizes[0])/(temp/num_edges));
				index = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index] + (j - ordered_basis_sizes[0])%ndofs[1];
			}
			else if (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) {
				int temp = ordered_basis_sizes[2] - ordered_basis_sizes[1];
				int temp_index = std::floor((j - ordered_basis_sizes[1])/(temp/num_faces));
				index = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index] + (j - ordered_basis_sizes[1])%ndofs[2];
			}
			else if (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]) {
				int temp = ordered_basis_sizes[4] - ordered_basis_sizes[3];
				int temp_index = std::floor((j - ordered_basis_sizes[3])/(temp));
				index = sizes[3] + ndofs[4] * simplex_sub_simplices[i][3][temp_index] + (j - ordered_basis_sizes[3])%ndofs[4];
			}
		}
		else if (j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) {
			int temp1 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
			int temp_index1 = std::floor((j - ordered_basis_sizes[2])/(temp1/num_faces));
			index = sizes[2] + ndofs[3] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[2])%ndofs[3];
		}
		else if (j >= ordered_basis_sizes[4]) {
			int temp1 = ordered_basis_sizes[7] - ordered_basis_sizes[4];
			int temp_index1 = std::floor((j - ordered_basis_sizes[4])/(temp1));
			index = sizes[4] + ndofs[5] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[4])%ndofs[5];
		}
	}

	else if (k == 2) {
		if (j < ordered_basis_sizes[0]) {
			index = simplex_sub_simplices[i][2][j];
		}
		else if (j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) {
			int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
			int temp_index1 = std::floor((j - ordered_basis_sizes[0])/(temp/num_faces));
			index = sizes[0] + ndofs[1] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[0])%ndofs[1];
		}
		else if (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[4]) {
			int temp1 = ordered_basis_sizes[4] - ordered_basis_sizes[1];
			int temp_index1 = std::floor((j - ordered_basis_sizes[1])/(temp1));
			index = sizes[1] + ndofs[2] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[1])%ndofs[2];
		}
		else if (j >= ordered_basis_sizes[4] && j < ordered_basis_sizes[5]) {
			int temp1 = ordered_basis_sizes[5] - ordered_basis_sizes[4];
			int temp_index1 = std::floor((j - ordered_basis_sizes[4])/(temp1));
			index = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[4])%ndofs[3];
		}
	}

 	return SUCCESS;
}


int FiniteElementExteriorCalculus::upsilon(EigVectorD &upsi,
										  VectorI &alpha,
										  int n,
										  VectorD &bary_coords,
										  DenMatD &grad_bary_coords) {

	upsi = EigVectorD::Zero(grad_bary_coords.cols());

	for (int i = 0; i < 4; ++i) {
		EigVectorD d = EigVectorD::Zero(grad_bary_coords.cols());

		VectorI local_indices;
		if (i == 0) {
			local_indices = {0, 1, 2};	
		}
		else if (i == 1) {
			local_indices = {0, 1, 3};	
		}
		else if (i == 2) {
			local_indices = {0, 2, 3};	
		}
		else if (i == 3) {
			local_indices = {1, 2, 3};	
		}

		Vector2I SIGMA;
		for (int k = 0; k < 3; ++k) {
			VectorI vec {local_indices[k%3], local_indices[(k+1)%3], local_indices[(k+2)%3]};
			SIGMA.push_back(vec);
		}
		size_t SIGMA_size = SIGMA.size();

		for (size_t k = 0; k < SIGMA_size; ++k) {
			VectorI sigma = SIGMA[k];
			EigVectorD t;
			t_ij(t,
				 sigma[1],
				 sigma[2],
				 grad_bary_coords);
			d += alpha[sigma[0]] * t;
		}

		d *= pow(-1, i) * (alpha[i] + 1);

		double B = 0;
		VectorI temp_alpha = alpha;
		temp_alpha[i] += 1;

		bernstein(B,
				  temp_alpha,
				  n + 1,
				  bary_coords);

		upsi += d * B;
	}

	return SUCCESS;
}



int FiniteElementExteriorCalculus::curl_psi_T(EigVectorD &curl_psi,
											  VectorI &alpha,
											  int n,
											  int l,
											  VectorD &bary_coords,
											  DenMatD &grad_bary_coords) {

	curl_psi = EigVectorD::Zero(grad_bary_coords.cols());

	for (int i = 0; i < 4; ++i) {
		EigVectorD d;
		double B = 0;

		VectorI temp_alpha = alpha;
		temp_alpha[l] -= 1;
		temp_alpha[i] -= 1;

		if (temp_alpha[l] < 0 || temp_alpha[i] < 0) {
			continue;
		}

		EigVectorD t;
		t_ij(t,
			 i,
			 l,
			 grad_bary_coords);

		d = n * (n + 1) * t;

		bernstein(B,
				  temp_alpha,
				  n - 1,
				  bary_coords);

		curl_psi += d * B;
	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::curl_phi_FT(EigVectorD &curl_phi,
										  	   VectorI &alpha,
										  	   int n,
											   VectorD &bary_coords,
											   DenMatD &grad_bary_coords,
											   VectorI &local_indices) {

	size_t size = alpha.size();

	if(size != bary_coords.size()) {
		return FAILURE;
	}

	Vector2I SIGMA;

	for (int k = 0; k < 3; ++k) {
		VectorI vec {local_indices[k%3], local_indices[(k+1)%3], local_indices[(k+2)%3]};
		SIGMA.push_back(vec);
	}

	size_t SIGMA_size = SIGMA.size();

	curl_phi = EigVectorD::Zero(grad_bary_coords.cols());

	for (size_t i = 0; i < SIGMA_size; ++i) {
		VectorI sigma = SIGMA[i];

		for (int k = 0; k < 3; ++k) {
			double B = 0;
			VectorI temp_alpha = alpha;
			temp_alpha[sigma[0]] += 1;
			temp_alpha[sigma[k]] -= 1;

			if (temp_alpha[sigma[k]] < 0) {
				continue;
			}

			EigVectorD t1;
			t_ij(t1,
				 sigma[k],
				 sigma[1],
				 grad_bary_coords);
			EigVectorD t2;
			t_ij(t2,
				 sigma[k],
				 sigma[2],
				 grad_bary_coords);

			EigVectorD d = (n + 1) * (alpha[sigma[0]] + 1) * (alpha[sigma[2]] * t1 - alpha[sigma[1]] * t2);

			bernstein(B,
					  temp_alpha,
					  n,
					  bary_coords);

			curl_phi += d * B;
		}
	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::chi_l(EigVectorD &chi,
										 VectorI &alpha,
										 VectorD &bary_coords,
										 DenMatD &grad_bary_coords,
										 VectorI &local_indices) {

	size_t size = alpha.size();

	if(size != bary_coords.size()) {
		return FAILURE;
	}

	Vector2I SIGMA;

	for (int k = 0; k < 3; ++k) {
		VectorI vec {local_indices[k%3], local_indices[(k+1)%3], local_indices[(k+2)%3]};
		SIGMA.push_back(vec);
	}

	size_t SIGMA_size = SIGMA.size();

	Vector2I e;
	compute_index_sets_o(e,
						 1,
						 1);

	chi = EigVectorD::Zero(grad_bary_coords.cols());

	for (size_t i = 0; i < SIGMA_size; ++i) {
		VectorI sigma = SIGMA[i];

		EigVectorD t;
		EigVectorD temp1 = grad_bary_coords.row(sigma[1]);
		EigVectorD temp2 = grad_bary_coords.row(sigma[2]);
		cross_product(t,
					  temp1,
					  temp2);

		double B = 0;
		bernstein(B,
				  e[sigma[0]],
				  1,
				  bary_coords);

		chi += t * B;
	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::curl_omega_ij(EigVectorD &curl_omega,
												 int i,
												 int j,
												 DenMatD &grad_bary_coords) {
	EigVectorD t;
	t_ij(t,
		 i, 
		 j,
		 grad_bary_coords);

	curl_omega = 2 * t;

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

	// print_vector(local_indices);

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

	// if(size != bary_coords.size()) {
	// 	return FAILURE;
	// }

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

  		set_increasing_ordering1(sets);
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


int FiniteElementExteriorCalculus::set_mass_matrices_to_null() {

    for (int i = 0; i < num_simplices[complex_dimension]; ++i) {
    	auto x = std::numeric_limits<SpMatD>::quiet_NaN();
    	mass_matrices.push_back(x);
    }

    return SUCCESS;
}


int FiniteElementExteriorCalculus::set_bb_mass_matrices_to_null() {

    for (int i = 0; i < 4; ++i) {
    	auto x = std::numeric_limits<SpMatD>::quiet_NaN();
    	bb_mass_matrices.push_back(x);
    }

    return SUCCESS;
}


int FiniteElementExteriorCalculus::compute_mass_matrix_k(int &k) {
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
	mass_matrices[k] = mass_matrix;


	return SUCCESS;
}


int FiniteElementExteriorCalculus::compute_mass_matrices() {
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
		mass_matrices[k] = mass_matrix;

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

	DenMatD gradients;
	barycentric_gradients(gradients,
						  pts);

	if (n == 0) {

		size_t gradients_size = gradients.rows();
		S.resize(gradients_size, gradients_size);

		for (size_t i = 0; i < gradients_size; ++i) {
			for (size_t j = i; j < gradients_size; ++j) {
				S.coeffRef(i,j) = gradients.row(i).dot(gradients.row(j));
				if (i != j) {
					S.coeffRef(j, i) = S.coeffRef(i, j);
				}
			}
		}
	}

	else if (n == 1) {
		Vector2I edge_indices;
		compute_index_sets_o(edge_indices,
							 2,
							 2);
		size_t edge_indices_size = edge_indices.size();

		S.resize(edge_indices_size, edge_indices_size);

		for (size_t i = 0; i < edge_indices_size; ++i) {
			EigVectorD t1;
			VectorEigVectorD temp_grad1;
			for (size_t k = 0; k < 4; ++k) {
				if (edge_indices[i][k] > 0) {
					temp_grad1.push_back(gradients.row(k));
				}
			}
			cross_product(t1,
						  temp_grad1[0],
						  temp_grad1[1]);

			for (size_t j = i; j < edge_indices_size; ++j) {
				EigVectorD t2;
				VectorEigVectorD temp_grad2;
				for (size_t k = 0; k < 4; ++k) {
					if (edge_indices[j][k] > 0) {
						temp_grad2.push_back(gradients.row(k));
					}
				}
				cross_product(t2,
							  temp_grad2[0],
							  temp_grad2[1]);

				S.coeffRef(i, j) = t1.dot(t2);
				if (i != j) {
					S.coeffRef(j, i) = S.coeffRef(i, j);
				}
			}
		}
	}

	double vol = get_simplex_volume(pts);
	S = S * vol;
		
	return SUCCESS;
}


int FiniteElementExteriorCalculus::bb_mass_matrix_H_div(DenMatD &mass_matrix,
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

	size_t size = alpha.size();
	mass_matrix.resize(size, size);

	DenMatD S_1;
	S_n(S_1,
		pts,
		1);

	Vector2I all_edges;
	compute_index_sets_o(all_edges,
						 2,
						 2);
	size_t all_edges_size = all_edges.size();

	Vector2I all_faces;
	compute_index_sets_o(all_faces,
						 3,
						 3);

	for (int i = 0; i < size; ++i) {
		for (int j = i; j < size; ++j) {

			double x = 0;
			
			// 36a
			if (i < ordered_basis_sizes[0] && j < ordered_basis_sizes[0]) {

				x = 0;

				Vector2I SIGMA_p;
				Vector2I SIGMA_q;
				VectorI local_indices_p;
				VectorI local_indices_q;

				for (int k = 0; k < 4; ++k) {
					if (alpha[i][k] > 0) {
						local_indices_p.push_back(k);
					}
					if (alpha[j][k] > 0) {
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

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (size_t s = 0; s < SIGMA_q_size; ++s) {

						double M;
						M_alpha_beta(M,
									 e[SIGMA_p[t][0]],
									 e[SIGMA_q[s][0]]);

						double S;
						S_ij_kl(S,
								SIGMA_p[t][1],
								SIGMA_p[t][2],
								SIGMA_q[s][1],
								SIGMA_q[s][2],
								S_1,
								all_edges);

						x += M * S;
					}
				}	
			}

			//36b
			else if (i < ordered_basis_sizes[0] &&
					(j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1])) {

				x = 0;

				Vector2I SIGMA_p;
				Vector2I SIGMA_q;
				VectorI local_indices_p;
				VectorI local_indices_q;

				for (int k = 0; k < 4; ++k) {
					if (alpha[i][k] > 0) {
						local_indices_p.push_back(k);
					}
				}

				int E_nF_size = ordered_basis_sizes[1] - ordered_basis_sizes[0];
				int face_index = std::floor((j - ordered_basis_sizes[0])/(E_nF_size/4));
				VectorI face = all_faces[face_index];

				for (int k = 0; k < 4; ++k) {
					if (face[k] > 0) {
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

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (size_t s = 0; s < SIGMA_q_size; ++s) {
						for (int k = 0; k < 3; ++k) {
							VectorI temp_alpha = alpha[j];
							temp_alpha[SIGMA_q[s][0]] += 1;
							temp_alpha[SIGMA_q[s][k]] -= 1;

							if (temp_alpha[SIGMA_q[s][k]] < 0) {
								continue;
							}

							double M;
							M_alpha_beta(M,
										 e[SIGMA_p[t][0]],
										 temp_alpha);

							double S1;
							S_ij_kl(S1,
									SIGMA_p[t][1],
									SIGMA_p[t][2],
									SIGMA_q[s][k],
									SIGMA_q[s][1],
									S_1,
									all_edges);
							double S2;
							S_ij_kl(S2,
									SIGMA_p[t][1],
									SIGMA_p[t][2],
									SIGMA_q[s][k],
									SIGMA_q[s][2],
									S_1,
									all_edges);

							x += M * (alpha[j][SIGMA_q[s][0]] + 1) 
							       * (alpha[j][SIGMA_q[s][2]] * S1 - alpha[j][SIGMA_q[s][1]] * S2);
						}
					}
				}

				x = x * (n + 1);
			}

			//36c
			else if (i < ordered_basis_sizes[0] && 
					(j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[4])) {

				x = 0;

				int l_q = 0;
				if (j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) {
					l_q = 1;
				}
				else if (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]) {
					l_q = 2;
				}

				Vector2I SIGMA_p;
				VectorI local_indices_p;

				for (int k = 0; k < 4; ++k) {
					if (alpha[i][k] > 0) {
						local_indices_p.push_back(k);
					}
				}

				for (int k = 0; k < 3; ++k) {
					VectorI vec_p {local_indices_p[k%3], local_indices_p[(k+1)%3], local_indices_p[(k+2)%3]};
					SIGMA_p.push_back(vec_p);
				}

				size_t SIGMA_p_size = SIGMA_p.size();

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (size_t k = 0; k < 4; ++k) {
						VectorI temp_alpha = alpha[j];
						temp_alpha[l_q] -= 1;
						temp_alpha[k] -= 1;

						if (temp_alpha[k] < 0 || temp_alpha[l_q] < 0) {
							continue;
						}

						double M;
						M_alpha_beta(M,
									 e[SIGMA_p[t][0]],
									 temp_alpha);

						double S1;
						S_ij_kl(S1,
								SIGMA_p[t][1],
								SIGMA_p[t][2],
								k,
								l_q,
								S_1,
								all_edges);


						x += M * S1;
					}
				}

				x = x * (n + 1) * (n + 2);
			}

			// 36d
			else if (i < ordered_basis_sizes[0] && 
					(j >= ordered_basis_sizes[4])) {

				x = 0;

				Vector2I SIGMA_p;
				VectorI local_indices_p;

				for (int k = 0; k < 4; ++k) {
					if (alpha[i][k] > 0) {
						local_indices_p.push_back(k);
					}
				}

				for (int k = 0; k < 3; ++k) {
					VectorI vec_p {local_indices_p[k%3], local_indices_p[(k+1)%3], local_indices_p[(k+2)%3]};
					SIGMA_p.push_back(vec_p);
				}

				size_t SIGMA_p_size = SIGMA_p.size();

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (size_t k = 0; k < 4; ++k) {
						Vector2I SIGMA_q;
						VectorI local_indices_q;
						if (k == 0) {
							local_indices_q = {0, 1, 2};
						}
						else if (k == 1) {
							local_indices_q = {0, 1, 3};
						}
						else if (k == 2) {
							local_indices_q = {0, 2, 3};
						}
						else if (k == 3) {
							local_indices_q = {1, 2, 3};
						}

						for (int v = 0; v < 3; ++v) {
							VectorI vec_q {local_indices_q[v%3], local_indices_q[(v+1)%3], local_indices_q[(v+2)%3]};
							SIGMA_q.push_back(vec_q);
						}

						size_t SIGMA_q_size = SIGMA_q.size();

						for (size_t s = 0; s < SIGMA_q_size; ++s) {
							VectorI temp_alpha = alpha[j];
							temp_alpha[k] += 1;

							double M;
							M_alpha_beta(M,
										 e[SIGMA_p[t][0]],
										 temp_alpha);

							double S1;
							S_ij_kl(S1,
									SIGMA_p[t][1],
									SIGMA_p[t][2],
									SIGMA_q[s][1],
									SIGMA_q[s][2],
									S_1,
									all_edges);

							x += pow(-1, k) * (alpha[j][k] + 1) * alpha[j][SIGMA_q[s][0]] * M * S1;
						}
					}
				}
			}

			//36e
			else if ((i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) &&
					 (j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1])) {

				x = 0;

				Vector2I SIGMA_p;
				Vector2I SIGMA_q;
				VectorI local_indices_p;
				VectorI local_indices_q;

				int E_nF_size = ordered_basis_sizes[1] - ordered_basis_sizes[0];
				int face_index_p = std::floor((i - ordered_basis_sizes[0])/(E_nF_size/4));
				VectorI face_p = all_faces[face_index_p];

				for (int k = 0; k < 4; ++k) {
					if (face_p[k] > 0) {
						local_indices_p.push_back(k);
					}
				}

				int face_index_q = std::floor((j - ordered_basis_sizes[0])/(E_nF_size/4));
				VectorI face_q = all_faces[face_index_q];

				for (int k = 0; k < 4; ++k) {
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

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (int l = 0; l < 3; ++l) {
						for (size_t s = 0; s < SIGMA_q_size; ++s) {
							for (int k = 0; k < 3; ++k) {
								VectorI temp_alpha1 = alpha[i];
								temp_alpha1[SIGMA_p[t][0]] += 1;
								temp_alpha1[SIGMA_p[t][l]] -= 1;

								if (temp_alpha1[SIGMA_p[t][l]] < 0) {
									continue;								
								}

								VectorI temp_alpha2 = alpha[j];
								temp_alpha2[SIGMA_q[s][0]] += 1;
								temp_alpha2[SIGMA_q[s][k]] -= 1;

								if (temp_alpha2[SIGMA_q[s][k]] < 0) {
									continue;								
								}

								double M;
								M_alpha_beta(M,
									     	 temp_alpha1,
									     	 temp_alpha2);

								double S1;
								S_ij_kl(S1,
										SIGMA_p[t][l],
										SIGMA_p[t][1],
										SIGMA_q[s][k],
										SIGMA_q[s][1],
										S_1,
										all_edges);
								double S2;
								S_ij_kl(S2,
										SIGMA_p[t][l],
										SIGMA_p[t][1],
										SIGMA_q[s][k],
										SIGMA_q[s][2],
										S_1,
										all_edges);
								double S3;
								S_ij_kl(S3,
										SIGMA_p[t][l],
										SIGMA_p[t][2],
										SIGMA_q[s][k],
										SIGMA_q[s][1],
										S_1,
										all_edges);
								double S4;
								S_ij_kl(S4,
										SIGMA_p[t][l],
										SIGMA_p[t][2],
										SIGMA_q[s][k],
										SIGMA_q[s][2],
										S_1,
										all_edges);

								x += M * (alpha[i][SIGMA_p[t][0]] + 1) * (alpha[j][SIGMA_q[s][0]] + 1)
									   * (  alpha[i][SIGMA_p[t][2]] * alpha[j][SIGMA_q[s][2]] * S1
									   	  - alpha[i][SIGMA_p[t][2]] * alpha[j][SIGMA_q[s][1]] * S2
									   	  - alpha[i][SIGMA_p[t][1]] * alpha[j][SIGMA_q[s][2]] * S3
									   	  + alpha[i][SIGMA_p[t][1]] * alpha[j][SIGMA_q[s][1]] * S4);
							}
						}
					}
				}
				x = x * pow(n+1, 2);
			}

			//36f
			else if ((i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) &&
					 (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[4])) {

				x = 0;

				int l_q = 0;
				if (j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) {
					l_q = 1;
				}
				else if (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]) {
					l_q = 2;
				}

				Vector2I SIGMA_p;
				VectorI local_indices_p;

				int E_nF_size = ordered_basis_sizes[1] - ordered_basis_sizes[0];
				int face_index_p = std::floor((i - ordered_basis_sizes[0])/(E_nF_size/4));
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

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (int l = 0; l < 3; ++l) {
						for (int k = 0; k < 4; ++k) {
							VectorI temp_alpha1 = alpha[i];
							temp_alpha1[SIGMA_p[t][0]] += 1;
							temp_alpha1[SIGMA_p[t][l]] -= 1;

							if (temp_alpha1[SIGMA_p[t][l]] < 0) {
								continue;
							}

							VectorI temp_alpha2 = alpha[j];
							temp_alpha2[l_q] -= 1;
							temp_alpha2[k] -= 1;

							if (temp_alpha2[k] < 0 || temp_alpha2[l_q] < 0) {
								continue;
							}

							double M;
							M_alpha_beta(M,
										 temp_alpha1,
										 temp_alpha2);

							double S1;
							S_ij_kl(S1,
									SIGMA_p[t][l],
									SIGMA_p[t][1],
									k,
									l_q,
									S_1,
									all_edges);
							double S2;
							S_ij_kl(S2,
									SIGMA_p[t][l],
									SIGMA_p[t][2],
									k,
									l_q,
									S_1,
									all_edges);

							x += M * (alpha[i][SIGMA_p[t][0]] + 1) 
							       * (alpha[i][SIGMA_p[t][2]] * S1 - alpha[i][SIGMA_p[t][1]] * S2);
						}
					}
				}
				x = x * (n + 2) * pow(n + 1, 2);
			}

			//36g
			else if ((i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) &&
					 (j >= ordered_basis_sizes[4])) {

				x = 0;

				Vector2I SIGMA_p;
				VectorI local_indices_p;

				int E_nF_size = ordered_basis_sizes[1] - ordered_basis_sizes[0];
				int face_index_p = std::floor((i - ordered_basis_sizes[0])/(E_nF_size/4));
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

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (int l = 0; l < 3; ++l) {
						for (int k = 0; k < 4; ++k) {
							Vector2I SIGMA_q;
							VectorI local_indices_q;
							if (k == 0) {
								local_indices_q = {0, 1, 2};
							}
							else if (k == 1) {
								local_indices_q = {0, 1, 3};
							}
							else if (k == 2) {
								local_indices_q = {0, 2, 3};
							}
							else if (k == 3) {
								local_indices_q = {1, 2, 3};
							}

							for (int v = 0; v < 3; ++v) {
								VectorI vec_q {local_indices_q[v%3], local_indices_q[(v+1)%3], local_indices_q[(v+2)%3]};
								SIGMA_q.push_back(vec_q);
							}

							size_t SIGMA_q_size = SIGMA_q.size();

							for (size_t s = 0; s < SIGMA_q_size; ++s) {

								VectorI temp_alpha1 = alpha[i];
								temp_alpha1[SIGMA_p[t][0]] += 1;
								temp_alpha1[SIGMA_p[t][l]] -= 1;

								if (temp_alpha1[SIGMA_p[t][l]] < 0) {
									continue;
								}

								VectorI temp_alpha2 = alpha[j];
								temp_alpha2[k] += 1;

								double M;
								M_alpha_beta(M,
											 temp_alpha1,
											 temp_alpha2);

								double S1;
								S_ij_kl(S1,
										SIGMA_p[t][l],
										SIGMA_p[t][1],
										SIGMA_q[s][1],
										SIGMA_q[s][2],
										S_1,
										all_edges);
								double S2;
								S_ij_kl(S2,
										SIGMA_p[t][l],
										SIGMA_p[t][2],
										SIGMA_q[s][1],
										SIGMA_q[s][2],
										S_1,
										all_edges);

								x += pow(-1, k) * M * (alpha[i][SIGMA_p[t][0]] + 1) * (alpha[j][k] + 1) * alpha[j][SIGMA_q[s][0]]
								       * (alpha[i][SIGMA_p[t][2]] * S1 - alpha[i][SIGMA_p[t][1]] * S2);
							}
						}
					}
				}
				x = x * (n + 1);
			}

			//36h
			else if ((i >= ordered_basis_sizes[1] && i < ordered_basis_sizes[4]) &&
					 (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[4])) {

				x = 0;

				int l_p = 0;
				if (i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) {
					l_p = 1;
				}
				else if (i >= ordered_basis_sizes[3] && i < ordered_basis_sizes[4]) {
					l_p = 2;
				}

				int l_q = 0;
				if (j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) {
					l_q = 1;
				}
				else if (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]) {
					l_q = 2;
				}

				for (int l = 0; l < 4; ++l) {
					for (int k = 0; k < 4; ++k) {
						VectorI temp_alpha1 = alpha[i];
						temp_alpha1[l_p] -= 1;
						temp_alpha1[l] -= 1;

						if (temp_alpha1[l_p] < 0 || temp_alpha1[l] < 0) {
							continue;
						}

						VectorI temp_alpha2 = alpha[j];
						temp_alpha2[l_q] -= 1;
						temp_alpha2[k] -= 1;

						if (temp_alpha2[l_q] < 0 || temp_alpha2[k] < 0) {
							continue;
						}

						double M;
						M_alpha_beta(M,
									 temp_alpha1,
									 temp_alpha2);

						double S;
						S_ij_kl(S,
								l,
								l_p,
								k,
								l_q,
								S_1,
								all_edges);

						x += M * S;
					}
				}
				x = x * pow(n+1, 2) * pow(n+2, 2);
			}

			//36i
			else if ((i >= ordered_basis_sizes[1] && i < ordered_basis_sizes[4]) && 
					 (j >= ordered_basis_sizes[4])) {

				x = 0;

				int l_p = 0;
				if (i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) {
					l_p = 1;
				}
				else if (i >= ordered_basis_sizes[3] && i < ordered_basis_sizes[4]) {
					l_p = 2;
				}

				for (size_t l = 0; l < 4; ++l) {
					for (size_t k = 0; k < 4; ++k) {
						Vector2I SIGMA_q;
						VectorI local_indices_q;
						if (k == 0) {
							local_indices_q = {0, 1, 2};
						}
						else if (k == 1) {
							local_indices_q = {0, 1, 3};
						}
						else if (k == 2) {
							local_indices_q = {0, 2, 3};
						}
						else if (k == 3) {
							local_indices_q = {1, 2, 3};
						}

						for (int v = 0; v < 3; ++v) {
							VectorI vec_q {local_indices_q[v%3], local_indices_q[(v+1)%3], local_indices_q[(v+2)%3]};
							SIGMA_q.push_back(vec_q);
						}

						size_t SIGMA_q_size = SIGMA_q.size();

						for (size_t s = 0; s < SIGMA_q_size; ++s) {

							VectorI temp_alpha1 = alpha[i];
							temp_alpha1[l_p] -= 1;
							temp_alpha1[l] -= 1;

							if (temp_alpha1[l_p] < 0 || temp_alpha1[l] < 0) {
								continue;
							}

							VectorI temp_alpha2 = alpha[j];
							temp_alpha2[k] += 1;

							double M;
							M_alpha_beta(M,
										 temp_alpha1,
										 temp_alpha2);

							double S1;
							S_ij_kl(S1,
									l,
									l_p,
									SIGMA_q[s][1],
									SIGMA_q[s][2],
									S_1,
									all_edges);

							x += pow(-1, k) * M * (alpha[j][k] + 1) * alpha[j][SIGMA_q[s][0]] * S1;
						}
					}
				}
				x = x * (n + 1) * (n + 2);
			}

			//36j
			else if (i >= ordered_basis_sizes[4] &&
					 j >= ordered_basis_sizes[4]) {

				x = 0;

				for (size_t l = 0; l < 4; ++l) {
					Vector2I SIGMA_p;
					VectorI local_indices_p;
					if (l == 0) {
						local_indices_p = {0, 1, 2};
					}
					else if (l == 1) {
						local_indices_p = {0, 1, 3};
					}
					else if (l == 2) {
						local_indices_p = {0, 2, 3};
					}
					else if (l == 3) {
						local_indices_p = {1, 2, 3};
					}

					for (int v = 0; v < 3; ++v) {
						VectorI vec_p {local_indices_p[v%3], local_indices_p[(v+1)%3], local_indices_p[(v+2)%3]};
						SIGMA_p.push_back(vec_p);
					}

					size_t SIGMA_p_size = SIGMA_p.size();

					for (size_t t = 0; t < SIGMA_p_size; ++t) {

						for (int k = 0; k < 4; ++k) {
							Vector2I SIGMA_q;
							VectorI local_indices_q;
							if (k == 0) {
								local_indices_q = {0, 1, 2};
							}
							else if (k == 1) {
								local_indices_q = {0, 1, 3};
							}
							else if (k == 2) {
								local_indices_q = {0, 2, 3};
							}
							else if (k == 3) {
								local_indices_q = {1, 2, 3};
							}

							for (int v = 0; v < 3; ++v) {
								VectorI vec_q {local_indices_q[v%3], local_indices_q[(v+1)%3], local_indices_q[(v+2)%3]};
								SIGMA_q.push_back(vec_q);
							}

							size_t SIGMA_q_size = SIGMA_q.size();

							for (size_t s = 0; s < SIGMA_q_size; ++s) {

								VectorI temp_alpha1 = alpha[i];
								temp_alpha1[l] += 1;

								VectorI temp_alpha2 = alpha[j];
								temp_alpha2[k] += 1;

								double M;
								M_alpha_beta(M,
											 temp_alpha1,
											 temp_alpha2);

								double S1;
								S_ij_kl(S1,
										SIGMA_p[t][1],
										SIGMA_p[t][2],
										SIGMA_q[s][1],
										SIGMA_q[s][2],
										S_1,
										all_edges);

								x += pow(-1, l+k) * M * (alpha[i][l] + 1) * alpha[i][SIGMA_p[t][0]] 
													  * (alpha[j][k] + 1) * alpha[j][SIGMA_q[s][0]]
													  * S1;
							}
						}
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

	for (int i = 0; i < size; ++i) {
		for (int j = i; j < size; ++j) {

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
					VectorI vec = alpha[j];
					vec[s] -= 1;

					bool flag = false;
					if (vec[s] < 0) {
						flag = true;
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

					VectorI vec = alpha[j];
					vec[SIGMA[k][0]] += 1;

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

				int l_q = 0;
				if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
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

					VectorI vec = alpha[j];
					vec[s] -= 1;

					bool flag = false;
					if (vec[s] < 0) {
						flag = true;
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

					VectorI vec1 = alpha[i];
					vec1[t] -= 1;

					bool flag1 = false;
					if (vec1[t] < 0) {
						flag1 = true;
					}

					if (flag1) {
						continue;
					}

					for (int s = 0; s < 4; ++s) {

						VectorI vec2 = alpha[j];
						vec2[s] -= 1;

						bool flag2 = false;
						if (vec2[s] < 0) {
							flag2 = true;
						}

						if (flag2) {
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

					VectorI vec1 = alpha[i];
					vec1[t] -= 1;

					bool flag1 = false;
					if (vec1[t] < 0) {
						flag1 = true;
					}

					if (flag1) {
						continue;
					}

					for (int k = 0; k < SIGMA_size; ++k) {

						VectorI vec2 = alpha[j];
						vec2[SIGMA[k][0]] += 1;

						double M;
						M_alpha_beta(M,
									 vec1,
									 vec2);

						int sigma_1 = SIGMA[k][0];
						int sigma_2 = SIGMA[k][1];
						int sigma_3 = SIGMA[k][2];

						x += M * (alpha[j][sigma_1] + 1) * (alpha[j][sigma_3] * S_0.coeffRef(t, sigma_2) - alpha[j][sigma_2] * S_0.coeffRef(t, sigma_3));
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

				int l_q = 0;
				if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
					l_q = 1;
				}
				else if (j >= ordered_basis_sizes[6] && j < ordered_basis_sizes[7]) {
					l_q = 2;
				}

				for (int t = 0; t < 4; ++t) {

					VectorI vec1 = alpha[i];
					vec1[t] -= 1;

					bool flag1 = false;
					if (vec1[t] < 0) {
						flag1 = true;
					}

					if (flag1) {
						continue;
					}

					for (int s = 0; s < 4; ++s) {
						int delta = 0;
						if (s == l_q) {
							delta = 1;
						}

						VectorI vec2 = alpha[j];
						vec2[s] -= 1;

						bool flag2 = false;
						if (vec2[s] < 0) {
							flag2 = true;
						}

						if (flag2) {
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

					VectorI vec1 = alpha[i];
					vec1[SIGMA_p[t][0]] += 1;

					for (int s = 0; s < SIGMA_q_size; ++s) {

						VectorI vec2 = alpha[j];
						vec2[SIGMA_q[s][0]] += 1;

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

				int l_q = 0;
				if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
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

					VectorI vec1 = alpha[i];
					vec1[SIGMA_p[t][0]] += 1;

					for (int s = 0; s < 4; ++s) {
						int delta = 0;
						if (s == l_q) {
							delta = 1;
						}

						VectorI vec2 = alpha[j];
						vec2[s] -= 1;

						bool flag2 = false;
						if (vec2[s] < 0) {
							flag2 = true;
						}

						if(flag2) {
							continue;
						}

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

				int l_q = 0;
				int l_p = 0;
				if (i >= ordered_basis_sizes[5] && i < ordered_basis_sizes[6]) {
					l_p = 1;
				}
				else if (i >= ordered_basis_sizes[6] && i < ordered_basis_sizes[7]) {
					l_p = 2;
				}
				if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
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

					VectorI vec1 = alpha[i];
					vec1[t] -= 1;

					bool flag1 = false;
					if (vec1[t] < 0) {
						flag1 = true;
					}

					if (flag1) {
						continue;
					}

					for (int s = 0; s < 4; ++s) {
						int delta_s = 0;
						if (l_q == s) {
							delta_s = 1;
						}

						VectorI vec2 = alpha[j];
						vec2[s] -= 1;

						bool flag2 = false;
						if (vec2[s] < 0) {
							flag2 = true;
						}

						if (flag2) {
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

			// 33k
			else if(((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) ||
				   (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) || 
				   (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4])) && 
				   (i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3])) {
				
				x = 0;
				
				int E_nF_size = ordered_basis_sizes[3] - ordered_basis_sizes[2];
				int face_index = std::floor((i - ordered_basis_sizes[2])/(E_nF_size/4));
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

					VectorI vec1 = alpha[j];
					vec1[t] -= 1;

					bool flag1 = false;
					if (vec1[t] < 0) {
						flag1 = true;
					}

					if (flag1) {
						continue;
					}

					for (int k = 0; k < SIGMA_size; ++k) {

						VectorI vec2 = alpha[i];
						vec2[SIGMA[k][0]] += 1;

						double M;
						M_alpha_beta(M,
									 vec1,
									 vec2);

						int sigma_1 = SIGMA[k][0];
						int sigma_2 = SIGMA[k][1];
						int sigma_3 = SIGMA[k][2];

						x += M * (alpha[i][sigma_1] + 1) * (alpha[i][sigma_3] * S_0.coeffRef(t, sigma_2) - alpha[i][sigma_2] * S_0.coeffRef(t, sigma_3));
					}
				}

				x = (n + 1) * x;
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
													  Vector2I &index_sets,
													  int d) {

	if (n < 0) {
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

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
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

			#ifdef MULTICORE
				#pragma omp critical
			#endif
			mass_matrix.coeffRef(i, j) = num/den;
		}
	}

	return SUCCESS;	
}


int FiniteElementExteriorCalculus::bb_stiffness_matrix_H_curl(DenMatD &stiffness_matrix,
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

	size_t size = alpha.size();
	stiffness_matrix.resize(size, size);

	DenMatD S_1;
	S_n(S_1,
		pts,
		1);

	Vector2I all_edges;
	compute_index_sets_o(all_edges,
						 2,
						 2);
	size_t all_edges_size = all_edges.size();

	Vector2I all_faces;
	compute_index_sets_o(all_faces,
						 3,
						 3);

	for (int i = 0; i < size; ++i) {
		for (int j = i; j < size; ++j) {

			double x = 0;
			
			// 34a
			if (i < ordered_basis_sizes[0] && j < ordered_basis_sizes[0]) {

				x = 0;

				VectorI local_indices_p;
				VectorI local_indices_q;

				for (int k = 0; k < 4; ++k) {
					if (alpha[i][k] > 0) {
						local_indices_p.push_back(k);
					}
					if (alpha[j][k] > 0) {
						local_indices_q.push_back(k);
					}
				}

				double S;
				S_ij_kl(S,
						local_indices_p[0],
						local_indices_p[1],
						local_indices_q[0],
						local_indices_q[1],
						S_1,
						all_edges);

				x = 4 * S;
			}

			//34c
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

				for (int s = 0; s < SIGMA_size; ++s) {
					double S1;
					S_ij_kl(S1,
							index_p[0],
							index_p[1],
							SIGMA[s][0],
							SIGMA[s][1],
							S_1,
							all_edges);
					double S2;
					S_ij_kl(S2,
							index_p[0],
							index_p[1],
							SIGMA[s][2],
							SIGMA[s][1],
							S_1,
							all_edges);
					double S3;
					S_ij_kl(S3,
							index_p[0],
							index_p[1],
							SIGMA[s][0],
							SIGMA[s][2],
							S_1,
							all_edges);
					double S4;
					S_ij_kl(S4,
							index_p[0],
							index_p[1],
							SIGMA[s][1],
							SIGMA[s][2],
							S_1,
							all_edges);

					x += (alpha[j][SIGMA[s][0]] + 1) * 
						 (alpha[j][SIGMA[s][2]] * (S1 + S2) - alpha[j][SIGMA[s][1]] * (S3 + S4));
				}
				x = x * 12/((n+2)*(n+3));
			}

			//34e
			else if ((i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) &&
					 (j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3])) {

				x = 0;

				Vector2I SIGMA_p;
				Vector2I SIGMA_q;
				VectorI local_indices_p;
				VectorI local_indices_q;

				int E_nF_size = ordered_basis_sizes[3] - ordered_basis_sizes[2];
				int face_index_p = std::floor((i - ordered_basis_sizes[2])/(E_nF_size/4));
				VectorI face_p = all_faces[face_index_p];

				for (int k = 0; k < 4; ++k) {
					if (face_p[k] > 0) {
						local_indices_p.push_back(k);
					}
				}

				int face_index_q = std::floor((j - ordered_basis_sizes[2])/(E_nF_size/4));
				VectorI face_q = all_faces[face_index_q];

				for (int k = 0; k < 4; ++k) {
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

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (int l = 0; l < 3; ++l) {
						for (size_t s = 0; s < SIGMA_q_size; ++s) {
							for (int k = 0; k < 3; ++k) {
								VectorI temp_alpha1 = alpha[i];
								temp_alpha1[SIGMA_p[t][0]] += 1;
								temp_alpha1[SIGMA_p[t][l]] -= 1;

								if (temp_alpha1[SIGMA_p[t][l]] < 0) {
									continue;								
								}

								VectorI temp_alpha2 = alpha[j];
								temp_alpha2[SIGMA_q[s][0]] += 1;
								temp_alpha2[SIGMA_q[s][k]] -= 1;

								if (temp_alpha2[SIGMA_q[s][k]] < 0) {
									continue;								
								}

								double M;
								M_alpha_beta(M,
									     	 temp_alpha1,
									     	 temp_alpha2);

								double S1;
								S_ij_kl(S1,
										SIGMA_p[t][l],
										SIGMA_p[t][1],
										SIGMA_q[s][k],
										SIGMA_q[s][1],
										S_1,
										all_edges);
								double S2;
								S_ij_kl(S2,
										SIGMA_p[t][l],
										SIGMA_p[t][1],
										SIGMA_q[s][k],
										SIGMA_q[s][2],
										S_1,
										all_edges);
								double S3;
								S_ij_kl(S3,
										SIGMA_p[t][l],
										SIGMA_p[t][2],
										SIGMA_q[s][k],
										SIGMA_q[s][1],
										S_1,
										all_edges);
								double S4;
								S_ij_kl(S4,
										SIGMA_p[t][l],
										SIGMA_p[t][2],
										SIGMA_q[s][k],
										SIGMA_q[s][2],
										S_1,
										all_edges);

								x += M * (alpha[i][SIGMA_p[t][0]] + 1) * (alpha[j][SIGMA_q[s][0]] + 1)
									   * (  alpha[i][SIGMA_p[t][2]] * alpha[j][SIGMA_q[s][2]] * S1
									   	  - alpha[i][SIGMA_p[t][2]] * alpha[j][SIGMA_q[s][1]] * S2
									   	  - alpha[i][SIGMA_p[t][1]] * alpha[j][SIGMA_q[s][2]] * S3
									   	  + alpha[i][SIGMA_p[t][1]] * alpha[j][SIGMA_q[s][1]] * S4);
							}
						}
					}
				}
				x = x * pow(n+1, 2);
			}

			//34f
			else if ((i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) &&
					 (j >= ordered_basis_sizes[4])) {

				x = 0;

				int l_q = 0;
				if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
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

				for (size_t t = 0; t < SIGMA_p_size; ++t) {
					for (int l = 0; l < 3; ++l) {
						for (int k = 0; k < 4; ++k) {
							VectorI temp_alpha1 = alpha[i];
							temp_alpha1[SIGMA_p[t][0]] += 1;
							temp_alpha1[SIGMA_p[t][l]] -= 1;

							if (temp_alpha1[SIGMA_p[t][l]] < 0) {
								continue;
							}

							VectorI temp_alpha2 = alpha[j];
							temp_alpha2[l_q] -= 1;
							temp_alpha2[k] -= 1;

							if (temp_alpha2[k] < 0 || temp_alpha2[l_q] < 0) {
								continue;
							}

							double M;
							M_alpha_beta(M,
										 temp_alpha1,
										 temp_alpha2);

							double S1;
							S_ij_kl(S1,
									SIGMA_p[t][l],
									SIGMA_p[t][1],
									k,
									l_q,
									S_1,
									all_edges);
							double S2;
							S_ij_kl(S2,
									SIGMA_p[t][l],
									SIGMA_p[t][2],
									k,
									l_q,
									S_1,
									all_edges);

							x += M * (alpha[i][SIGMA_p[t][0]] + 1) 
							       * (alpha[i][SIGMA_p[t][2]] * S1 - alpha[i][SIGMA_p[t][1]] * S2);
						}
					}
				}
				x = x * (n + 2) * pow(n + 1, 2);
			}

			//34g
			else if (i >= ordered_basis_sizes[4] && j >= ordered_basis_sizes[4]) {

				x = 0;

				int l_p = 0;
				if (i >= ordered_basis_sizes[5] && i < ordered_basis_sizes[6]) {
					l_p = 1;
				}
				else if (i >= ordered_basis_sizes[6] && i < ordered_basis_sizes[7]) {
					l_p = 2;
				}

				int l_q = 0;
				if (j >= ordered_basis_sizes[5] && j < ordered_basis_sizes[6]) {
					l_q = 1;
				}
				else if (j >= ordered_basis_sizes[6] && j < ordered_basis_sizes[7]) {
					l_q = 2;
				}

				for (int l = 0; l < 4; ++l) {
					for (int k = 0; k < 4; ++k) {
						VectorI temp_alpha1 = alpha[i];
						temp_alpha1[l_p] -= 1;
						temp_alpha1[l] -= 1;

						if (temp_alpha1[l_p] < 0 || temp_alpha1[l] < 0) {
							continue;
						}

						VectorI temp_alpha2 = alpha[j];
						temp_alpha2[l_q] -= 1;
						temp_alpha2[k] -= 1;

						if (temp_alpha2[l_q] < 0 || temp_alpha2[k] < 0) {
							continue;
						}

						double M;
						M_alpha_beta(M,
									 temp_alpha1,
									 temp_alpha2);

						double S;
						S_ij_kl(S,
								l,
								l_p,
								k,
								l_q,
								S_1,
								all_edges);

						x += M * S;
					}
				}
				x = x * pow(n+1, 2) * pow(n+2, 2);
			}

			stiffness_matrix.coeffRef(i, j) = x;
			if (i != j) {
				stiffness_matrix.coeffRef(j, i) = x;
			}
		}
	}

	return SUCCESS;
}


int FiniteElementExteriorCalculus::bb_stiffness_matrix_H_div(DenMatD &stiffness_matrix,
														    Vector2D &pts,
														    int n,
														    Vector2I &alpha,
														    VectorI &ordered_basis_sizes,
														    DenMatD &grad_bary_coords) {

	if (n < 0) {
		return FAILURE;
	}

	Vector2I e;
	compute_index_sets_o(e,
						 1,
						 1);

	size_t size = alpha.size();
	stiffness_matrix.resize(size, size);

	for (int i = 0; i < size; ++i) {
		for (int j = i; j < size; ++j) {

			double x = 0;

			if (i < ordered_basis_sizes[0] && j < ordered_basis_sizes[0]) {
				double e;
				e_ijk(e,
					  1,
					  2,
					  3,
					  grad_bary_coords);

				x = pow(-1, i + j) * 9 * pow(e, 2) * get_simplex_volume(pts); 
			}

			else if(i >= ordered_basis_sizes[4] && j >= ordered_basis_sizes[4]) {

				x = 0;

				for (int u = 0; u < 4; ++u) {
					for (int l = 0; l < 4; ++l) {
						for (int v = 0; v < 4; ++v) {
							for (int k = 0; k < 4; ++k) {
								
								int delta_ul = 0;
								if (u == l) {
									delta_ul = 1;
								}
								int delta_vk = 0;
								if (v == k) {
									delta_vk = 1;
								}

								VectorI temp_alpha1 = alpha[i];
								temp_alpha1[u] += 1;
								temp_alpha1[l] -= 1;
								if (temp_alpha1[l] < 0) {
									continue;
								}

								VectorI temp_alpha2 = alpha[j];
								temp_alpha2[v] += 1;
								temp_alpha2[k] -= 1;
								if (temp_alpha2[k] < 0) {
									continue;
								}

								double M;
								M_alpha_beta(M,
											 temp_alpha1,
											 temp_alpha2);
								double e;
								e_ijk(e,
									  1,
									  2,
									  3,
									  grad_bary_coords);

								x += pow(-1, u + v) * M
									  * (alpha[i][u] + 1) * (delta_ul * n - alpha[i][l])
									  * (alpha[j][v] + 1) * (delta_vk * n - alpha[j][k])
									  * pow(e, 2) * get_simplex_volume(pts);
							}
						}
					}
				}
				x = x * pow(n + 1, 2);
			}

			stiffness_matrix.coeffRef(i, j) = x;
			if (i != j) {
				stiffness_matrix.coeffRef(j, i) = x;
			}
		}
	}

	return SUCCESS;
}

double FiniteElementExteriorCalculus::bb_error_H_curl_1d_quad(int n,
														      Vector3I &simplices,
															  Vector2D &vertices,
															  VectorI &num_simplices,
															  int q_order) {

	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif

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
												   int q_order) {

	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif

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

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {

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
												      int q_order) {

	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif

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
			temp_alpha.clear();
			compute_index_sets_o(temp_alpha,
								 n + 2,
								 i);
			
			size_t temp_alpha_size = temp_alpha.size();

			for (int l = 0; l < 2; ++l) {	
				if (temp_alpha_size != 0) {
					alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
				}
				ordered_basis_sizes.push_back(total + temp_alpha_size);
				total += temp_alpha_size;
			}

			size_t counter = 0; 
			for (size_t j = 0; j < temp_alpha_size; ++j) {
				if (temp_alpha[j][2] == 1) {
					alpha.push_back(temp_alpha[j]);
					++counter;
				}
			}

			ordered_basis_sizes.push_back(total + counter);	
			total += counter;
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
	for(size_t s = 0; s < num_simplices[N-1]; ++s) {
		double e = 0;

		Vector2D pts;
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplices[N-1][s][k]]);
		}
		double vol = get_simplex_volume(pts);

		DenMatD grad_bary_coords;
		barycentric_gradients(grad_bary_coords,
							  pts);

		VectorDenMatD basis_elements;
		for(int i = 0; i < (int)alpha_size; ++i) {
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
				inner_product += weights[k] * f.dot(basis_elements[j].row(k));
			}

			b.coeffRef(j) = vol * inner_product/sum_weights;
		}

		EigVectorD coeffs = M.colPivHouseholderQr().solve(b);

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {

			EigVectorD f_dash = EigVectorD::Zero(embed_dim);
			for (size_t j = 0; j < alpha_size; ++j) {
				f_dash += coeffs.coeffRef(j) * basis_elements[j].row(node_index);
			}

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


double FiniteElementExteriorCalculus::bb_error_H_div(int n,
												     int q_order) {

	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif

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
						 3,
						 3);
	ordered_basis_sizes.push_back(alpha.size());

	size_t total = alpha.size();

	temp_alpha.clear();
	compute_index_sets_p(temp_alpha,
						 n,
						 3);
	size_t temp_alpha_size = temp_alpha.size();
	if (temp_alpha_size != 0) {
		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	}
	ordered_basis_sizes.push_back(total + temp_alpha_size);
	total += temp_alpha_size;

	temp_alpha.clear();
	compute_index_sets_o(temp_alpha,
						 n + 2,
						 4);
	temp_alpha_size = temp_alpha.size();
	for (int l = 0; l < 2; ++l) {	
		if (temp_alpha_size != 0) {
			alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
		}
		ordered_basis_sizes.push_back(total + temp_alpha_size);
		total += temp_alpha_size;
	}
	int counter = 0; 
	for (size_t j = 0; j < temp_alpha_size; ++j) {
		if (temp_alpha[j][2] == 1) {
			alpha.push_back(temp_alpha[j]);
			++counter;
		}
	}
	ordered_basis_sizes.push_back(total + counter);	
	total += counter;

	temp_alpha.clear();
	compute_index_sets_p(temp_alpha,
						 n,
						 4);
	temp_alpha_size = temp_alpha.size();
	if (temp_alpha_size != 0) {
		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	}
	ordered_basis_sizes.push_back(total + temp_alpha_size);
	total += temp_alpha_size;

	size_t alpha_size = alpha.size();

	Vector2I all_faces;
	compute_index_sets_o(all_faces,
						 3,
						 3);

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t s = 0; s < num_simplices[N-1]; ++s) {
		double e = 0;

		Vector2D pts;
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplices[N-1][s][k]]);
		}
		double vol = get_simplex_volume(pts);

		DenMatD grad_bary_coords;
		barycentric_gradients(grad_bary_coords,
							  pts);

		VectorDenMatD basis_elements;
		for(int i = 0; i < (int)alpha_size; ++i) {
			DenMatD temp_basis_elements(nodes_size, embed_dim);
			
			for(size_t j = 0; j < nodes_size; ++j) {
				if (i < ordered_basis_sizes[0]) {
					EigVectorD chi;

					VectorI local_indices;
					for (int j = 0; j < 4; ++j) {
						if (alpha[i][j] > 0) {
							local_indices.push_back(j);
						}
					}
					
					chi_l(chi,
						  alpha[i],
						  nodes[j],
						  grad_bary_coords,
						  local_indices);
					
					temp_basis_elements.row(j) = chi;
				}
				else if (i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) {
					EigVectorD curl_phi;

					int E_nF_size = ordered_basis_sizes[1] - ordered_basis_sizes[0];
					int face_index = std::floor((i - ordered_basis_sizes[0])/(E_nF_size/4));
					VectorI face = all_faces[face_index];
					
					VectorI temp_local_indices;

					for (int k = 0; k < 4; ++k) {
						if (face[k] > 0) {
							temp_local_indices.push_back(k);
						}
					}

					curl_phi_FT(curl_phi,
							   alpha[i],
							   n,
							   nodes[j],
							   grad_bary_coords,
							   temp_local_indices);

					temp_basis_elements.row(j) = curl_phi;
				}
				else if (i >= ordered_basis_sizes[1] && i < ordered_basis_sizes[2]) {
					EigVectorD curl_psi;

					curl_psi_T(curl_psi,
							   alpha[i],
							   n + 1,
							   0,
							   nodes[j],
							   grad_bary_coords);

					temp_basis_elements.row(j) = curl_psi;
				}
				else if (i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) {
					EigVectorD curl_psi;

					curl_psi_T(curl_psi,
							   alpha[i],
							   n + 1,
							   1,
							   nodes[j],
							   grad_bary_coords);

					temp_basis_elements.row(j) = curl_psi;
				}
				else if (i >= ordered_basis_sizes[3] && i < ordered_basis_sizes[4]) {
					EigVectorD curl_psi;

					curl_psi_T(curl_psi,
							   alpha[i],
							   n + 1,
							   2,
							   nodes[j],
							   grad_bary_coords);

					temp_basis_elements.row(j) = curl_psi;
				}
				else if (i >= ordered_basis_sizes[4] && i < ordered_basis_sizes[5]) {
					EigVectorD upsi;

					upsilon(upsi,
						    alpha[i],
						    n,
						    nodes[j],
						    grad_bary_coords);

					temp_basis_elements.row(j) = upsi;
				}
			}

			basis_elements.push_back(temp_basis_elements);
		}

		DenMatD M;
		bb_mass_matrix_H_div(M,
							 pts,
						 	 n,
						 	 alpha,
						 	 ordered_basis_sizes);

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
				inner_product += weights[k] * f.dot(basis_elements[j].row(k));
			}

			b.coeffRef(j) = vol * inner_product/sum_weights;
		}

		std::cout<<b<<"\n\n";
		EigVectorD coeffs = M.colPivHouseholderQr().solve(b);

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {
		
			EigVectorD f_dash = EigVectorD::Zero(embed_dim);
			for (size_t j = 0; j < alpha_size; ++j) {
				f_dash += coeffs.coeffRef(j) * basis_elements[j].row(node_index);
			}

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


double FiniteElementExteriorCalculus::bb_error_stiffness_H_curl(int n,
												      			int q_order) {

	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif

	return SUCCESS;
}


int FiniteElementExteriorCalculus::compute_bb_mass_matrices(int x,
															int n) {

	int num_edges;
	binomialCoeff(num_edges,
				  complex_dimension + 1,
				  2);
	int num_faces;
	binomialCoeff(num_faces,
				  complex_dimension + 1,
				  3);

	if (x == 0) {
		VectorTripletD triplet;

		Vector2I alpha;
		Vector2I temp_alpha;
		VectorI ordered_basis_sizes;
		VectorI sizes;

		compute_index_sets_o(alpha,
						 	 1,
						 	 1);

		int d = alpha[0].size();
		ordered_basis_sizes.push_back(alpha.size());

		for(int i = 2; i <= d; ++i) {
			temp_alpha.clear();
			compute_index_sets_o(temp_alpha,
							 	 n,
							 	 i);

			alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
			ordered_basis_sizes.push_back(alpha.size());
		}
		size_t alpha_size = alpha.size();

		VectorI ndofs;
		size_t size = 0;
		for (size_t i = 0; i < complex_dimension + 1; ++i) {
			int temp;
			binomialCoeff(temp,
						  std::max(0, n-1),
						  i);
			ndofs.push_back(temp);
			sizes.push_back(size + temp * num_simplices[i]);
			size += temp * num_simplices[i];
		}

		bb_mass_matrices[0].resize(size, size);

		DenMatD temp_mass_matrix;
		bb_mass_matrix_H_1(temp_mass_matrix,
					   	   n,
					   	   alpha);

		for (size_t i = 0; i < num_simplices[complex_dimension]; ++i) {
			Vector2D pts;
			for(size_t k = 0; k < complex_dimension + 1; ++k) {
				pts.push_back(vertices[simplices[complex_dimension][i][k]]);
			}
			double vol = get_simplex_volume(pts);

			DenMatD mass_matrix = temp_mass_matrix * vol;

			for (int j = 0; j < (int)alpha_size; ++j) {
				for (int k = j; k < (int)alpha_size; ++k) {
					if (k < ordered_basis_sizes[0] && j < ordered_basis_sizes[0]) {
						triplet.push_back(TripletD(simplex_sub_simplices[i][0][j], simplex_sub_simplices[i][0][k], mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(simplex_sub_simplices[i][0][k], simplex_sub_simplices[i][0][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j < ordered_basis_sizes[0]) && 
							(k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1])) {
						int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
						int temp_index = std::floor((k - ordered_basis_sizes[0])/(temp/num_edges));
						int index = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index] + (k - ordered_basis_sizes[0])%ndofs[1];
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][0][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][0][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j < ordered_basis_sizes[0]) && 
							(k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2])) {
						int temp = ordered_basis_sizes[2] - ordered_basis_sizes[1];
						int temp_index = std::floor((k - ordered_basis_sizes[1])/(temp/num_faces));
						int index = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index] + (k - ordered_basis_sizes[1])%ndofs[2];
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][0][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][0][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j < ordered_basis_sizes[0]) && 
							(k >= ordered_basis_sizes[2] && k < ordered_basis_sizes[3])) {
						int temp = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp_index = std::floor((k - ordered_basis_sizes[2])/(temp));
						int index = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index] + (k - ordered_basis_sizes[2])%ndofs[3];
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][0][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][0][j], mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) && 
							(k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1])) {
						int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
						int temp_index1 = std::floor((j - ordered_basis_sizes[0])/(temp/num_edges));
						int temp_index2 = std::floor((k - ordered_basis_sizes[0])/(temp/num_edges));
						int index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index1] + (j - ordered_basis_sizes[0])%ndofs[1];
						int index2 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index2] + (k - ordered_basis_sizes[0])%ndofs[1];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) && 
							(k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2])) {
						int temp1 = ordered_basis_sizes[1] - ordered_basis_sizes[0];
						int temp2 = ordered_basis_sizes[2] - ordered_basis_sizes[1];
						int temp_index1 = std::floor((j - ordered_basis_sizes[0])/(temp1/num_edges));
						int temp_index2 = std::floor((k - ordered_basis_sizes[1])/(temp2/num_faces));
						int index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index1] + (j - ordered_basis_sizes[0])%ndofs[1];
						int index2 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index2] + (k - ordered_basis_sizes[1])%ndofs[2];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) && 
							(k >= ordered_basis_sizes[2] && k < ordered_basis_sizes[3])) {
						int temp1 = ordered_basis_sizes[1] - ordered_basis_sizes[0];
						int temp2 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp_index1 = std::floor((j - ordered_basis_sizes[0])/(temp1/num_edges));
						int temp_index2 = std::floor((k - ordered_basis_sizes[2])/(temp2));
						int index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index1] + (j - ordered_basis_sizes[0])%ndofs[1];
						int index2 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[2])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) && 
							(k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2])) {
						int temp1 = ordered_basis_sizes[2] - ordered_basis_sizes[1];
						int temp2 = ordered_basis_sizes[2] - ordered_basis_sizes[1];
						int temp_index1 = std::floor((j - ordered_basis_sizes[1])/(temp1/num_faces));
						int temp_index2 = std::floor((k - ordered_basis_sizes[1])/(temp2/num_faces));
						int index1 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[1])%ndofs[2];
						int index2 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index2] + (k - ordered_basis_sizes[1])%ndofs[2];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) && 
							(k >= ordered_basis_sizes[2] && k < ordered_basis_sizes[3])) {
						int temp1 = ordered_basis_sizes[2] - ordered_basis_sizes[1];
						int temp2 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp_index1 = std::floor((j - ordered_basis_sizes[1])/(temp1/num_faces));
						int temp_index2 = std::floor((k - ordered_basis_sizes[2])/(temp2));
						int index1 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[1])%ndofs[2];
						int index2 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[2])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) && 
							(k >= ordered_basis_sizes[2] && k < ordered_basis_sizes[3])) {
						int temp1 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp2 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp_index1 = std::floor((j - ordered_basis_sizes[2])/(temp1));
						int temp_index2 = std::floor((k - ordered_basis_sizes[2])/(temp2));
						int index1 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[2])%ndofs[3];
						int index2 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[2])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
				}
			}
		}

		bb_mass_matrices[0].setFromTriplets(triplet.begin(), triplet.end());
		bb_mass_matrices[0].makeCompressed();
	}

	else if (x == 1) {
		VectorTripletD triplet;

		Vector2I alpha;
		Vector2I temp_alpha;
		VectorI ordered_basis_sizes;
		VectorI sizes;

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
				temp_alpha.clear();
				compute_index_sets_o(temp_alpha,
									 n + 2,
									 i);
				
				size_t temp_alpha_size = temp_alpha.size();

				for (int l = 0; l < 2; ++l) {	
					if (temp_alpha_size != 0) {
						alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
					}
					ordered_basis_sizes.push_back(total + temp_alpha_size);
					total += temp_alpha_size;
				}

				size_t counter = 0; 
				for (size_t j = 0; j < temp_alpha_size; ++j) {
					if (temp_alpha[j][2] == 1) {
						alpha.push_back(temp_alpha[j]);
						++counter;
					}
				}

				ordered_basis_sizes.push_back(total + counter);	
				total += counter;
			}
		}

		size_t alpha_size = alpha.size();

		VectorI ndofs;
		size_t size = 0;
		for (size_t i = 0; i < complex_dimension + 1; ++i) {
			int temp;
			binomialCoeff(temp,
						  n,
						  i);
			ndofs.push_back(temp);
			sizes.push_back(size + temp * num_simplices[std::max(1, (int)i)]);
			size += temp * num_simplices[std::max(1, (int)i)];

			if (i == 2) {
				temp = 0;
				binomialCoeff(temp,
							  n+2,
							  2);
				temp = temp - 1;
				ndofs.push_back(temp);
				sizes.push_back(size + temp * num_simplices[std::max(1, (int)i)]);
				size += temp * num_simplices[std::max(1, (int)i)];
			}

			else if (i == 3) {
				int temp1;
				int temp2;
				binomialCoeff(temp1,
							  n+1,
							  3);
				binomialCoeff(temp2,
							  n,
							  2);
				temp = 2*temp1 + temp2;
				ndofs.push_back(temp);
				sizes.push_back(size + temp * num_simplices[std::max(1, (int)i)]);
				size += temp * num_simplices[std::max(1, (int)i)];
			}
		}

		bb_mass_matrices[1].resize(size, size);

		for (size_t i = 0; i < num_simplices[complex_dimension]; ++i) {
			Vector2D pts;
			VectorI simplex = simplices[complex_dimension][i];
			std::sort(simplex.begin(), simplex.end());
			for(size_t k = 0; k < complex_dimension + 1; ++k) {
				pts.push_back(vertices[simplex[k]]);
			}

			DenMatD mass_matrix;
			bb_mass_matrix_H_curl(mass_matrix,
								  pts,
							   	  n,
							   	  alpha,
							   	  ordered_basis_sizes);

			for (int j = 0; j < (int)alpha_size; ++j) {
				for (int k = j; k < (int)alpha_size; ++k) {
					if (k < ordered_basis_sizes[0] && j < ordered_basis_sizes[0]) {
						triplet.push_back(TripletD(simplex_sub_simplices[i][1][j], simplex_sub_simplices[i][1][k], mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(simplex_sub_simplices[i][1][k], simplex_sub_simplices[i][1][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if (j < ordered_basis_sizes[0] && 
							((k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1]) ||
							 (k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2]) ||
							 (k >= ordered_basis_sizes[3] && k < ordered_basis_sizes[4]))) {
						int index = -1;
						if (k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1]) {
							int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
							int temp_index = std::floor((k - ordered_basis_sizes[0])/(temp/num_edges));
							index = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index] + (k - ordered_basis_sizes[0])%ndofs[1];
						}
						else if (k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2]) {
							int temp = ordered_basis_sizes[2] - ordered_basis_sizes[1];
							int temp_index = std::floor((k - ordered_basis_sizes[1])/(temp/num_faces));
							index = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index] + (k - ordered_basis_sizes[1])%ndofs[2];
						}
						else if (k >= ordered_basis_sizes[3] && k < ordered_basis_sizes[4]) {
							int temp = ordered_basis_sizes[4] - ordered_basis_sizes[3];
							int temp_index = std::floor((k - ordered_basis_sizes[3])/(temp));
							index = sizes[3] + ndofs[4] * simplex_sub_simplices[i][3][temp_index] + (k - ordered_basis_sizes[3])%ndofs[4];
						}	
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][1][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][1][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j < ordered_basis_sizes[0]) && 
							(k >= ordered_basis_sizes[2] && k < ordered_basis_sizes[3])) {
						int temp = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp_index = std::floor((k - ordered_basis_sizes[2])/(temp/num_faces));
						int index = sizes[2] + ndofs[3] * simplex_sub_simplices[i][2][temp_index] + (k - ordered_basis_sizes[2])%ndofs[3];
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][1][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][1][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if (j < ordered_basis_sizes[0] && 
							 k >= ordered_basis_sizes[4]) {
						int temp = ordered_basis_sizes[7] - ordered_basis_sizes[4];
						int temp_index = std::floor((k - ordered_basis_sizes[4])/(temp));
						int index = sizes[4] + ndofs[5] * simplex_sub_simplices[i][3][temp_index] + (k - ordered_basis_sizes[4])%ndofs[5];
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][1][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][1][j], mass_matrix.coeffRef(k, j)));
						}
					}

					else if (((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) ||
							  (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) ||
							  (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4])) &&
							 ((k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1]) ||
							  (k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2]) ||
							  (k >= ordered_basis_sizes[3] && k < ordered_basis_sizes[4]))) {
						int index1 = -1;
						int index2 = -1;

						if (j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) {
							int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
							int temp_index = std::floor((j - ordered_basis_sizes[0])/(temp/num_edges));
							index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index] + (j - ordered_basis_sizes[0])%ndofs[1];
						}
						else if (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) {
							int temp = ordered_basis_sizes[2] - ordered_basis_sizes[1];
							int temp_index = std::floor((j - ordered_basis_sizes[1])/(temp/num_faces));
							index1 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index] + (j - ordered_basis_sizes[1])%ndofs[2];
						}
						else if (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]) {
							int temp = ordered_basis_sizes[4] - ordered_basis_sizes[3];
							int temp_index = std::floor((j - ordered_basis_sizes[3])/(temp));
							index1 = sizes[3] + ndofs[4] * simplex_sub_simplices[i][3][temp_index] + (j - ordered_basis_sizes[3])%ndofs[4];
						}

						if (k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1]) {
							int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
							int temp_index = std::floor((k - ordered_basis_sizes[0])/(temp/num_edges));
							index2 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index] + (k - ordered_basis_sizes[0])%ndofs[1];
						}
						else if (k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2]) {
							int temp = ordered_basis_sizes[2] - ordered_basis_sizes[1];
							int temp_index = std::floor((k - ordered_basis_sizes[1])/(temp/num_faces));
							index2 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index] + (k - ordered_basis_sizes[1])%ndofs[2];
						}
						else if (k >= ordered_basis_sizes[3] && k < ordered_basis_sizes[4]) {
							int temp = ordered_basis_sizes[4] - ordered_basis_sizes[3];
							int temp_index = std::floor((k - ordered_basis_sizes[3])/(temp));
							index2 = sizes[3] + ndofs[4] * simplex_sub_simplices[i][3][temp_index] + (k - ordered_basis_sizes[3])%ndofs[4];
						}
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if (((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) ||
							  (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) ||
							  (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4])) &&
							  (k >= ordered_basis_sizes[2] && k < ordered_basis_sizes[3])) {
						int index1 = -1;

						if (j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) {
							int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
							int temp_index = std::floor((j - ordered_basis_sizes[0])/(temp/num_edges));
							index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index] + (j - ordered_basis_sizes[0])%ndofs[1];
						}
						else if (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) {
							int temp = ordered_basis_sizes[2] - ordered_basis_sizes[1];
							int temp_index = std::floor((j - ordered_basis_sizes[1])/(temp/num_faces));
							index1 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index] + (j - ordered_basis_sizes[1])%ndofs[2];
						}
						else if (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]) {
							int temp = ordered_basis_sizes[4] - ordered_basis_sizes[3];
							int temp_index = std::floor((j - ordered_basis_sizes[3])/(temp));
							index1 = sizes[3] + ndofs[4] * simplex_sub_simplices[i][3][temp_index] + (j - ordered_basis_sizes[3])%ndofs[4];
						}

						int temp2 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp_index2 = std::floor((k - ordered_basis_sizes[2])/(temp2/num_faces));
						int index2 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][2][temp_index2] + (k - ordered_basis_sizes[2])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if (((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) ||
							  (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) ||
							  (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4])) &&
							  (k >= ordered_basis_sizes[4])) {

						int index1 = -1;

						if (j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) {
							int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
							int temp_index = std::floor((j - ordered_basis_sizes[0])/(temp/num_edges));
							index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index] + (j - ordered_basis_sizes[0])%ndofs[1];
						}
						else if (j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[2]) {
							int temp = ordered_basis_sizes[2] - ordered_basis_sizes[1];
							int temp_index = std::floor((j - ordered_basis_sizes[1])/(temp/num_faces));
							index1 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index] + (j - ordered_basis_sizes[1])%ndofs[2];
						}
						else if (j >= ordered_basis_sizes[3] && j < ordered_basis_sizes[4]) {
							int temp = ordered_basis_sizes[4] - ordered_basis_sizes[3];
							int temp_index = std::floor((j - ordered_basis_sizes[3])/(temp));
							index1 = sizes[3] + ndofs[4] * simplex_sub_simplices[i][3][temp_index] + (j - ordered_basis_sizes[3])%ndofs[4];
						}

						int temp2 = ordered_basis_sizes[7] - ordered_basis_sizes[4];
						int temp_index2 = std::floor((k - ordered_basis_sizes[4])/(temp2));
						int index2 = sizes[4] + ndofs[5] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[4])%ndofs[5];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) && 
							(k >= ordered_basis_sizes[2] && k < ordered_basis_sizes[3])) {
						int temp1 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp2 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp_index1 = std::floor((j - ordered_basis_sizes[2])/(temp1/num_faces));
						int temp_index2 = std::floor((k - ordered_basis_sizes[2])/(temp2/num_faces));
						int index1 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[2])%ndofs[3];
						int index2 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][2][temp_index2] + (k - ordered_basis_sizes[2])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) && 
							(k >= ordered_basis_sizes[4])) {
						int temp1 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp2 = ordered_basis_sizes[7] - ordered_basis_sizes[4];
						int temp_index1 = std::floor((j - ordered_basis_sizes[2])/(temp1/num_faces));
						int temp_index2 = std::floor((k - ordered_basis_sizes[4])/(temp2));
						int index1 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[2])%ndofs[3];
						int index2 = sizes[4] + ndofs[5] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[4])%ndofs[5];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[4]) && 
							(k >= ordered_basis_sizes[4])) {
						int temp1 = ordered_basis_sizes[7] - ordered_basis_sizes[4];
						int temp2 = ordered_basis_sizes[7] - ordered_basis_sizes[4];
						int temp_index1 = std::floor((j - ordered_basis_sizes[4])/(temp1));
						int temp_index2 = std::floor((k - ordered_basis_sizes[4])/(temp2));
						int index1 = sizes[4] + ndofs[5] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[4])%ndofs[5];
						int index2 = sizes[4] + ndofs[5] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[4])%ndofs[5];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[2] && j < ordered_basis_sizes[3]) && 
							((k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1]) ||
							 (k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2]) ||
							 (k >= ordered_basis_sizes[3] && k < ordered_basis_sizes[4]))) {
						int index2 = -1;

						if (k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1]) {
							int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
							int temp_index = std::floor((k - ordered_basis_sizes[0])/(temp/num_edges));
							index2 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][1][temp_index] + (k - ordered_basis_sizes[0])%ndofs[1];
						}
						else if (k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[2]) {
							int temp = ordered_basis_sizes[2] - ordered_basis_sizes[1];
							int temp_index = std::floor((k - ordered_basis_sizes[1])/(temp/num_faces));
							index2 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][2][temp_index] + (k - ordered_basis_sizes[1])%ndofs[2];
						}
						else if (k >= ordered_basis_sizes[3] && k < ordered_basis_sizes[4]) {
							int temp = ordered_basis_sizes[4] - ordered_basis_sizes[3];
							int temp_index = std::floor((k - ordered_basis_sizes[3])/(temp));
							index2 = sizes[3] + ndofs[4] * simplex_sub_simplices[i][3][temp_index] + (k - ordered_basis_sizes[3])%ndofs[4];
						}

						int temp1 = ordered_basis_sizes[3] - ordered_basis_sizes[2];
						int temp_index1 = std::floor((j - ordered_basis_sizes[2])/(temp1/num_faces));
						int index1 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[2])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
				}
			}
		}

		bb_mass_matrices[1].setFromTriplets(triplet.begin(), triplet.end());
		bb_mass_matrices[1].makeCompressed();
	}

	else if (x == 2) {
		VectorTripletD triplet;

		Vector2I alpha;
	Vector2I temp_alpha;
	VectorI ordered_basis_sizes;
	VectorI sizes;

	compute_index_sets_o(alpha,
						 3,
						 3);
	ordered_basis_sizes.push_back(alpha.size());

	size_t total = alpha.size();

	temp_alpha.clear();
	compute_index_sets_p(temp_alpha,
						 n,
						 3);
	size_t temp_alpha_size = temp_alpha.size();
	if (temp_alpha_size != 0) {
		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	}
	ordered_basis_sizes.push_back(total + temp_alpha_size);
	total += temp_alpha_size;

	temp_alpha.clear();
	compute_index_sets_o(temp_alpha,
						 n + 2,
						 4);
	temp_alpha_size = temp_alpha.size();
	for (int l = 0; l < 2; ++l) {	
		if (temp_alpha_size != 0) {
			alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
		}
		ordered_basis_sizes.push_back(total + temp_alpha_size);
		total += temp_alpha_size;
	}
	int counter = 0; 
	for (size_t j = 0; j < temp_alpha_size; ++j) {
		if (temp_alpha[j][2] == 1) {
			alpha.push_back(temp_alpha[j]);
			++counter;
		}
	}
	ordered_basis_sizes.push_back(total + counter);	
	total += counter;

	temp_alpha.clear();
	compute_index_sets_p(temp_alpha,
						 n,
						 4);
	temp_alpha_size = temp_alpha.size();
	if (temp_alpha_size != 0) {
		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	}
	ordered_basis_sizes.push_back(total + temp_alpha_size);
	total += temp_alpha_size;

	size_t alpha_size = alpha.size();

	VectorI ndofs;
	size_t size = 0;

	int temp = 1;
	ndofs.push_back(temp);
	sizes.push_back(size + temp * num_simplices[2]);
	size += temp * num_simplices[2];

	binomialCoeff(temp,
				  n + 2,
				  2);
	temp = temp - 1;
	ndofs.push_back(temp);
	sizes.push_back(size + temp * num_simplices[2]);
	size += temp * num_simplices[2];

	int temp1;
	int temp2;
	binomialCoeff(temp1,
				  n + 1,
				  3);
	binomialCoeff(temp2,
				  n,
				  2);
	temp = 2*temp1 + temp2;
	ndofs.push_back(temp);
	sizes.push_back(size + temp * num_simplices[3]);
	size += temp * num_simplices[3];

	binomialCoeff(temp,
				  n + 3,
				  3);
	temp = temp - 1;
	ndofs.push_back(temp);
	sizes.push_back(size + temp * num_simplices[3]);
	size += temp * num_simplices[3];

		bb_mass_matrices[2].resize(size, size);

		for (size_t i = 0; i < num_simplices[complex_dimension]; ++i) {
			Vector2D pts;
			VectorI simplex = simplices[complex_dimension][i];
			std::sort(simplex.begin(), simplex.end());
			for(size_t k = 0; k < complex_dimension + 1; ++k) {
				pts.push_back(vertices[simplex[k]]);
			}

			DenMatD mass_matrix;
			bb_mass_matrix_H_div(mass_matrix,
								 pts,
						   	   	 n,
						   	   	 alpha,
						   	   	 ordered_basis_sizes);

			for (size_t j = 0; j < alpha_size; ++j) {
				for (size_t k = j; k < alpha_size; ++k) {
					if (j < ordered_basis_sizes[0] && k < ordered_basis_sizes[0]) {
						triplet.push_back(TripletD(simplex_sub_simplices[i][2][j], simplex_sub_simplices[i][2][k], mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(simplex_sub_simplices[i][2][k], simplex_sub_simplices[i][2][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j < ordered_basis_sizes[0]) && 
							(k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1])) {
						int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
						int temp_index = std::floor((k - ordered_basis_sizes[0])/(temp/num_faces));
						int index = sizes[0] + ndofs[1] * simplex_sub_simplices[i][2][temp_index] + (k - ordered_basis_sizes[0])%ndofs[1];
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][2][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][2][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j < ordered_basis_sizes[0]) && 
							(k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[4])) {
						int temp = ordered_basis_sizes[4] - ordered_basis_sizes[1];
						int temp_index = std::floor((k - ordered_basis_sizes[1])/(temp));
						int index = sizes[1] + ndofs[2] * simplex_sub_simplices[i][3][temp_index] + (k - ordered_basis_sizes[1])%ndofs[2];
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][2][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][2][j], mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j < ordered_basis_sizes[0]) && 
							(k >= ordered_basis_sizes[4] && k < ordered_basis_sizes[5])) {
						int temp = ordered_basis_sizes[5] - ordered_basis_sizes[4];
						int temp_index = std::floor((k - ordered_basis_sizes[4])/(temp));
						int index = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index] + (k - ordered_basis_sizes[4])%ndofs[3];
						
						triplet.push_back(TripletD(simplex_sub_simplices[i][2][j], index, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index, simplex_sub_simplices[i][2][j], mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) && 
							(k >= ordered_basis_sizes[0] && k < ordered_basis_sizes[1])) {
						int temp = ordered_basis_sizes[1] - ordered_basis_sizes[0];
						int temp_index1 = std::floor((j - ordered_basis_sizes[0])/(temp/num_faces));
						int temp_index2 = std::floor((k - ordered_basis_sizes[0])/(temp/num_faces));
						int index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[0])%ndofs[1];
						int index2 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][2][temp_index2] + (k - ordered_basis_sizes[0])%ndofs[1];
						
						if(index1 >=22 && index1<27 && index2 >=22 && index2<27) {
							std::cout<<index1<<" "<<index2<<"\n";
						}
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) && 
							(k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[4])) {
						int temp1 = ordered_basis_sizes[1] - ordered_basis_sizes[0];
						int temp2 = ordered_basis_sizes[4] - ordered_basis_sizes[1];
						int temp_index1 = std::floor((j - ordered_basis_sizes[0])/(temp1/num_faces));
						int temp_index2 = std::floor((k - ordered_basis_sizes[1])/(temp2));
						int index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[0])%ndofs[1];
						int index2 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[1])%ndofs[2];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j >= ordered_basis_sizes[0] && j < ordered_basis_sizes[1]) && 
							(k >= ordered_basis_sizes[4] && k < ordered_basis_sizes[5])) {
						int temp1 = ordered_basis_sizes[1] - ordered_basis_sizes[0];
						int temp2 = ordered_basis_sizes[5] - ordered_basis_sizes[4];
						int temp_index1 = std::floor((j - ordered_basis_sizes[0])/(temp1/num_faces));
						int temp_index2 = std::floor((k - ordered_basis_sizes[4])/(temp2));
						int index1 = sizes[0] + ndofs[1] * simplex_sub_simplices[i][2][temp_index1] + (j - ordered_basis_sizes[0])%ndofs[1];
						int index2 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[4])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[4]) && 
							(k >= ordered_basis_sizes[1] && k < ordered_basis_sizes[4])) {
						int temp1 = ordered_basis_sizes[4] - ordered_basis_sizes[1];
						int temp2 = ordered_basis_sizes[4] - ordered_basis_sizes[1];
						int temp_index1 = std::floor((j - ordered_basis_sizes[1])/(temp1));
						int temp_index2 = std::floor((k - ordered_basis_sizes[1])/(temp2));
						int index1 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[1])%ndofs[2];
						int index2 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[1])%ndofs[2];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
					else if ((j >= ordered_basis_sizes[1] && j < ordered_basis_sizes[4]) && 
							(k >= ordered_basis_sizes[4] && k < ordered_basis_sizes[5])) {
						int temp1 = ordered_basis_sizes[4] - ordered_basis_sizes[1];
						int temp2 = ordered_basis_sizes[5] - ordered_basis_sizes[4];
						int temp_index1 = std::floor((j - ordered_basis_sizes[1])/(temp1));
						int temp_index2 = std::floor((k - ordered_basis_sizes[4])/(temp2));
						int index1 = sizes[1] + ndofs[2] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[1])%ndofs[2];
						int index2 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[4])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}

					else if ((j >= ordered_basis_sizes[4] && j < ordered_basis_sizes[5]) && 
							(k >= ordered_basis_sizes[4] && k < ordered_basis_sizes[5])) {
						int temp1 = ordered_basis_sizes[5] - ordered_basis_sizes[4];
						int temp2 = ordered_basis_sizes[5] - ordered_basis_sizes[4];
						int temp_index1 = std::floor((j - ordered_basis_sizes[4])/(temp1));
						int temp_index2 = std::floor((k - ordered_basis_sizes[4])/(temp2));
						int index1 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index1] + (j - ordered_basis_sizes[4])%ndofs[3];
						int index2 = sizes[2] + ndofs[3] * simplex_sub_simplices[i][3][temp_index2] + (k - ordered_basis_sizes[4])%ndofs[3];
						
						triplet.push_back(TripletD(index1, index2, mass_matrix.coeffRef(j, k)));
						if (j != k) {
							triplet.push_back(TripletD(index2, index1, mass_matrix.coeffRef(k, j)));
						}
					}
				}
			}
		}

		bb_mass_matrices[2].setFromTriplets(triplet.begin(), triplet.end());
		bb_mass_matrices[2].makeCompressed();
	}

	return SUCCESS;
}


double FiniteElementExteriorCalculus::bb_error_H_1_global(int n,
												   		  int q_order) {

	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif

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
	VectorI sizes;

	compute_index_sets_o(alpha,
					 	 1,
					 	 1);

	int d = alpha[0].size();
	ordered_basis_sizes.push_back(alpha.size());

	for(int i = 2; i <= std::min(n, d); ++i) {
		temp_alpha.clear();
		compute_index_sets_o(temp_alpha,
						 	 n,
						 	 i);

		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
		ordered_basis_sizes.push_back(alpha.size());
	}
	size_t alpha_size = alpha.size();

	VectorI ndofs;
	size_t size = 0;
	for (size_t i = 0; i < complex_dimension + 1; ++i) {
		int temp;
		binomialCoeff(temp,
					  std::max(0, n-1),
					  i);
		ndofs.push_back(temp);
		sizes.push_back(size + temp * num_simplices[i]);
		size += temp * num_simplices[i];
	}

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

	compute_bb_mass_matrices(0, n);

	EigVectorD b = EigVectorD::Zero(size);

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t i = 0; i < num_simplices[N-1]; ++i) {

		Vector2D pts;
		VectorI simplex = simplices[N-1][i];
		// std::sort(simplex.begin(), simplex.end());
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplex[k]]);
		}
		double vol = get_simplex_volume(pts);

		for(size_t j = 0; j < alpha_size; ++j) {
			double inner_product = 0.0;

			for(size_t k = 0; k < nodes_size; ++k) {
				VectorD vec(embed_dim, 0.0);

				for(size_t v = 0; v < N; ++v) {
					for(size_t l = 0; l < embed_dim; ++l) {
						vec[l] += pts[v][l] * nodes[k][v];
					}
				}

				inner_product += weights[k] * get_analytical_soln(vec) * basis_elements.coeffRef(j, k);
			}

			int index;
			get_global_index(index,
							 0,
							 i,
							 j,
							 ordered_basis_sizes,
							 ndofs,
							 sizes,
							 simplex_sub_simplices,
							 complex_dimension);
			#ifdef MULTICORE
				#pragma omp critical
			#endif
			b.coeffRef(index) += vol * inner_product/sum_weights;
		}
	}

	Eigen::SimplicialLLT<SpMatD> solver;
	EigVectorD coeffs = solver.compute(bb_mass_matrices[0]).solve(b);

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t i = 0; i < num_simplices[N-1]; ++i) {
		double e = 0;

		Vector2D pts;
		VectorI simplex = simplices[N-1][i];
		// std::sort(simplex.begin(), simplex.end());
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplex[k]]);
		}
		double vol = get_simplex_volume(pts);

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {
			double f_dash = 0;
			for (size_t j = 0; j < alpha_size; ++j) {
				int index;
				get_global_index(index,
								 0,
								 i,
								 j,
								 ordered_basis_sizes,
								 ndofs,
								 sizes,
								 simplex_sub_simplices,
								 complex_dimension);
				f_dash += coeffs[index] * basis_elements.coeffRef(j, node_index);
			}
			
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


double FiniteElementExteriorCalculus::bb_error_H_curl_global(int n,
												   		     int q_order) {

	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif

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
	VectorI sizes;

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
			temp_alpha.clear();
			compute_index_sets_o(temp_alpha,
								 n + 2,
								 i);
			
			size_t temp_alpha_size = temp_alpha.size();

			for (int l = 0; l < 2; ++l) {	
				if (temp_alpha_size != 0) {
					alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
				}
				ordered_basis_sizes.push_back(total + temp_alpha_size);
				total += temp_alpha_size;
			}

			size_t counter = 0; 
			for (size_t j = 0; j < temp_alpha_size; ++j) {
				if (temp_alpha[j][2] == 1) {
					alpha.push_back(temp_alpha[j]);
					++counter;
				}
			}

			ordered_basis_sizes.push_back(total + counter);	
			total += counter;
		}
	}

	size_t alpha_size = alpha.size();

	VectorI ndofs;
	size_t size = 0;
	for (size_t i = 0; i < complex_dimension + 1; ++i) {
		int temp;
		binomialCoeff(temp,
					  n,
					  i);
		ndofs.push_back(temp);
		sizes.push_back(size + temp * num_simplices[std::max(1, (int)i)]);
		size += temp * num_simplices[std::max(1, (int)i)];

		if (i == 2) {
			temp = 0;
			binomialCoeff(temp,
						  n+2,
						  2);
			temp = temp - 1;
			ndofs.push_back(temp);
			sizes.push_back(size + temp * num_simplices[std::max(1, (int)i)]);
			size += temp * num_simplices[std::max(1, (int)i)];
		}

		else if (i == 3) {
			int temp1;
			int temp2;
			binomialCoeff(temp1,
						  n+1,
						  3);
			binomialCoeff(temp2,
						  n,
						  2);
			temp = 2*temp1 + temp2;
			ndofs.push_back(temp);
			sizes.push_back(size + temp * num_simplices[std::max(1, (int)i)]);
			size += temp * num_simplices[std::max(1, (int)i)];
		}
	}

	compute_bb_mass_matrices(1, n);

	EigVectorD b = EigVectorD::Zero(size);
	Vector2DenMatD global_basis_elements;

	Vector2I all_faces;
	compute_index_sets_o(all_faces,
						 3,
						 3);

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t s = 0; s < num_simplices[N-1]; ++s) {

		Vector2D pts;
		VectorI simplex = simplices[N-1][s];
		std::sort(simplex.begin(), simplex.end());
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplex[k]]);
		}
		double vol = get_simplex_volume(pts);

		DenMatD grad_bary_coords;
		barycentric_gradients(grad_bary_coords,
							  pts);

		VectorDenMatD basis_elements;
		for(int i = 0; i < (int)alpha_size; ++i) {
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

		global_basis_elements.push_back(basis_elements);

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
				inner_product += weights[k] * f.dot(basis_elements[j].row(k));
			}

			int index;
			get_global_index(index,
							 1,
							 s,
							 j,
							 ordered_basis_sizes,
							 ndofs,
							 sizes,
							 simplex_sub_simplices,
							 complex_dimension);
			#ifdef MULTICORE
				#pragma omp critical
			#endif
			b.coeffRef(index) += vol * inner_product/sum_weights;
		}
	}

	Eigen::SimplicialLLT<SpMatD> solver;
	EigVectorD coeffs = solver.compute(bb_mass_matrices[1]).solve(b);

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t s = 0; s < num_simplices[N-1]; ++s) {
		double e = 0;

		Vector2D pts;
		VectorI simplex = simplices[N-1][s];
		std::sort(simplex.begin(), simplex.end());
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplex[k]]);
		}
		double vol = get_simplex_volume(pts);

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {

			EigVectorD f_dash = EigVectorD::Zero(embed_dim);
			for (size_t j = 0; j < alpha_size; ++j) {
				int index;
				get_global_index(index,
								 1,
								 s,
								 j,
								 ordered_basis_sizes,
								 ndofs,
								 sizes,
								 simplex_sub_simplices,
								 complex_dimension);
				f_dash += coeffs.coeffRef(index) * global_basis_elements[s][j].row(node_index);
			}

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


double FiniteElementExteriorCalculus::bb_error_H_div_global(int n,
												   		     int q_order) {

	#ifdef PYTHON
		pybind11::gil_scoped_acquire acquire;
	#endif

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
	VectorI sizes;

	compute_index_sets_o(alpha,
						 3,
						 3);
	ordered_basis_sizes.push_back(alpha.size());

	size_t total = alpha.size();

	temp_alpha.clear();
	compute_index_sets_p(temp_alpha,
						 n,
						 3);
	size_t temp_alpha_size = temp_alpha.size();
	if (temp_alpha_size != 0) {
		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	}
	ordered_basis_sizes.push_back(total + temp_alpha_size);
	total += temp_alpha_size;

	temp_alpha.clear();
	compute_index_sets_o(temp_alpha,
						 n + 2,
						 4);
	temp_alpha_size = temp_alpha.size();
	for (int l = 0; l < 2; ++l) {	
		if (temp_alpha_size != 0) {
			alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
		}
		ordered_basis_sizes.push_back(total + temp_alpha_size);
		total += temp_alpha_size;
	}
	int counter = 0; 
	for (size_t j = 0; j < temp_alpha_size; ++j) {
		if (temp_alpha[j][2] == 1) {
			alpha.push_back(temp_alpha[j]);
			++counter;
		}
	}
	ordered_basis_sizes.push_back(total + counter);	
	total += counter;

	temp_alpha.clear();
	compute_index_sets_p(temp_alpha,
						 n,
						 4);
	temp_alpha_size = temp_alpha.size();
	if (temp_alpha_size != 0) {
		alpha.insert(alpha.end(), temp_alpha.begin(), temp_alpha.end());
	}
	ordered_basis_sizes.push_back(total + temp_alpha_size);
	total += temp_alpha_size;

	size_t alpha_size = alpha.size();

	VectorI ndofs;
	size_t size = 0;

	int temp = 1;
	ndofs.push_back(temp);
	sizes.push_back(size + temp * num_simplices[2]);
	size += temp * num_simplices[2];

	binomialCoeff(temp,
				  n + 2,
				  2);
	temp = temp - 1;
	ndofs.push_back(temp);
	sizes.push_back(size + temp * num_simplices[2]);
	size += temp * num_simplices[2];

	int temp1;
	int temp2;
	binomialCoeff(temp1,
				  n + 1,
				  3);
	binomialCoeff(temp2,
				  n,
				  2);
	temp = 2*temp1 + temp2;
	ndofs.push_back(temp);
	sizes.push_back(size + temp * num_simplices[3]);
	size += temp * num_simplices[3];

	binomialCoeff(temp,
				  n + 3,
				  3);
	temp = temp - 1;
	ndofs.push_back(temp);
	sizes.push_back(size + temp * num_simplices[3]);
	size += temp * num_simplices[3];

	compute_bb_mass_matrices(2, n);

	EigVectorD b = EigVectorD::Zero(size);
	Vector2DenMatD global_basis_elements;

	Vector2I all_faces;
	compute_index_sets_o(all_faces,
						 3,
						 3);

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t s = 0; s < num_simplices[N-1]; ++s) {

		Vector2D pts;
		VectorI simplex = simplices[N-1][s];
		std::sort(simplex.begin(), simplex.end());
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplex[k]]);
		}
		double vol = get_simplex_volume(pts);

		DenMatD grad_bary_coords;
		barycentric_gradients(grad_bary_coords,
							  pts);

		VectorDenMatD basis_elements;
		for(int i = 0; i < (int)alpha_size; ++i) {
			DenMatD temp_basis_elements(nodes_size, embed_dim);
			
			for(size_t j = 0; j < nodes_size; ++j) {
				if (i < ordered_basis_sizes[0]) {
					EigVectorD chi;

					VectorI local_indices;
					for (int j = 0; j < 4; ++j) {
						if (alpha[i][j] > 0) {
							local_indices.push_back(j);
						}
					}
					
					chi_l(chi,
						  alpha[i],
						  nodes[j],
						  grad_bary_coords,
						  local_indices);
					
					temp_basis_elements.row(j) = chi;
				}
				else if (i >= ordered_basis_sizes[0] && i < ordered_basis_sizes[1]) {
					EigVectorD curl_phi;

					int E_nF_size = ordered_basis_sizes[1] - ordered_basis_sizes[0];
					int face_index = std::floor((i - ordered_basis_sizes[0])/(E_nF_size/4));
					VectorI face = all_faces[face_index];
					
					VectorI temp_local_indices;

					for (int k = 0; k < 4; ++k) {
						if (face[k] > 0) {
							temp_local_indices.push_back(k);
						}
					}

					curl_phi_FT(curl_phi,
							   alpha[i],
							   n,
							   nodes[j],
							   grad_bary_coords,
							   temp_local_indices);

					temp_basis_elements.row(j) = curl_phi;
				}
				else if (i >= ordered_basis_sizes[1] && i < ordered_basis_sizes[2]) {
					EigVectorD curl_psi;

					curl_psi_T(curl_psi,
							   alpha[i],
							   n + 1,
							   0,
							   nodes[j],
							   grad_bary_coords);

					temp_basis_elements.row(j) = curl_psi;
				}
				else if (i >= ordered_basis_sizes[2] && i < ordered_basis_sizes[3]) {
					EigVectorD curl_psi;

					curl_psi_T(curl_psi,
							   alpha[i],
							   n + 1,
							   1,
							   nodes[j],
							   grad_bary_coords);

					temp_basis_elements.row(j) = curl_psi;
				}
				else if (i >= ordered_basis_sizes[3] && i < ordered_basis_sizes[4]) {
					EigVectorD curl_psi;

					curl_psi_T(curl_psi,
							   alpha[i],
							   n + 1,
							   2,
							   nodes[j],
							   grad_bary_coords);

					temp_basis_elements.row(j) = curl_psi;
				}
				else if (i >= ordered_basis_sizes[4] && i < ordered_basis_sizes[5]) {
					EigVectorD upsi;

					upsilon(upsi,
						    alpha[i],
						    n,
						    nodes[j],
						    grad_bary_coords);

					temp_basis_elements.row(j) = upsi;
				}
			}

			basis_elements.push_back(temp_basis_elements);
		}

		#ifdef MULTICORE
			#pragma omp critical
		#endif
		global_basis_elements.push_back(basis_elements);

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
				inner_product += weights[k] * f.dot(basis_elements[j].row(k));
			}

			int index;
			get_global_index(index,
							 2,
							 s,
							 j,
							 ordered_basis_sizes,
							 ndofs,
							 sizes,
							 simplex_sub_simplices,
							 complex_dimension);
			#ifdef MULTICORE
				#pragma omp critical
			#endif
			b.coeffRef(index) += vol * inner_product/sum_weights;
		}
	}

	Eigen::SimplicialLLT<SpMatD> solver;
	EigVectorD coeffs = solver.compute(bb_mass_matrices[2]).solve(b);
	// std::cout<<coeffs<<"\n";

	#ifdef MULTICORE
		#pragma omp parallel for
	#endif
	for(size_t s = 0; s < num_simplices[N-1]; ++s) {
		double e = 0;

		Vector2D pts;
		VectorI simplex = simplices[N-1][s];
		std::sort(simplex.begin(), simplex.end());
		for(size_t k = 0; k < N; ++k) {
			pts.push_back(vertices[simplex[k]]);
		}
		double vol = get_simplex_volume(pts);

		for(size_t node_index = 0; node_index < nodes_size; ++node_index) {

			EigVectorD f_dash = EigVectorD::Zero(embed_dim);
			for (size_t j = 0; j < alpha_size; ++j) {
				int index;
				get_global_index(index,
								 2,
								 s,
								 j,
								 ordered_basis_sizes,
								 ndofs,
								 sizes,
								 simplex_sub_simplices,
								 complex_dimension);
				f_dash += coeffs.coeffRef(index) * global_basis_elements[s][j].row(node_index);
			}

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

    set_mass_matrices_to_null();
    set_bb_mass_matrices_to_null();
    
}


FiniteElementExteriorCalculus::FiniteElementExteriorCalculus(SimplicialComplex sc) : GeometryComplex(sc) {

    set_mass_matrices_to_null();
    set_bb_mass_matrices_to_null();
}

FiniteElementExteriorCalculus::~FiniteElementExteriorCalculus() {}