#ifndef FINITE_ELEMENT_EXTERIOR_CALCULUS_H
#define FINITE_ELEMENT_EXTERIOR_CALCULUS_H

#include "simplicial_complex.h"
#include "geometry.h"
#include "definitions.h"
#include "discrete_exterior_calculus.h"

class FiniteElementExteriorCalculus: public GeometryComplex {

 public:
    VectorSpmatD hodge_stars;

    bool all_hodge_stars;

 public:
    int get_hodge_star(int k, SpMatD &mat);

 public:
    FiniteElementExteriorCalculus();

    FiniteElementExteriorCalculus(SimplicialComplex sc);

    ~FiniteElementExteriorCalculus();

    int compute_hodge_stars();

    int compute_hodge_star_k(int &k);

    int set_hodge_stars_to_null();

    int bb_mass_matrix_H_1(DenMatD &mass_matrix,
    					   int n, 
    					   int m,
    					   Vector2I &index_sets,
    					   int d = 3);

    int bb_mass_matrix_H_curl(DenMatD &mass_matrix,
    						  Vector2D &pts,
							  int n,
							  Vector2I &alpha,
							  VectorI &ordered_basis_sizes);

    int M_alpha_beta(double &M,
			  		VectorI &alpha,
			  		VectorI &beta);

    int S_n(DenMatD &S,
			Vector2D &points,
	  		int n);
    
    int bernstein(double &bernstein_poly,
				 VectorI &alpha,
				 int n,
				 VectorD &bary_coords);

    int omega_ij(EigVectorD &omega,
				VectorD &bary_coords,
				DenMatD &grad_bary_coords);

    int grad_B(EigVectorD &grad_b,
			  VectorI &alpha,
			  int n,
			  VectorD &bary_coords,
			  DenMatD &grad_bary_coords);

    int phi_FT(EigVectorD &phi,
			  VectorI &alpha,
			  int n,
			  VectorD &bary_coords,
			  DenMatD &grad_bary_coords,
			  VectorI &local_indices);

    int psi_T(EigVectorD &psi,
			  VectorI &alpha,
			  int n,
			  int l,
			  VectorD &bary_coords,
			  DenMatD &grad_bary_coords);

    int compute_index_sets_o(Vector2I &sets,
				   			 int sum,
				   			 int dim,
				   			 int d = 3);

    int compute_index_sets_t(Vector2I &sets,
				   			 int sum,
				   			 int dim,
				   			 int d = 3);

    int compute_index_sets_p(Vector2I &sets,
				   			 int sum,
				   			 int dim,
				   			 int d = 3);

    double bb_error_H_1(int n,
				    	Vector3I &simplices,
						Vector2D &vertices,
						VectorI &num_simplices,
						int q_order = 4);

    double bb_error_H_curl(int n,
				      	   Vector3I &simplices,
					  	   Vector2D &vertices,
						   VectorI &num_simplices,
						   int q_order = 4);

    double bb_error_H_curl_1d_quad(int n,
							       Vector3I &simplices,
								   Vector2D &vertices,
								   VectorI &num_simplices,
								   int q_order = 4);

    int test_basis_functions(int n);
    double test_mass_matrix(int n,
    					 int q_order = 4);
};

#endif
