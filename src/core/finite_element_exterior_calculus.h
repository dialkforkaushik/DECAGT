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

    int mass_matrix_bb_0(DenMatD &mass_matrix, 
    					 int n, 
    					 int m,
    					 int d = 3);
    
    int bb_basis(double &bernstein_poly,
				 VectorI &alpha,
				 int n,
				 VectorD &bary_coords);

    int omega_ij(double &omega,
				VectorD &bary_coords,
				VectorD &grad_bary_coords);

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
};

#endif
