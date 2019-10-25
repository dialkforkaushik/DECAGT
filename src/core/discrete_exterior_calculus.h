#ifndef DISCRETE_EXTERIOR_CALCULUS_H
#define DISCRETE_EXTERIOR_CALCULUS_H

#include "simplicial_complex.h"
#include "geometry.h"
#include "definitions.h"

class DiscreteExteriorCalculus: public GeometryComplex {

 public:
    VectorSpmatD hodge_stars;
    
    bool all_hodge_stars;

 public:

    int compute_hodge_stars();

    int compute_hodge_star_k(int &k);

    DiscreteExteriorCalculus();

    DiscreteExteriorCalculus(SimplicialComplex sc);
    
    ~DiscreteExteriorCalculus();

    int set_hodge_stars_to_null();
};

#endif
