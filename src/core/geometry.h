#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "simplicial_complex.h"
#include "definitions.h"
 
class GeometryComplex : public SimplicialComplex {
public:
    Vector2D primal_volume;
    Vector2D dual_volume;

public:
    GeometryComplex();

    GeometryComplex(SimplicialComplex sc);

    ~GeometryComplex();

public:

	int compute_primal_volumes();

  int compute_dual_volumes();

  double compute_primal_volume_k(int &dim,
                              int &k);

  int compute_dual_volume_k(int &dim,
                            int &k);

  int set_volumes_to_null();

  std::tuple<Vector2D, DenMatD> simplex_quivers(VectorD form);
    
};

#endif

