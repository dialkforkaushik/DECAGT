#ifndef SIMPLICIAL_COMPLEX_H
#define SIMPLICIAL_COMPLEX_H

#include "definitions.h"


class SimplicialComplex {
public:
    Vector2D vertices;
    Vector2I simplex;
    Vector2I simplex_sorted;
    size_t complex_dimension;
    size_t embedding_dimension;
 
    VectorI num_simplices;
    
    VectorSpmatI boundary_matrices;
    VectorSpmatIC boundary_matrices_col_major;
    Vector3I simplices;
    MapVectorI elements;
    Vector3I simplex_sub_simplices;
    VectorMap3I adjacency1d;
    VectorMap3I adjacency2d;

public:
    SimplicialComplex();

    SimplicialComplex(Vector2D &v, 
                      Vector2I &f);
	
    ~SimplicialComplex();

    int read_files(int &rows,
                   VectorS &files);
        
    int get_complex_dim();
    
    int get_embedding_dim();
    
    int compute_simplices();
    
    int compute_boundary_matrices();

    int compute_adjacency1d();

    int compute_adjacency2d();

    int compute_elements();

    int compute_simplex_sub_simplices();

    int build_complex();

    Vector2D circumcenter(int dim);
};

#endif




