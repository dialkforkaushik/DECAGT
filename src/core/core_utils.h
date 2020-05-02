#ifndef UTILITIES_H
#define UTILITIES_H

#include "definitions.h"

int factorial (int &n);

int binomialCoeff(int &nCk,
				  int n,
				  int k); 

int count_columns(std::string &line,
				  int &columns);

int get_combinations_simplex(VectorI &simplex,
							 Vector2I &vector,
							 int k = -1);

int get_combinations_simplex(Vector2I &simplex,
							 Vector3I &vector,
							 int k = -1);

int get_combinations_mesh(Vector2I &simplex,
						  Vector2I &vector,
						  int k = -1);

int is_rotated(VectorI v1,
               VectorI v2);

int signed_volume(Vector2D &pts,
				  double &vol);

int unsigned_volume(Vector2D &pts,
					double &vol);

int calculate_dual_volume(Vector2D &dual_volume,
						  Vector2D &temp_centers,
						  Vector2D &vertices,
						  Vector3I &simplices,
						  VectorMap3I &adjacency1d,
						  VectorI &pts,
						  VectorD &highest_dim_circumcenter,
						  int dim,
						  int index,
						  size_t complex_dimension);

int get_circumcenter(VectorD &center,
				     double &radius,
				     Vector2D &pts);

int k_cochain_norm(VectorD U,
				   std::string weight,
				   VectorSpmatD hodge_stars,
				   int k);

int get_closest_simplex_to_point(VectorD &point,
							 VectorI &simplex,
							 Vector2D &vertices,
							 Vector3I &simplices,
							 VectorMap3I &adjacency1d,
							 size_t embedding_dimension,
							 size_t complex_dimension);

double get_analytical_soln(VectorD &vec);

double get_simplex_volume(Vector2D &vertices);

int get_index_sets(Vector2I &sets,
				   int sum,
				   int dim,
				   int d = 3);

double error_0(VectorD &U,
				int q_order,
				Vector3I &simplices,
				Vector2D &vertices,
				VectorI &num_simplices);

double quadratic_error_0(VectorD &U,
						int q_order,
						Vector3I &simplices,
						Vector2D &vertices,
						VectorI &num_simplices);

double quadratic_error_0_bb(VectorD &U,
							int q_order,
							Vector3I &simplices,
							Vector2D &vertices,
							VectorI &num_simplices);

double quadratic_error_0_bb_mass(VectorD &U,
								int q_order,
								Vector3I &simplices,
								Vector2D &vertices,
								VectorI &num_simplices);

int print_vector(Vector2D &vec);

int print_vector(Vector2I &vec);

int print_vector(Vector3I &vec);

int print_vector(VectorD &vec);

int print_vector(VectorI &vec);

int print_vector(VectorMap3I &vec);

#endif
