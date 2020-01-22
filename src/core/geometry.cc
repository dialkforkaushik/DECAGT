#include "definitions.h"
#include "simplicial_complex.h"
#include "core_utils.h"
#include "geometry.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <Eigen/Eigen>
#include <set>



int barycentric_gradients(Vector2D &pts,
						  DenMatD &X) {

	int cols = pts[0].size();
	int rows = pts.size();

	DenMatD V(rows - 1, cols);
	MapEigVectorD v0(pts[0].data(), cols);
	
	for (int i = 1; i < rows; ++i) {
		MapEigVectorD v(pts[i].data(), cols);
		V.row(i - 1) = v - v0;
	}

	DenMatD temp_grads = V.completeOrthogonalDecomposition().pseudoInverse();
	DenMatD grads = temp_grads.transpose();

	DenMatD mat = DenMatD::Ones(1, grads.rows());
	EigVectorD vec = (mat * grads).row(0);

	X = DenMatD::Zero(grads.rows() + 1, grads.cols());
	X.row(0) = -1 * vec;
	X.bottomRows(grads.rows()) = grads;

	return SUCCESS;
}


int GeometryComplex::set_volumes_to_null() {

    for (int i = 0; i < complex_dimension + 1; ++i) {
    	primal_volume.push_back(VectorD());
    	dual_volume.push_back(VectorD());

    	for (int j = 0; j < num_simplices[i]; ++j) {
    		primal_volume[i].push_back(0.0);
    		dual_volume[i].push_back(0.0);
    	}
    }

    return SUCCESS;
}


double GeometryComplex::compute_primal_volume_k(int &dim,
											    int &k) {
	double primal_volume_k;

	if (dim == 0) {
		primal_volume_k = 1.0;
		return SUCCESS;
	}

	int row = num_simplices[dim];
	int col = simplices[dim][0].size();

	Vector2D pts;
	pts.reserve(col);
	for (int i = 0; i < col; ++i) {
		pts.push_back(vertices[simplices[dim][k][i]]);
	}

	if (row - 1 == col) {
		signed_volume(pts, primal_volume_k);
	}
	else {
		unsigned_volume(pts, primal_volume_k);
	}

	primal_volume[dim][k] = primal_volume_k;

	return SUCCESS;
}



int GeometryComplex::compute_dual_volume_k(int &dim,
										   int &k) {

	size_t N = complex_dimension + 1;

	Vector2D temp_centers;

	calculate_dual_volume(dual_volume,
						  temp_centers,
						  vertices,
						  simplices,
						  adjacency1d,
						  simplices[N-1][k],
						  N-1,
						  k);
	
	temp_centers.clear();
	
	return SUCCESS;
}


int GeometryComplex::compute_primal_volumes() {
	size_t N = complex_dimension + 1;

	primal_volume.push_back(VectorD());
	for (int i = 0; i < num_simplices[0]; ++i) {
		primal_volume[0].push_back(1.0);
	}

	for (int i = 1; i < N; ++i) {
		primal_volume.push_back(VectorD());
		int row = num_simplices[i];
		int col = simplices[i][0].size();

		for (int j = 0; j < row; ++j) {
			Vector2D pts;
			for (int k = 0; k < col; ++k) {
				pts.push_back(vertices[simplices[i][j][k]]);
			}
			double vol;
			if (row - 1 == col) {
				signed_volume(pts, vol);
			}
			else {
				unsigned_volume(pts, vol);
			}
			
			primal_volume[i].push_back(vol);
		}
	}

	return SUCCESS;
}



int GeometryComplex::compute_dual_volumes() {
	size_t N = complex_dimension + 1;

	Vector2D temp_centers;

	#pragma omp parallel for private(temp_centers)
	for (int j = 0; j < num_simplices[N-1]; ++j) {
		calculate_dual_volume(dual_volume,
							  temp_centers,
							  vertices,
							  simplices,
							  adjacency1d,
							  simplices[N-1][j],
							  N-1,
							  j);
		
		temp_centers.clear();
	}
	
	return SUCCESS;
}


std::tuple<Vector2D, DenMatD> GeometryComplex::simplex_quivers(VectorD form) {

	Vector2D quiver_bases;
	DenMatD quiver_dirs;

	VectorD sum(embedding_dimension, 0);
	for (int i = 0; i < num_simplices[complex_dimension]; ++i) {
		std::fill(sum.begin(), sum.end(), 0);
		for (int j = 0; j < complex_dimension + 1; ++j) {
			// sumvertices[simplices[complex_dimension][i][j]];
			for (int k = 0; k < embedding_dimension; ++k) {
				sum[k] += vertices[simplices[complex_dimension][i][j]][k];
			}
		}
		for (int j = 0; j < embedding_dimension; ++j) {
			sum[j] = sum[j]/(complex_dimension + 1);
		}

		quiver_bases.push_back(sum);
	}

    quiver_dirs = DenMatD::Zero(num_simplices[complex_dimension],embedding_dimension);

    for (int i = 0; i < num_simplices[complex_dimension]; ++i) {
    	DenMatD d_lambda;
    	Vector2D pts;
    	for (int j = 0; j < complex_dimension + 1; ++j) {
    		pts.push_back(vertices[simplex_sorted[i][j]]);
    	}
    	barycentric_gradients(pts, d_lambda);

    	Vector2I edges;
    	get_combinations_simplex(simplex_sorted[i], edges, 2);

    	VectorI indices;
    	for (int j = 0; j < edges.size(); ++j) {
    		indices.push_back(elements[edges[j]]);
    	}

    	VectorD values;
    	for (int j = 0; j < indices.size(); ++j) {
    		values.push_back(form[indices[j]]);
    	}

    	Vector2I combinations;
    	VectorI range;
    	for (int j = 0; j < simplex_sorted[i].size(); ++j) {
    		range.push_back(j);
    	}
    	get_combinations_simplex(range, combinations, 2);
    	for (int j = 0; j < values.size(); ++j) {
    		quiver_dirs.row(i) += values[j] * (d_lambda.row(combinations[j][1]) - d_lambda.row(combinations[j][0]));
    	}

    }

    quiver_dirs /= (complex_dimension + 1);

    return std::make_tuple(quiver_bases, quiver_dirs);
}


GeometryComplex::GeometryComplex() {
	
	build_complex();
	set_volumes_to_null();
}

GeometryComplex::GeometryComplex(SimplicialComplex sc) : SimplicialComplex(sc.vertices,
																		   sc.simplex) {
	build_complex();
	set_volumes_to_null();
}

GeometryComplex::~GeometryComplex() {
}