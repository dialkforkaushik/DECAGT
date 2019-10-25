#include "definitions.h"
#include "simplicial_complex.h"
#include "core_utils.h"
#include "discrete_exterior_calculus.h"
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
#include <cmath>
#include <limits>


int DiscreteExteriorCalculus::set_hodge_stars_to_null() {

    for (int i = 0; i < complex_dimension + 1; ++i) {
        SpMatD matrix (simplices[i].size(), simplices[i].size());

        for (int j = 0; j < simplices[i].size(); j++) {
            matrix.coeffRef(j, j) = std::numeric_limits<double>::quiet_NaN();
        }

        hodge_stars.push_back(matrix);
    }

    return SUCCESS;
}


int DiscreteExteriorCalculus::compute_hodge_stars() {

    for (int i = 0; i < complex_dimension + 1; ++i) {
    	for (int j = 0; j < simplices[i].size(); j++) {
    		if (std::isnan(hodge_stars[i].coeffRef(j, j)))
    			hodge_stars[i].coeffRef(j, j) = dual_volume[i][j] / primal_volume[i][j];
    	}
    }

    return SUCCESS;
}


int DiscreteExteriorCalculus::compute_hodge_star_k(int &k) {

    for (int j = 0; j < num_simplices[k]; j++) {
        if (std::isnan(hodge_stars[k].coeffRef(j, j)))
            hodge_stars[k].coeffRef(j, j) = dual_volume[k][j] / primal_volume[k][j];
    }

    return SUCCESS;
}



DiscreteExteriorCalculus::DiscreteExteriorCalculus() {

    set_hodge_stars_to_null();

}


DiscreteExteriorCalculus::DiscreteExteriorCalculus(SimplicialComplex sc) : GeometryComplex(sc) {

    set_hodge_stars_to_null();
}


DiscreteExteriorCalculus::~DiscreteExteriorCalculus() {}
