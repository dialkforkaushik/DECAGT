#include "simplicial_complex.h"
#include "geometry.h"
#include "core_utils.h"
#include "discrete_exterior_calculus.h"
#include "definitions.h"
#include "finite_element_exterior_calculus.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	SimplicialComplex sc;
	FiniteElementExteriorCalculus fem(sc);

	fem.compute_bb_mass_matrices(0, 4);
	std::cout<<fem.bb_mass_matrices[0]<<"\n";

	return SUCCESS;
}