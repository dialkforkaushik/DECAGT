#include "simplicial_complex.h"
#include "geometry.h"
#include "core_utils.h"
#include "discrete_exterior_calculus.h"
#include "finite_element_exterior_calculus.h"
#include "definitions.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	double error;
	SimplicialComplex sc;
	sc.build_complex();

	FiniteElementExteriorCalculus fem(sc);

	// fem.test_basis_functions(3);
	error = fem.test_mass_matrix(3, 1);
	std::cout<<"\nerror: "<<error<<"\n";


	return SUCCESS;
}