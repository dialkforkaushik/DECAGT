#include "simplicial_complex.h"
#include "core_utils.h"
#include "finite_element_exterior_calculus.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	double error;
	SimplicialComplex sc;
	sc.build_complex();

	FiniteElementExteriorCalculus fem(sc);

	// for(int i = 0; i < 4; ++i) {
	// 	U.push_back(sc.vertices[i]);
	// }

	error = fem.bb_error_stiffness_H_curl(1,1);
	std::cout<<error<<"\n\n";

	return 0;
}



