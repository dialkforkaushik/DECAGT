#include "simplicial_complex.h"
#include "core_utils.h"
#include "finite_element_exterior_calculus.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	double error;
	VectorD U;
	SimplicialComplex sc;
	sc.build_complex();

	FiniteElementExteriorCalculus fem(sc);

	// for(int i = 0; i < sc.vertices.size(); ++i) {
	// 	U.push_back(get_analytical_soln(sc.vertices[i]));
	// }

	// DenMatD M;
	// fem.mass_matrix_bb_0(M,5,5);
	// std::cout<<M<<"\n\n";

	for(int i = 1; i < 21; ++i) {
		error = fem.bb_error(i,
							sc.simplices,
							sc.vertices,
							sc.num_simplices,
							4);

		std::cout<<error<<"\n";
	}

	return 0;
}



