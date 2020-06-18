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

	// DenMatD M;
	// fem.bb_mass_matrix_H_curl(M,U,1);
	// std::cout<<M<<"\n\n";

	for (int i = 1; i < 11; ++i) {
		std::cout<<"\nPolynomial Degree: "<<i<<"\n";
		error = fem.bb_error_H_curl(i,
							     sc.simplices,
							     sc.vertices,
							     sc.num_simplices,
							     4);

		std::cout<<"\nError: "<<error<<"\n";
	}

	for (int i = 1; i < 11; ++i) {
		std::cout<<"\nPolynomial Degree: "<<i<<"\n";
		error = fem.bb_error_H_curl(i,
							     sc.simplices,
							     sc.vertices,
							     sc.num_simplices,
							     i);

		std::cout<<"\nError: "<<error<<"\n";
	}

	for (int i = 1; i < 11; ++i) {
		std::cout<<"\nPolynomial Degree: "<<i<<"\n";
		error = fem.bb_error_H_curl(i,
							     sc.simplices,
							     sc.vertices,
							     sc.num_simplices,
							     i+1);

		std::cout<<"\nError: "<<error<<"\n";
	}

	return 0;
}



