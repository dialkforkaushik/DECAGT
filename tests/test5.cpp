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

	for (int i = 0; i < 11; ++i) {
		std::cout<<"\nPolynomial Degree: "<<i<<"\n";
		error = fem.bb_error_H_div(i,
							       i+3);

		std::cout<<"\nError: "<<error<<"\n";
	}

	// EigVectorD c;
	// EigVectorD v1(3);
	// v1 << -1, -1, -1;
	// EigVectorD v2(3);
	// v2 << 0, 1, 0;

	// cross_product(c, v1, v2);
	// std::cout<<c;

	return SUCCESS;
}