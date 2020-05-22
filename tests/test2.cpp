#include "simplicial_complex.h"
#include "core_utils.h"
#include "definitions.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	double error;
	VectorD U;
	SimplicialComplex sc;
	sc.build_complex();

	for(int i = 0; i < sc.vertices.size(); ++i) {
		U.push_back(get_analytical_soln(sc.vertices[i]));
	}

	error = quadratic_error_0(U,
			1,
			sc.simplices,
			sc.vertices,
			sc.num_simplices);

	std::cout<<error<<"\n";

	return 0;
}



