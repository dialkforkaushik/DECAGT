#include "simplicial_complex.h"
#include "core_utils.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	SimplicialComplex sc;
	sc.compute_simplices();
	sc.compute_boundary_matrices();

	return 0;
}



