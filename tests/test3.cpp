#include "simplicial_complex.h"
#include "core_utils.h"
#include "finite_element_exterior_calculus.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	Vector2I sets;
	SimplicialComplex sc;
	sc.build_complex();
	
	FiniteElementExteriorCalculus fem(sc);
	DenMatD mass_matrix;
	fem.mass_matrix_bb_0(mass_matrix,
						1,
						1);
	return 0;
}



