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
	fem.compute_index_sets_p(sets,
							 4,
							 4);
	print_vector(sets);
	std::cout<<sets.size();
	return 0;
}