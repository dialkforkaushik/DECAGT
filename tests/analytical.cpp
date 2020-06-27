#include "simplicial_complex.h"
#include "core_utils.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	VectorD f {1.1,2.54,4.3};
	VectorD g;
	get_analytical_soln_vec(g,
							f);

	print_vector(g);
	std::cout<<g[0];

	return 0;
}



