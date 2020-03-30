#include "simplicial_complex.h"
#include "core_utils.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	VectorD f {1,2,3};
	std::cout<<get_analytical_soln(f);

	return 0;
}



