#include "simplicial_complex.h"
#include "geometry.h"
#include "core_utils.h"
#include "discrete_exterior_calculus.h"
#include "definitions.h"
#include "finite_element_exterior_calculus.h"
#include <iostream>
#include <vector>

int main (int argc, char const *argv[]) {
	DiscreteExteriorCalculus dec;
	dec.compute_dual_volumes();

	print_vector(dec.dual_volume);

	return SUCCESS;
}