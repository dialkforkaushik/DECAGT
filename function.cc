#include <iostream>
#include <vector>
#include <cmath>

#define _USE_MATH_DEFINES

extern "C" double function(std::vector < double > f) {
	double out = 0.0;
	for(int i = 0; i < f.size(); ++i) {
		out += sin(f[i]*M_PI);
		// out += f[i]*f[i];
		// out += f[i];
	}

	// out += exp(-1 * (pow(f[0] - 0.5, 2) * pow(f[1] - 0.5, 2) * pow(f[2] - 0.5, 2)));
	// out = f[0]*f[1];
	// out = sqrt(out);

	return out;
}