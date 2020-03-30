#include <iostream>
#include <vector>
#include <cmath>

#define _USE_MATH_DEFINES

extern "C" double function(std::vector < double > f) {
	double x = 0.0;
	for(int i = 0; i < f.size(); ++i) {
		// x += sin(f[i]*M_PI);
		// x += f[i]*f[i];
		x += (i+1)*f[i];
	}

	// x += exp(-1 * (pow(f[0] - 0.5, 2) * pow(f[1] - 0.5, 2) * pow(f[2] - 0.5, 2)));
	x = f[0]*f[1];
	// x = sqrt(x);

	return x;
}