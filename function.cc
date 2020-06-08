#include <iostream>
#include <vector>
#include <cmath>

#define _USE_MATH_DEFINES

// extern "C" double function(std::vector < double > f) {
// 	double out = 0.0;
// 	for(int i = 0; i < f.size(); ++i) {
// 		// out += sin(f[i]*M_PI);
// 		// out += f[i]*f[i];
// 		out += f[i];
// 		// out += f[i]*f[i]*f[i];
// 	}

// 	return out;
// }

extern "C" std::vector<double> function(double f[]) {

	std::vector<double> out;

	// out.push_back(1);
	// out.push_back(1);
	// out.push_back(1);

	out.push_back(f[1] - f[2]);
	out.push_back(f[2] - f[0]);
	out.push_back(f[0] - f[1]);

	return out;
}
