#include "GA-Binary/ga_binary.hpp"
#include "GA-Real/ga_real.hpp"
#include "Nelder-Mead/nelder_mead.hpp"
#include "PSO/pso.hpp"

int main(int argc, char const *argv[])
{
	size_t N = 1;
	if (argc > 1) {
		N = std::stoi(argv[1]);
	}
	else {
		return 0;
	}

	size_t maxIterations = 100000 * N;
	double acc = 0.001;
	double real_val = -9.10385; // for N = 5

	size_t dim = 3 * N;

	/* // Nelder-Mead
	double **p = new double *[dim + 1];
	for (size_t t = 0; t < dim + 1; t++) {
		p[t] = new double[dim];
		for (size_t z = 0; z < dim; z++) {
			p[t][z] = RandomGenerator::generateDouble(-2.5, 2.5);
		}
	}

	NelderMead obj = NelderMead(dim, p, &ELJ);
	obj.setParameters(maxIterations, acc, real_val);
	obj.minimize();

	for (size_t t = 0; t < dim + 1; t++) {
		if (p[t] != nullptr) {
			delete[] p[t];
		}
	}
	delete[] p;
	*/

	PSO obj = PSO(dim, &ELJ, -2.5, 2.5);
	obj.setParameters(maxIterations, acc, real_val, dim);
	obj.minimize(true);
	return 0;
}
