#include "Nelder-Mead/nelder_mead.hpp"

int main(int argc, char const *argv[])
{
	int N = 1;
	if (argc > 1) {
		N = std::stoi(argv[1]);
	}
	else {
		return 0;
	}

	int maxIterations = 1000;
	float acc = 0.001f;
	float real_val = -9.1038f;

	int dim = 3 * N;
	float **p = new float *[dim + 1];
	for (size_t t = 0; t < dim + 1; t++) {
		p[t] = new float[dim];
		for (size_t z = 0; z < dim; z++) {
			p[t][z] = RandomGenerator::generateFloat(-2.5f, 2.5f);
		}
	}

	NelderMead obj = NelderMead(dim, p, &ELJ);
	obj.initialize(maxIterations, acc, real_val);
	obj.minimize();

	for (size_t t = 0; t < dim + 1; t++) {
		if (p[t] != nullptr) {
			delete[] p[t];
		}
	}
	delete[] p;
	return 0;
}
