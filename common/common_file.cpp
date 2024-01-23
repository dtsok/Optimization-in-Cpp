#include "common_file.hpp"
#include <chrono>

std::default_random_engine RandomGenerator::generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
double RandomGenerator::generateDouble(double l, double h)
{
	std::uniform_real_distribution<double> dist(l, h);
	return dist(generator);
}

int RandomGenerator::generateInt(double l, double h)
{
	std::uniform_int_distribution<int> dist(l, h);
	return dist(generator);
}

double l2_norm(const double *X, const size_t i, const size_t j)
{
	double val = 0;
	for (size_t z = 0; z < 3; z++) {
		val += (X[i + z] - X[j + z]) * (X[i + z] - X[j + z]);
	}

	return std::sqrt(val);
}

double l2_norm_gen(double *x, double *y, const size_t N)
{
	double val = 0;
	for (size_t i = 0; i < N; i++) {
		val += (x[i] - y[i]) * (x[i] - y[i]);
	}
	return std::sqrt(val);
}

double ELJ(const size_t N, const double *X)
{
	double sum = 0;
	for (size_t i = 0; i < N - 3; i += 3) {
		double temp = 0;
		for (size_t j = i + 3; j < N; j += 3) {
			double rij = l2_norm(X, i, j);
			double factor = std::pow((1 / rij), 6);
			temp += factor * factor - factor;
		}
		sum += temp;
	}

	return 4 * sum;
}

/*double ELJ(const int N, const double *X, const double epsilon = 1, const double sigma = 1)
{
	double sum = 0;
	for (size_t i = 0; i < N - 3; i += 3) {
		double temp = 0;
		for (size_t j = i + 3; j < N; j += 3) {
			double rij = l2_norm(X, i, j);
			double factor = std::pow((sigma / rij), 6);
			temp += factor * factor - factor;
		}
		sum += temp;
	}

	return 4 * epsilon * sum;
}*/