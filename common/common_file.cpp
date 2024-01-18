#include "common_file.hpp"
#include <chrono>

std::default_random_engine RandomGenerator::generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());
float RandomGenerator::generateFloat(float l, float h)
{
	std::uniform_real_distribution<float> dist(l, h);
	return dist(generator);
}

float l2_norm(const float *X, const int i, const int j)
{
	float val = 0;
	for (int z = 0; z < 3; z++) {
		val += (X[i + z] - X[j + z]) * (X[i + z] - X[j + z]);
	}

	return std::sqrt(val);
}

float l2_norm_gen(float *x, float *y, const int N)
{
	float val = 0;
	for (size_t i = 0; i < N; i++) {
		val += (x[i] - y[i]) * (x[i] - y[i]);
	}
	return std::sqrt(val);
}

float ELJ(const int N, const float *X)
{
	float sum = 0;
	for (size_t i = 0; i < N - 3; i += 3) {
		float temp = 0;
		for (size_t j = i + 3; j < N; j += 3) {
			float rij = l2_norm(X, i, j);
			float factor = std::pow((1 / rij), 6);
			temp += factor * factor - factor;
		}
		sum += temp;
	}

	return 4 * sum;
}

/*float ELJ(const int N, const float *X, const float epsilon = 1, const float sigma = 1)
{
	float sum = 0;
	for (size_t i = 0; i < N - 3; i += 3) {
		float temp = 0;
		for (size_t j = i + 3; j < N; j += 3) {
			float rij = l2_norm(X, i, j);
			float factor = std::pow((sigma / rij), 6);
			temp += factor * factor - factor;
		}
		sum += temp;
	}

	return 4 * epsilon * sum;
}*/